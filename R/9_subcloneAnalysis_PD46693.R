# DGE of majorCl vs subCl in PD46693
setwd('~/lustre_mt22/CN_methods/revision_2204_v2/')

# Library
library(tidyverse)
library(Seurat)
library(DESeq2)

outdir = '~/lustre_mt22/CN_methods/revision_2204_v2/PD46693_clonalAnalysis/'
if(!dir.exists(outdir)){
  dir.create(outdir)
}

# Import PD46693 sratObj
pd46693 = readRDS('alleleIntegrator_output_completed/NB/PD46693_probAbb.rds')
DimPlot(pd46693,group.by = 'cell_type')
FeaturePlot(pd46693,features = 'chr4.probAbberant')


# Subset for Tumour population only ####
tumour = subset(pd46693,subset = cell_type == 'Tumour')
tumour$clonalType = ifelse(tumour$chr4.probAbberant > 0.99,'Minor',ifelse(tumour$chr4.probAbberant <0.01,'Major','uncalled'))
Idents(tumour) = 'clonalType'
majCellIDs = rownames(tumour@meta.data[tumour@meta.data$clonalType == 'Major',])
subCellIDs = rownames(tumour@meta.data[tumour@meta.data$clonalType == 'Minor',])

#######========================================================================================================

# Get expression matrices for major clone and subclone
major.mtx = tumour@assays$RNA@counts[,colnames(tumour@assays$RNA@counts) %in% majCellIDs]
minor.mtx = tumour@assays$RNA@counts[,colnames(tumour@assays$RNA@counts) %in% subCellIDs]
m=match(rownames(major.mtx),rownames(minor.mtx))
sum(is.na(m))

# Pseudo Bulk analysis for all Genes ####
cts = cbind(major.mtx,minor.mtx)
colnames(cts) = c(rep('majCl',ncol(major.mtx)),rep('minCl',ncol(minor.mtx)))
coldata = data.frame(condition = c(rep('majCl',ncol(major.mtx)),rep('minCl',ncol(minor.mtx))))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
# Remove genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# Run DESeq2
dds <- DESeq(dds)
results <- results(dds)
results
plotMA(results)

res.filtered = results %>% as.data.frame() %>% filter(!is.na(padj))
res.filtered$sig = ifelse(res.filtered$padj < 0.05,T,F)
sigGenesAll = rownames(res.filtered)[res.filtered$padj <0.05]


# DEG for Transcription Factor only ####
# Get all TF gene names
TF = read.delim('~/Homo_sapiens_transcription_factors_gene_list.txt',sep = '\t')
gene.df = read.csv('~/Hsap_gene_list.csv')
# Match gene names and EnsID
TF = TF[TF$Ensembl.ID %in% gene.df$ensembl_id,]
m <- match(TF$Ensembl.ID, gene.df$ensembl_id)
sum(is.na(m))
TF$gene_name <- gene.df$gene[m]
TF$chr <- gene.df$chr_name[m]


tf.res = res.filtered[rownames(res.filtered) %in% TF$gene_name,]
tf.res$padj2 = p.adjust(tf.res$pvalue,method='BH')
tf.res$sig = ifelse(tf.res$padj2 < 0.05,TRUE,FALSE)
tf.sigGenes = rownames(tf.res[tf.res$padj2 < 0.05,])


## MA plot ##
pseudoCnt = data.frame(majCnt = rowSums(major.mtx),minCnt = rowSums(minor.mtx)[m])
# Normalize by read depth
pseudoCnt = sweep(pseudoCnt,2,colSums(pseudoCnt),'/')
# Transform the count
pseudoCnt$logMaj = log2(pseudoCnt$majCnt)
pseudoCnt$logMin = log2(pseudoCnt$minCnt)
pseudoCnt$M = pseudoCnt$logMaj - pseudoCnt$logMin
pseudoCnt$A = log((pseudoCnt$majCnt + pseudoCnt$minCnt)/2)

# Match all Significantly DE genes
pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c('NTRK1', 'BCL11A', 'TH', 'CHGB', 'HMX1'),'red',
                       ifelse(rownames(pseudoCnt) %in% c(sigGenesAll,tf.sigGenes),'yes','no'))
pseudoCnt$gene = rownames(pseudoCnt)
write_delim(pseudoCnt,file.path(outdir,'PD46693_dataForMAplot.txt'),delim = '\t',col_names = T)

res.filtered$gene = rownames(res.filtered)
tf.res$gene = rownames(tf.res)
write_delim(res.filtered,file.path(outdir,'PD46693_allGenes_pseudobulkDESeq2_all.txt'),delim = '\t',col_names = T)
write_delim(tf.res[tf.res$sig==T,],file.path(outdir,'PD46693_TF_pseudobulkDESeq2_sigOnly.txt'),delim = '\t',col_names = T)


pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c(sigGenes,tf.sig),'yes','no')
pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c(sigGenes,tf.sig),'yes','no')
pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c('NTRK1', 'BCL11A', 'TH', 'CHGB', 'HMX1'),'red',as.character(pseudoCnt$sig))
sigcol = c(yes = 'black',no='lightgrey',red = 'red')



pb.out$gene = rownames(pb.out)
write.csv(pb.out,'~/CN_method_tmp/figures/GOSH25_MAplot_allGenes_pseudobulkDESeq2.csv')


