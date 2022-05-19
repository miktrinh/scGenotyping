# DGE of majorCl vs subCl in PD46693

# Load Libraries
library(tidyverse)
library(Seurat)

# Import PD46693 sratObj
pd46693 = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/NB/PD46693_probAbb.rds')
DimPlot(pd46693,group.by = 'cell_type')
FeaturePlot(pd46693,features = 'chr4.probAbberant')


# Subset for Tumour population only ####
tumour = subset(pd46693,subset = cell_type == 'Tumour')

# Clustering
tumour = NormalizeData(tumour)
tumour = FindVariableFeatures(tumour)
tumour = ScaleData(tumour, features = rownames(tumour))
tumour = RunPCA(tumour, npcs = 75)
ElbowPlot(tumour, ndims = 75)
tumour = FindNeighbors(tumour, dims=1:60)
tumour = FindClusters(tumour,resolution = 3)
tumour = RunUMAP(tumour, dims=1:55,n.neighbors = 30,min.dist = 0.6)

sum(is.na(tumour$chr4.probAbberant))
tumour$chr4.probAbberant = as.numeric(tumour$chr4.probAbberant)
FeaturePlot(tumour,features = 'chr4.probAbberant',cols = brewer.pal(5,'RdBu')[4:1],pt.size = 0.4)+ggtitle('PD46693 - Tumour')

table(tumour$clonalType,tumour$clonalType)
tumour$clonalType = ifelse(tumour$chr4.probAbberant > 0.99,'Minor',ifelse(tumour$chr4.probAbberant <0.01,'Major','uncalled'))


# Allele Integrator Dot plots for the major and minor clones ####

# Import Manifest
projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = 'NB_mani')
outdir = '/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/'

PDID = 'PD46693'
donorMani = projMani[projMani$PDID == PDID,]

####------------------ Generate Battenberg CN summary file ----------------####
# Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
# read BTB 
btb.fp = unique(donorMani$battenbergFp)
#----- Processing Battenberg data -------#
dna.data = annotateBTB(btb.fp,subCl.minSegLen = 1e7,PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'allelicRatio')  

# Remove X chromosome 
dna.data = dna.data[dna.data$Chr != 23,]

dna.data$clusterID = 'Tumour'

subCl.dna = dna.data[dna.data$type == 'sub',]
majCl.dna = dna.data[dna.data$type == 'maj',]

####------------------ Get MAF output ----------------####
f = list.files(outdir,pattern = paste0(PDID,"tmp_gCnts_allhSNPs.RDS"))

gCnts = readRDS(file.path(outdir,f))
m = match(gCnts$cellID,clonal@meta.data$cellID)
sum(is.na(m))
gCnts$clonalType = NA
gCnts$clonalType[!is.na(m)] = clonal@meta.data$clonalType[m[!is.na(m)]]

for(i in c('Major','Minor')){
  message(sprintf('Analyzing %s clone', i))
  gCnts.tmp = gCnts[!is.na(gCnts$clonalType)]
  gCnts.tmp = gCnts.tmp[gCnts.tmp$clonalType == i]
  
  ctmp = aggregateByLists(gCnts.tmp, assays = c("matCount", "patCount"), gCnts.tmp$clusterID)  
  
  colnames(ctmp) = gsub("^cellID$", "clusterID", colnames(ctmp))
  ctmp$totCount = ctmp$patCount + ctmp$matCount
  ctmp$MAF = ctmp$matCount/ctmp$totCount
  ctmp$chr = sapply(strsplit(ctmp$regionID,split = ':'),'[',1)
  ctmp$chr = gsub('X',23,ctmp$chr)
  ctmp$chr = as.numeric(as.character(ctmp$chr))
  ctmp$pos = sapply(strsplit(ctmp$regionID,split = ':'),'[',2)
  ctmp$pos = sapply(strsplit(ctmp$pos,split = '_'),'[',1)
  ctmp$pos = as.numeric(ctmp$pos)
  
  # Get absolute genomic position
  ctmp$abspos_kb = ctmp$pos/1000 # if chromosome 1, abspos = pos
  for(r in 1:nrow(ctmp)){
    chrom = as.numeric(ctmp$chr[r])
    if (chrom > 1){
      ctmp$abspos_kb[r] = ctmp$abspos_kb[r] + (chromInfo[(chromInfo$chrom == (chrom-1)) & (chromInfo$arm == 'q'),]$abspos_kb)
    }
  }
  
  ctmp = ctmp[order(c(ctmp$chr,ctmp$abspos_kb),decreasing = F),]
  ctmp$regionID = paste0(ctmp$clusterID,'_',ctmp$regionID)
  
  
  #### Aggregating by read coverage 
  cov = 500
  
  out = ctmp %>% arrange(clusterID,abspos_kb) %>% group_by(clusterID) %>% summarise(totCnt.cumsum = cumsum(totCount),regionID=regionID)
  m=match(out$regionID,ctmp$regionID)
  sum(is.na(m))
  out = cbind(out[,-c(1,3)],ctmp[m,])
  
  out$readCovBin = floor(out$totCnt.cumsum / cov) + 1
  out2 = out %>% arrange(clusterID,readCovBin,abspos_kb) %>% 
    group_by(clusterID,readCovBin) %>% 
    mutate(patCnt.cumsum = cumsum(patCount),matCnt.cumsum = cumsum(matCount)) %>% 
    filter(totCnt.cumsum == max(totCnt.cumsum))
  
  startPos = out %>% arrange(clusterID,readCovBin,abspos_kb) %>% 
    group_by(clusterID,readCovBin) %>% 
    filter(totCnt.cumsum == min(totCnt.cumsum))
  startPos$start = startPos$abspos_kb
  
  out2 = merge(out2,startPos[,c(2,11,12)],by=c('clusterID','readCovBin'))
  
  out2$midPos = (out2$abspos_kb + out2$start)/2
  out2$MAF.readCovBin = out2$matCnt.cumsum / (out2$matCnt.cumsum+out2$patCnt.cumsum)
  
  out2 = arrange(out2,clusterID,readCovBin)
  
  #----- Plotting AlleleIntegrator MAF results! -------#
  # Remove X chromosome 
  dna.data = dna.data[dna.data$Chr != 23,]
  if(i == 'Major'){
    dna.data[dna.data$Chr == 4,]$patNum = 1
    dna.data[dna.data$Chr == 4,]$tumFrac = 1/2
    dna.data[dna.data$Chr == 4,]$config = '2:1'
  }else if(i == 'Minor'){
    dna.data[dna.data$Chr == 4,]$patNum = 0
    dna.data[dna.data$Chr == 4,]$tumFrac = 1
    dna.data[dna.data$Chr == 4,]$config = '2:0'
  }
  
  dna.data$clusterID = 'Tumour'
  
  
  
  out2 = out2[out2$chr != 23,]
  chromInfo2 = chromInfo[chromInfo$chrom != 23,]
  
  plotFun = function(noFrame=T,noPlot=FALSE){
    celltype='Tumour'
    
    par(mar=c(0.5,0.8,1.5,0.1),xpd=TRUE)
    tmp = out2[out2$clusterID == celltype,]
    
    ncells_total = nrow(nb.srat@meta.data[nb.srat@meta.data$PD_ID == PDID & nb.srat@meta.data$finalAnn == celltype,])  
    ncells = length(unique(gCnts.tmp$cellID))
    
    
    # Plot main frame
    plot(out2$midPos*1000, out2$MAF.readCovBin,
         las=1,
         type='n',xaxt='n',yaxt='n',
         main = PDID, cex.main = 0.8,
         #xlim=c(-15.5,15),
         ylim=c(-0.1,1.3),
         frame.plot=F)
    
    text(x=1e6,y=1.32,paste0(i,' clone (',ncells,'/',ncells_total,')'),cex=0.7,family = 'Helvetica',font=1,adj = 0)
    #text(x=2.86e9,y=1.32,paste0('n=',ncells),cex=1,family = 'Helvetica',font=1,adj = 1)
    axis(2,at=c(0,1/2,1),labels = c(0,'1/2',1),las=1,pos = 0,tck = -.02,lwd = 0.3,cex.axis=0.7,hadj = 0.3,padj = 0.5)
    
    #Plot background chromosome
    xleft = c(0,chromInfo2[chromInfo2$arm == 'q' & chromInfo2$chrom!=22,]$abspos_kb*1000)
    xright = c(chromInfo2[chromInfo2$arm == 'q',]$abspos_kb*1000)
    
    col = replicate(c('white','lightgrey'),n = 22/2)
    rect(xleft=xleft,
         xright=xright,
         ybottom=-0.1,
         ytop=1.18,
         col = col,
         lty = 'blank')
    #Black surrounding border
    rect(xleft=min(xleft),
         xright=max(xright),
         ybottom=-0.1,
         ytop=1.18,
         col = colAlpha('white',0.0001),
         border = 'black',lwd = 0.4)
    
    #text(x=(xleft+xright)/2,y = 1.095,labels = c(1:22),cex = c(rep(0.7,10),rep(0.62,4),rep(0.53,4),rep(0.31,4)),font = 1)
    
    
    # Plot AI MAF
    points(tmp$midPos*1000,tmp$MAF.readCovBin,
           pch=19,
           cex=0.1,col=colAlpha('black',0.8))
    
    # Plot ground truth
    if(i=='Major'){
      lines(x=majCl.dna$abspos_kb*1000,majCl.dna$tumFrac,col='#b02c46',lwd=1.5)  
    }else if(i=='Minor'){
      lines(x=subCl.dna$abspos_kb*1000,subCl.dna$tumFrac,col='#b02c46',lwd=1.5)  
    }
    
  }
  
  saveFig(file.path(res,paste0('PD46693_',i,'clone_',cov,'_new')),plotFun,width = 6,height = 1.4,res=500)
}





#######========================================================================================================
tumour = subset(pd46693,subset = cell_type == 'Tumour')
tumour$clonalType = ifelse(tumour$chr4.probAbberant > 0.99,'Minor',ifelse(tumour$chr4.probAbberant <0.01,'Major','uncalled'))
Idents(tumour) = 'clonalType'
majCellIDs = rownames(tumour@meta.data[tumour@meta.data$clonalType == 'Major',])
subCellIDs = rownames(tumour@meta.data[tumour@meta.data$clonalType == 'Minor',])

# Get expression matrices for major clone and subclone
major.mtx = tumour@assays$RNA@counts[,colnames(tumour@assays$RNA@counts) %in% majCellIDs]
minor.mtx = tumour@assays$RNA@counts[,colnames(tumour@assays$RNA@counts) %in% subCellIDs]
m=match(rownames(major.mtx),rownames(minor.mtx))
sum(is.na(m))

# Pseudo Bulk analysis for all Genes ####
library(DESeq2)
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
#res.filtered$TF.sig = ifelse(rownames(res.filtered) %in% tf.sigGenes,T,F)

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
write_delim(pseudoCnt,'../../PD46693_dataForMAplot.txt',delim = '\t',col_names = T)

res.filtered$gene = rownames(res.filtered)
tf.res$gene = rownames(tf.res)
write_delim(res.filtered,'../../PD46693_allGenes_pseudobulkDESeq2_all.txt',delim = '\t',col_names = T)
write_delim(tf.res[tf.res$sig==T,],'../../PD46693_TF_pseudobulkDESeq2_sigOnly.txt',delim = '\t',col_names = T)

write_delim(res.filtered,'~/CN_method_tmp/PD46693_allGenes_pseudobulkDESeq2_all.txt',delim = '\t',col_names = T)
write_delim(tf.res[tf.res$sig==T,],'~/CN_method_tmp/PD46693_TF_pseudobulkDESeq2_sigOnly.txt',delim = '\t',col_names = T)

pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c(sigGenes,tf.sig),'yes','no')
pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c(sigGenes,tf.sig),'yes','no')
pseudoCnt$sig = ifelse(rownames(pseudoCnt) %in% c('NTRK1', 'BCL11A', 'TH', 'CHGB', 'HMX1'),'red',as.character(pseudoCnt$sig))
sigcol = c(yes = 'black',no='lightgrey',red = 'red')


plotFun = function(noFrame=FALSE,noPlot=FALSE){
  par(mfrow=c(1,1),mar=c(1,1,0.2,0.1))
  plot(pseudoCnt$A,pseudoCnt$M,
       las=1,pch=19,cex=0.02,col=sigcol[pseudoCnt$sig],
       #xlim =c(-3.4,-2.6),ylim=c(-0.01,0.015),
       #xlim =c(-3,16),ylim=c(-4.5,5.5),
       #type='n',
       xlab= 'Mean log2 normalized count',ylab='log2 FC',
       #xaxt='n',yaxt='n',
       frame.plot=T)
  #mtext('Mean log2 normalized count',side = 1,line = 0.1,cex = 0.7)
  #mtext('log2 FC',side = 2,line = 0.1,cex = 0.7)
  
  segments(x0=-30,x1 = 0, 
           y0=0, y1=0,
           col = 'black')
  points(pseudoCnt$A,pseudoCnt$M,pch=19,cex=0.02,col=sigcol[pseudoCnt$sig])
  points(pseudoCnt[pseudoCnt$sig == 'red',]$A,pseudoCnt[pseudoCnt$sig == 'red',]$M,pch=19,cex=0.02,col=sigcol[pseudoCnt[pseudoCnt$sig == 'red',]$sig])
  text(pseudoCnt[pseudoCnt$lab == 'yes',]$M~pseudoCnt[pseudoCnt$lab == 'yes',]$A,
       labels=rownames(pseudoCnt[pseudoCnt$lab == 'yes',]),
       data=pseudoCnt, cex=0.4, font=1,pos=1)
}

saveFig(file.path(resd,paste0('GOSH25_MAplot')),plotFun,width = 3.5,height = 2.4,res=500)
pb.out$gene = rownames(pb.out)
write.csv(pb.out,'~/CN_method_tmp/figures/GOSH25_MAplot_allGenes_pseudobulkDESeq2.csv')


