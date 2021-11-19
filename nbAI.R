# Run AlleleIntegrator on Neuroblastoma samples

#############
# Libraries #
#############
# Load libraries
library(alleleIntegrator)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(GenomicFeatures)
library(Seurat)
library(readxl)
library(tidyverse)

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
minSegLen=1e6

refGenome = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/DNA/genomeDNA.fa'
refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
gtf = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/gtf10X_GRCh38_120.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)
nParallel=60

nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')

nb.srat@meta.data$AI_sc_call3 = NA
nb.srat@meta.data$postProb_abbFrac3 = as.numeric(NA)
nb.srat@meta.data$postProb_normFrac3 = as.numeric(NA)

#----------- Run AlleleIntegrator on NB dataset (5 samples)
#----------------------------------------------------------- 
setwd('/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/NB/')
# 1. Preprocessing / QCing GOSH NB data (5 donors) ####
# Import Manifest
projMani = read_excel("../../../projectManifest.xlsx",sheet = "NB_mani")
#View(projMani)
# Import chromInfo
chromInfo_fp = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt'
chromInfo = read_delim(chromInfo_fp,delim = '\t',col_names = T)



PDID = unique(projMani$PDID)[4]
for(PDID in unique(projMani$PDID)){
  message(sprintf('Running AlleleIntegrator for Sample %s',PDID))
  outDir = file.path('./',PDID)
  if(!file.exists(outDir)){
    print('Making new Dir')
    dir.create(outDir)
  }
  
  donorMani = projMani[projMani$PDID == PDID,]
  
  # Set Sample specific params
  tumourDNA = unique(donorMani$tumourDNA)
  patientDNA = unique(donorMani$patientDNA)
  
  if(length(tumourDNA) != 1 || length(patientDNA) != 1){
    message(sprintf('[DonorID: %s ]\nIncorrect file paths: \ntumourDNA: %s\npatientDNA: %',PDID,tumourDNA,patientDNA))
    next()
  }
  
  bams10X = unique(donorMani$bams10X)
  names(bams10X) = gsub('^4602','',donorMani$sangerSampleID)
  
  for(i in seq_along(bams10X)){
    if(!grepl(names(bams10X)[i],bams10X[i])){
      message(sprintf('[DonorID: %s ] Incorrect scRNAseq BAM file, please check: %s, %s',PDID,names(bams10X)[i],bams10X[i]))
    }
  }
  
  #Get annotated cells
  srat = subset(nb.srat, subset = orig.ident %in% names(bams10X))
  normCells = srat@meta.data$cellID[which(srat@meta.data$cell_type == 'Leukocytes')]
  message(sprintf('[PDID %s] %d/%d are specified as Normal Cells',PDID,length(normCells),ncol(srat)))
  
  

  ####------------------ Generate Battenberg CN summary file ----------------####
  # Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
  btb.fp = unique(donorMani$battenbergFp)
  #----- Processing Battenberg data -------#
  segs = annotateBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = 2e7,PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,longFormat = F,method = 'allelicRatio')  
  
  segs = GRanges(segs$Chr,IRanges(segs$Start,segs$Stop),
                     totCN=segs$tumTot, idx=segs$Idx, chr=segs$Chr, 
                     patNum=segs$patNum, matNum=segs$matNum,
                     tumFrac = segs$tumFrac,clonalType=segs$type)
  names(segs) = seqnames(segs)
  message(sprintf('[PDID %s] final number of segments %d',PDID,nrow(mcols(segs))))
  if(nrow(mcols(segs)) == 0) {next()}
  
  
  #############################
  # Check genotype consistency
  #Are all the BAMs you're going to use from the same individual?  Check before you start
  genoCheck = matchBAMs(BAMs = c(norm=patientDNA,tum=tumourDNA,bams10X),
                        refGenomes = rep(c(refGenome,refGenome10X),c(2,length(bams10X))),
                        outputs = file.path(outDir,paste0(c(basename(c(patientDNA,tumourDNA)),names(bams10X)),'_genotypeCheck.tsv')),
                        liftOvers=rep(c(NA,liftChain),c(2,length(bams10X))),
                        is10X=rep(c(FALSE,TRUE),c(2,length(bams10X))),
                        nParallel=nParallel)
  #If anything is less than 0.8 and you should be concerned...
  message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))
  
  
  ######################
  # Call and phase SNPs
  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  
  #Use tumour DNA to phase them.  Set up plot area so we can inspect the CN changes
  if(PDID == 'PD42752-2'){
    segs = segs[names(segs) != 1]
  }
  par(mfrow=c(1,1))
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tumourDNA,outPath=file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.tsv')),nParallel=nParallel,plotMixtures=F)
  #Liftover to GRCh38
  phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
  #Annotate SNPs using GTF
  phSNPs38 = annotateSNPs(phSNPs38,gtf)
  
  ########################
  # Integrate with 10X.  
  #If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...

  phCnts = getAllelicExpression(loci=phSNPs38,refGenome = refGenome10X,bams = bams10X,
                                outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                                nParallel=nParallel)
  
  ##############
  # Filter data
  #You don't **have** to provide cluster information here.  You could just tell filterCells which bacodes represent cells by specificying "passCellIDs".  But clustering information helps a lot with interpretation
  clusterIDs = setNames(srat@meta.data$cell_type,srat@meta.data$cellID)
  #If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
  aggToClust=FALSE
  if(aggToClust){
    gCnts = aggregateByClusters(phCnts,clusterIDs)
    gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=c('Leukocytes'))
  }else{
    gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=c('Leukocytes'))
  }
  
  
  ##################
  # Calibrate model
  #Specify the error rate
  gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]
  
  #Detect allele specific expression
  if(PDID %in% c('PD42752-1','PD42752-2','PD43255')){
    # priorKappa = 60 #Check for other samples)
    gCnts = calcASE(gCnts,normalCells = normCells,priorKappa=60) 
  }else{
    gCnts = calcASE(gCnts,normalCells = normCells)  
  }
  
  #Get over-dispersion
  od = calcOverDispersion(gCnts)
  
  
  ############
  # Inference
  #dropping sub-clone segments when calculating GenomeWide probs
  subCl.segs = unique(segs[segs$clonalType == 'sub',]$chr)
  pp = abbSegProb(gCnts,od,segs = gCnts@metadata$segs[!gCnts@metadata$segs$chr %in% subCl.segs],abbFrac = 'tumFrac',globalOnoly=TRUE)  
  
  #############
  # Validation
  dat = plotRawData(gCnts,returnData=TRUE)
  plotPosteriorHeatmap(pp,'nLL')
  
  #Integrate global call with Seurat, dropping chr4 with sub-clone
  m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,srat@meta.data$cellID)
  sum(is.na(m))
  if(sum(is.na(m)) > 0){
    stop(sprintf('Sample %s - Cell_ID mismatches detected',PDID))
  }
  srat@meta.data$call = NA
  srat@meta.data$call[m] = ifelse(pp[seqnames(pp) == 'genomeWide',]$maxPostProb>0.99,pp[seqnames(pp) == 'genomeWide',]$mostLikelyState,'Uncalled')
  DimPlot(srat,group.by='call')
  table(srat@meta.data$call,srat@meta.data$finalAnn)
  table(nb.srat@meta.data$AI_sc_call2,nb.srat@meta.data$finalAnn,nb.srat@meta.data$PDID)
  # Add to the big all sample NB sratObject
  m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,nb.srat@meta.data$cellID)
  if(sum(is.na(m)) > 0){
    stop(sprintf('%s failed at final step...',PDID))
  }
  
  nb.srat@meta.data$AI_sc_call3[m] = ifelse(pp[seqnames(pp) == 'genomeWide',]$maxPostProb>0.99,pp[seqnames(pp) == 'genomeWide',]$mostLikelyState,'Uncalled')
  nb.srat@meta.data$postProb_abbFrac3[m] = pp[seqnames(pp) == 'genomeWide',]$postProb_abbFrac
  nb.srat@meta.data$postProb_normFrac3[m] = pp[seqnames(pp) == 'genomeWide',]$postProb_normFrac
}
  
  ##############
  # For BAF dot plot
  # Add...
  phCnts$altIsMum = ifelse(phCnts$altCountTum > phCnts$refCountTum,TRUE,FALSE)
  phCnts$matCount = ifelse(phCnts$altIsMum,phCnts$altCount,phCnts$refCount)
  phCnts$patCount = ifelse(phCnts$altIsMum,phCnts$refCount,phCnts$altCount)
  phCnts$informative = NA
  
  clusterIDs = setNames(srat@meta.data$cell_type,srat@meta.data$cellID)
  gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=c('Leukocytes'),dropUninformative = F)
  saveRDS(gCnts,file.path('/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/',paste0(PDID,'tmp_gCnts_allhSNPs.RDS')))
  
}


DimPlot(nb.srat,group.by = 'AI_sc_call')
saveRDS(nb.srat,'/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')



####===== Pink/Grey barplot ======####
# Extract AI.NB output
nb.AIdata = nb.srat@meta.data %>% group_by(PD_ID,cell_type,AI_sc_call) %>% summarise(nCells = n())
colnames(nb.AIdata) = c('PD_ID','cell_type','AI_output','nCells')
nb.AIdata$PD_ID = factor(nb.AIdata$PD_ID,levels = c("PD42184","PD42752-1","PD42752-2","PD46693","PD43255"))
nb.AIdata$AI_output = ifelse(nb.AIdata$AI_output == 'abbFrac','Aneuploid',
                             ifelse(nb.AIdata$AI_output == 'normFrac','Diploid','Uncalled'))
#CKdata$CKpred.normREF.default.80perc.2397 = factor(CKdata$CKpred.normREF.default.80perc.2397,levels=c('NA','diploid','aneuploid'))
pdf('~/work/AIplots/oct21/NB_AI_anpdip.pdf',width = 4,height = 4)
ggplot(nb.AIdata,aes(y=nCells,fill=AI_output,x=PD_ID))+
  geom_bar(position="fill", stat="identity",width = 0.7)+
  facet_grid(vars(cell_type)) + 
  scale_fill_manual(values = c('#DE006F','#D3D3D3','grey'))+
  theme_bw(base_size = 5)+ labs(fill="AI output")+ xlab('') + ylab('') +
  scale_y_continuous(n.breaks = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        axis.ticks.x = element_blank())
#axis.text.x = element_text(angle =90,hjust = 1,colour = 'black'),
dev.off()


####--------- GOHS025 --------#####
gosh25 = subset(nb.srat,subset = gosh_ID == 'GOSH25')
m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,gosh25@meta.data$cellID)
if(sum(is.na(m)) > 0){
  stop(sprintf('%s failed at final step...',PDID))
}
#gosh25@meta.data$probs[m] = pp[seqnames(pp) == 'genomeWide',]$maxPostProb

# Add chr4 probs
pp = abbSegProb(gCnts,od,segs = gCnts@metadata$segs,abbFrac = 'tumFrac',globalOnoly=F)  
pp.subCl = pp[names(pp) == '4.1']
m = match(pp.subCl$cellID,gosh25@meta.data$cellID)
sum(is.na(m))
gosh25@meta.data$chr4.probAbberant2[m] = pp.subCl$probAbberant

FeaturePlot(gosh25, features = 'chr4.probAbberant')

saveRDS(gosh25,'/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/NB/GOSH25_probAbb.rds')

# just tumour cells
tum.gosh25 = subset(gosh25,subset = cell_type == 'Tumour')
# Recluster
# Clustering
tum.gosh25 = NormalizeData(tum.gosh25)
tum.gosh25 = FindVariableFeatures(tum.gosh25)
tum.gosh25 = ScaleData(tum.gosh25, features = rownames(tum.gosh25))
tum.gosh25 = RunPCA(tum.gosh25, npcs = 75)
ElbowPlot(tum.gosh25, ndims = 75)
tum.gosh25 = FindNeighbors(tum.gosh25, dims=1:40)
tum.gosh25 = FindClusters(tum.gosh25,resolution = 1)
tum.gosh25 = RunUMAP(tum.gosh25, dims=1:40)
tum.gosh25$AI_sc_call
sum(is.na(tum.gosh25$chr4.probAbberant))
tum.gosh25$chr4.probAbberant = as.numeric(tum.gosh25$chr4.probAbberant)
pdf('~/work/AIplots/oct21/gosh25_Tumour_umap.pdf',width = 4,height = 4)
FeaturePlot(tum.gosh25,features = 'chr4.probAbberant',cols = brewer.pal(5,'RdBu')[4:1],pt.size = 0.4)+
  ggtitle('GOSH25 - Tumour Subclone')+
  theme(text = element_text(size=10))
dev.off()
DimPlot(tum.gosh25,group.by = 'AI_sc_call')



