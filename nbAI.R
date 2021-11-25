# Run AlleleIntegrator on Neuroblastoma samples
setwd('/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/NB/')

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
source('../../scripts/finalScripts/misc.R')


#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
minSegLen=1e6
subCl.minSegLen=2e7

# Import Manifest
projMani = read_excel("../../../projectManifest.xlsx",sheet = "NB_mani")
# Import chromInfo
chromInfo_fp = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt'
chromInfo = read_delim(chromInfo_fp,delim = '\t',col_names = T)

refGenome = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/DNA/genomeDNA.fa'
refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
gtf = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/gtf10X_GRCh38_120.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)
nParallel=60

nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')

nb.srat@meta.data$AI_sc_call = NA
nb.srat@meta.data$postProb_abbFrac = as.numeric(NA)
nb.srat@meta.data$postProb_normFrac = as.numeric(NA)

#----------- Run AlleleIntegrator on NB dataset (5 samples)
#----------------------------------------------------------- 
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
  segs = annotateBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,longFormat = F,method = 'allelicRatio')  
  
  if(PDID=='PD46693'){
    segs = segs[segs$Idx !=6,]
    segs[segs$Idx == 7,]$Start = 1
  }
  segs = GRanges(segs$Chr,IRanges(segs$Start,segs$Stop),
                       totCN=segs$tumTot, idx=segs$idx, chr=segs$Chr, 
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
  if(PDID == 'PD46693'){
    pp.subCl = abbSegProb(gCnts,od,segs = gCnts@metadata$segs,abbFrac = 'tumFrac',globalOnoly=TRUE)  
  }
  
  #dropping sub-clone segments when calculating GenomeWide probs
  subCl.segs = unique(segs[segs$clonalType == 'sub',]$chr)
  idx.toRm = c()
  for(chr in subCl.segs){
    maj.seg = segs[seqnames(segs) == chr & segs$clonalType == 'maj']
    min.seg = segs[seqnames(segs) == chr & segs$clonalType == 'sub']
    if(length(findOverlaps(maj.seg,min.seg))>0){
      if(unique(maj.seg$tumFrac) != unique(min.seg$tumFrac)){
        idx.toRm = c(idx.toRm,maj.seg$idx,min.seg$idx)  
      }
    }else{
      # Remove subclone segment
      idx.toRm = c(idx.toRm,min.seg$idx)  
    }
  }
  
  pp = abbSegProb(gCnts,od,segs = gCnts@metadata$segs[!gCnts@metadata$segs$idx %in% idx.toRm],abbFrac = 'tumFrac',globalOnoly=TRUE)  
  
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
  
  # Add to the big all sample NB sratObject
  m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,nb.srat@meta.data$cellID)
  if(sum(is.na(m)) > 0){
    stop(sprintf('%s failed at final step...',PDID))
  }
  
  nb.srat@meta.data$AI_sc_call[m] = ifelse(pp[seqnames(pp) == 'genomeWide',]$maxPostProb>0.99,pp[seqnames(pp) == 'genomeWide',]$mostLikelyState,'Uncalled')
  nb.srat@meta.data$postProb_abbFrac[m] = pp[seqnames(pp) == 'genomeWide',]$postProb_abbFrac
  nb.srat@meta.data$postProb_normFrac[m] = pp[seqnames(pp) == 'genomeWide',]$postProb_normFrac

  
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



####--------- PD46693 --------#####
pd46693 = subset(nb.srat,subset = PD_ID == 'PD46693')
# Add chr4 probs
m = match(pp.subCl[names(pp.subCl) == '4.1',]$cellID,pd46693@meta.data$cellID)

if(sum(is.na(m)) > 0){
  stop(sprintf('%s - Mismatches in cellIDs',PDID))
}

pd46693@meta.data$chr4.probAbberant[m] = pp.subCl[names(pp.subCl) == '4.1',]$probAbberant

FeaturePlot(pd46693, features = 'chr4.probAbberant')

saveRDS(pd46693,'/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/NB/PD46693_probAbb.rds')



