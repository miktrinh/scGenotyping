# Run alleleIntegrator on tumour samples
setwd('~/lustre_mt22/CN_methods/')
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
source('scripts/finalScripts/misc.R')

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
minSegLen=1e6
subCl.minSegLen=5e6
skipIfExists = T
normREF = T
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204'

# Import chromInfo
chromInfo_fp = '../chrom_abspos_kb.txt'
chromInfo = read_delim(chromInfo_fp,delim = '\t',col_names = T)

refGenome = '/lustre/scratch119/realdata/mdt1/team78pipelines/reference/Human/GRCH37d5/genome.fa'
refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
gtf = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/gtf10X_GRCh38_120.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)
nParallel=60



#----------- Run AlleleIntegrator on NB dataset (5 samples)
#----------------------------------------------------------- 
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('RCC','NB')){
    # Import seurat object
    #srat = readRDS(file.path('sc_seuratObjects',tumourType,PDID,paste0(PDID,'_processed_v1.RDS')))
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    srat@meta.data$cellID = rownames(srat@meta.data)
    #srat@meta.data$AIcall_normREF_orig = srat@meta.data$AIcall_normREF
    if(normREF){
      output_col = 'AIcall_normREF_v2'
    }else{
      output_col = 'AIcall_v2'
    }
    
    srat@meta.data[output_col] = NA
    
    for(PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      
      message(sprintf('Running AlleleIntegrator for Sample %s - tumourType: %s',PDID,tumourType))
      outDir = file.path(mainDir,'alleleIntegrator_output',tumourType,PDID)
      outDir2 = file.path(mainDir,'../revision_2204_v2/alleleIntegrator_output',tumourType,PDID)
      if(!file.exists(outDir2)){
        print('Making new Dir')
        dir.create(outDir2,recursive = T)
      }
      donorMani = projMani[projMani$PDID == PDID,]
      # Set Sample specific params
      tumourDNA = unique(donorMani$tumourDNA)[!is.na(unique(donorMani$tumourDNA))]
      patientDNA = unique(donorMani$patientDNA)[!is.na(unique(donorMani$patientDNA))]
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
      
      
      ####------------------ Generate Battenberg CN summary file ----------------####
      # Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
      btb.fp = unique(donorMani$battenbergFp[!grepl('^n_',donorMani$SampleID)])
      #----- Processing Battenberg data -------#
      #segs = annotateBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=tgtChrs,removeBalancedSegs=F,longFormat = F,method = 'allelicRatio')  
      
      # For normal adrenal sample PD42184: all CN segments with BAF != 0.5 are clonal --> use major clone config
      if(PDID == "PD42184"){
        segs = processBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,longFormat = F,keepClonalSegs = T,method = 'allelicRatio')  
        segs = segs[segs$clonalType == 'maj']
      }else{
        segs = processBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,longFormat = F,keepClonalSegs = F,method = 'allelicRatio')  
      }
      
      
      #if(PDID=='PD46693'){
      #  segs = segs[segs$idx !=8,]
        #segs[segs$Idx == 7,]$Start = 1
      #}
      #segs = GRanges(segs$Chr,IRanges(segs$Start,segs$Stop),
      #               totCN=segs$tumTot, idx=segs$idx, chr=segs$Chr, 
      #               patNum=segs$patNum, matNum=segs$matNum,
      #               tumFrac = segs$tumFrac,clonalType=segs$type)
      
      #message(sprintf('[PDID %s] final number of segments %d',PDID,nrow(mcols(segs))))
      if(nrow(mcols(segs)) == 0) {next()}
      
      
      
      #if(PDID == 'PD42752-2'){ # NB
      #  segs = segs[seqnames(segs) != 1]
      #}else if(PDID == 'PD47706'){
        #segs = segs[!segs$idx %in% c(9,4)] # remove a subCl seg because this segment should have been 2+2 in majCl, but it's too short and so have been merged into 2+3 config (same as the rest of the chromosome)
      #  segs = segs[!seqnames(segs) %in% c(4,8,13,18)]
      #}else if(PDID == 'PD35918'){
      #  segs = segs[seqnames(segs) != 8]
      #}else if(PDID == 'PD37104'){
      #  segs = segs[seqnames(segs) != 9]
      #}#else if(PDID == 'PD47705'){ # ATRT
      #start(segs) = 1
    #}#else if(PDID == "PD42181"){ # Ewings
      # segs = segs[names(segs) != 7]
    #}
      
      ###### Do we want to skip to the last step? #######
      #if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
      if(skipIfExists & file.exists(file.path('~/lustre_mt22/CN_methods/revision_2204_v2/alleleIntegrator_output',tumourType,PDID,paste0(PDID,'_phCnts.RDS')))){
        #phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))
        phCnts = readRDS(file.path('~/lustre_mt22/CN_methods/revision_2204_v2/alleleIntegrator_output',tumourType,PDID,paste0(PDID,'_phCnts.RDS')))
        
      }else{
        #############################
        # Check genotype consistency
        #Are all the BAMs you're going to use from the same individual?  Check before you start
        if(tumourType %in% c('RCC','NB')){
          outputs = file.path(outDir,'../../../../alleleIntegrator_output',tumourType,PDID,paste0(c(basename(c(patientDNA,tumourDNA)),names(bams10X)),'_genotypeCheck.tsv'))
        }else{
          outputs = file.path(outDir,paste0(c(basename(c(patientDNA,tumourDNA)),names(bams10X)),'_genotypeCheck.tsv'))
        }
        
        if(PDID != "PD37228"){
          genoCheck = matchBAMs(BAMs = c(norm=patientDNA,tum=tumourDNA,bams10X),
                                refGenomes = rep(c(refGenome,refGenome10X),c(2,length(bams10X))),
                                outputs = outputs,
                                liftOvers=rep(c(NA,liftChain),c(2,length(bams10X))),
                                is10X=rep(c(FALSE,TRUE),c(2,length(bams10X))),
                                nParallel=nParallel)
          #If anything is less than 0.8 and you should be concerned...
          message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))  
        }
        
        
        
        ######################
        # Call and phase SNPs
        hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
        hSNPs = findHetSNPs(patientDNA,refGenome,file.path('~/lustre_mt22/CN_methods/revision_2204/alleleIntegrator_output/RCC',paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
        
        #Expectation is that we'll find ~ 3 million of them
        message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
        
        
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
        
        # Save this object so that we can process it faster next time!
        #saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
        saveRDS(phCnts,file.path('~/lustre_mt22/CN_methods/revision_2204_v2/alleleIntegrator_output',tumourType,PDID,paste0(PDID,'_phCnts.RDS')))
      }
      
      
      
      ##############
      #Specify normal reference cells
      if(normREF){
        normIDs=c('Leukocytes')
        normCells = srat@meta.data$cellID[which(srat@meta.data$annot == 'Leukocytes' & srat@meta.data$PDID == PDID)]
        message(sprintf('[PDID %s] %d/%d are specified as Normal Cells',PDID,length(normCells),ncol(srat)))
        normCells = gsub('\\.','_',normCells)
        if(sum(grepl('^4602',normCells))>0){
          normCells = gsub('^4602','',normCells)
        }
        if(sum(normCells %in% phCnts$cellID) == 0){
          warning(sprintf('[PDID %s] Please check cellID, most likely mismatched',PDID))
        }else{
          message('Good normal cells!')
        }
      }else{
        normIDs = NULL
        normCells = c()
        
      }
      
      
      # Filter data
      #You don't **have** to provide cluster information here.  You could just tell filterCells which bacodes represent cells by specificying "passCellIDs".  But clustering information helps a lot with interpretation
      if(all(grepl('^4602',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))){
        passCellIDs = gsub('^4602','',rownames(srat@meta.data[srat@meta.data$PDID==PDID,]))
        clusterIDs = setNames(srat@meta.data$annot[srat@meta.data$PDID==PDID],
                              gsub('4602','',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))
      
      }else if(all(grepl('^CG.SB.',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))){
        passCellIDs = gsub('\\.','_',rownames(srat@meta.data[srat@meta.data$PDID==PDID,]))
        clusterIDs = setNames(srat@meta.data$annot[srat@meta.data$PDID==PDID],
                              gsub('\\.','_',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))
      }else{
        passCellIDs = rownames(srat@meta.data[srat@meta.data$PDID==PDID,])
        clusterIDs = setNames(srat@meta.data$annot[srat@meta.data$PDID==PDID],
                              srat@meta.data$cellID[srat@meta.data$PDID==PDID])
      }
      
      
      #If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
      aggToClust=FALSE
      if(aggToClust & !is.null(clusterIDs)){
        gCnts = aggregateByClusters(phCnts,clusterIDs)
        gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=normIDs)
      }else if (!aggToClust & !is.null(clusterIDs)){
        # Not aggToClust but using clusterInfo, including normCells being Leukocytes
        gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs)
      }else if (!aggToClust & is.null(clusterIDs)){
        # No annotation info is passed
        gCnts = filterCells(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs)
      }
      
      
      ##################
      # Calibrate model
      #Specify the error rate
      gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]
      
      #Detect allele specific expression
      if(!normREF){
        gCnts = calcASE(gCnts)    
      }else{
        if(PDID %in% c('PD42752-1','PD42752-2','PD43255','PD48777',"PD42181",'PD47706',"PD47705")){
          # priorKappa = 60 #Check for other samples)
          gCnts = calcASE(gCnts,normalCells = normCells,priorKappa=60) 
          #}else if(PDID == 'PD48777'){
          #  gCnts = calcASE(gCnts,priorKappa=60)  
        }else{
          gCnts = calcASE(gCnts,normalCells = normCells,expCut=300)    
          
        }
      }
      
      
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
            idx.toRm = c(idx.toRm,maj.seg$idx,min.seg$idx)}
        }else{
          # Remove subclone segment
          idx.toRm = c(idx.toRm,min.seg$idx)  
        }
      }
      
      if(length(idx.toRm) > 0){
        pp = abbSegProb(gCnts,od,segs = gCnts@metadata$segs[!gCnts@metadata$segs$idx %in% idx.toRm],abbFrac = 'tumFrac',globalOnoly=TRUE)    
      }else{
        pp = abbSegProb(gCnts,od,segs = gCnts@metadata$segs,abbFrac = 'tumFrac',globalOnoly=T)  
      }
      
      
      #############
      # Validation
      if(normREF){
        pdf(file.path(mainDir,'../revision_2204_v2/alleleIntegrator_output',paste0(PDID,'_',tumourType,'_rawAI_output.pdf')))  
      }else{
        pdf(file.path(mainDir,'alleleIntegrator_output',paste0(PDID,'_',tumourType,'_rawAI_output_noNormREF.pdf')))
      }
      
      dat = plotRawData(gCnts,returnData=TRUE)
      p = plotPosteriorHeatmap(pp,'nLL')
      print(p)
      dev.off()
      
      if(all(grepl('^4602',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))){
        pp = pp[pp$cellID %in% gsub('^4602','',srat@meta.data$cellID)]  
        m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,gsub('^4602','',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))
      }else if(all(grepl('^CG.SB.',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))){
        pp = pp[gsub('_','.',pp$cellID) %in% gsub('_','.',srat@meta.data$cellID[srat@meta.data$PDID==PDID])]  
        m = match(gsub('_','.',pp[seqnames(pp) == 'genomeWide',]$cellID),gsub('_','.',srat@meta.data$cellID[srat@meta.data$PDID==PDID]))
      }else{
        pp = pp[pp$cellID %in% srat@meta.data$cellID[srat@meta.data$PDID==PDID]]  
        m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,srat@meta.data$cellID[srat@meta.data$PDID==PDID])
      }
      
      sum(is.na(m))
      
      srat@meta.data[[output_col]][srat@meta.data$PDID==PDID][m] = ifelse(pp[seqnames(pp) == 'genomeWide',]$maxPostProb>0.99,pp[seqnames(pp) == 'genomeWide',]$mostLikelyState,'Uncalled')
      srat@meta.data[[output_col]] = ifelse(is.na(srat@meta.data[[output_col]]),'Uncalled',as.character(srat@meta.data[[output_col]]))
      #srat@meta.data$call[m] = ifelse(pp$maxPostProb>0.5,pp$mostLikelyState,'Uncalled')
      
    }
    DimPlot(srat,group.by=output_col)
    #table(srat@meta.data$call,srat@meta.data$seurat_clusters)
    #saveRDS(srat,file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    saveRDS(srat,file.path(mainDir,'../revision_2204_v2/sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    
  }
  
}

  
  
  