setwd('/lustre/scratch117/casm/team274/mt22/CN_methods/caveman_analysis/')


#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(VariantAnnotation)
library(alleleIntegrator)


# Set global parameters
refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
mainDir = '/lustre/scratch117/casm/team274/mt22'
outDir = './'
nParallel=24


###------- Extract SNVs coverage per cell per sample from CAVEMAN output ------###
###----------------------------------------------------------------------------###
# Run caveman analysis for all samples (NB + RCC)
all.out = data.frame()
for(sheet in c('NB_mani','RCC_mani')){
  # Import projMani
  projMani = read_excel(file.path(mainDir,"projectManifest.xlsx"),sheet = sheet)
  if(sheet == 'NB_mani'){
    srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')  
  }else if (sheet == 'RCC_mani'){
    srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')  
  }
  
  
  ###1. Caveman output
  #PDID=unique(projMani$PDID)[1]
  for(PDID in unique(projMani$PDID)){
    print(PDID)
    
    donorMani = projMani[projMani$PDID == PDID,]
    cm.vcf = unique(donorMani$caveman_matched)[!is.na(unique(donorMani$caveman_matched))]
    # Check if file exist
    if(!file.exists(cm.vcf)){
      warning(sprintf('[%s]  Unable to detect caveman file: %s', PDID, cm.vcf))
      next
    }
    # Read caveman output VCF file
    vcf = readVcf(cm.vcf)
    # Keep only SNVs that PASSed all Filters
    vcf.pass = vcf[vcf@fixed$FILTER == 'PASS',]
    snvs = rowRanges(vcf.pass)
    message(sprintf('[%s] Number of SNVs detected from Caveman output is %d',PDID,nrow(mcols(snvs))))
    # convert ALT from DNAStringSetList to Character
    snvs$ALT = unlist(unstrsplit(CharacterList(snvs$ALT), sep = ','))
    #Filter out anything that has more than one ALT (only keep substitution SNVs)
    snvs = snvs[lengths(snvs$ALT)==1,] # which is all of them!
    
    
    
    ### 2. Using AlleleCounter to look for these SNVs in scRNAseq data
    # Lift SNVs coordinate from refDNA(hg37) over to ref10X (hg38)
    #Liftover to GRCh38
    snvs38 = changeGenomeVersion(snvs,liftChain)
    # Get scRNAseq BAMS
    bams10X = unique(donorMani$bams10X)
    names(bams10X) = gsub('^4602','',donorMani$sangerSampleID)
    
    for(i in seq_along(bams10X)){
      if(!grepl(names(bams10X)[i],bams10X[i])){
        message(sprintf('[DonorID: %s ] Incorrect scRNAseq BAM file, please check: %s, %s',PDID,names(bams10X)[i],bams10X[i]))
      }
    }
    
    
    out = runAlleleCounter(loci=snvs38,refGenome = refGenome10X,bams = bams10X,
                           outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                           nParallel=nParallel)
    
    # Keep SNPs with >=1 read coverage only
    out.sub = out[out$Tot >=1,]
    # Aggregate: for each cell, number of altSNPs with coverage >=5
    sum.out = mcols(out.sub) %>% as.data.frame() %>% group_by(cellID) %>% summarise(nSNPs = n(), totalReads=sum(Tot))
    if(sheet == 'RCC_mani'){
      tmp = data.frame(cellID = srat@meta.data[srat@meta.data$PDID == PDID,]$cellID,
                       finalAnn = srat@meta.data[srat@meta.data$PDID == PDID,]$finalAnn2,
                       PDID=PDID,totalSNPs=nrow(mcols(snvs38)))  
    }else if(sheet == 'NB_mani'){
      tmp = data.frame(cellID = srat@meta.data[srat@meta.data$PDID == PDID,]$cellID,
                       finalAnn = srat@meta.data[srat@meta.data$PDID == PDID,]$cell_type,
                       PDID=PDID,totalSNPs=nrow(mcols(snvs38)))
    }
    
    tmp2 = merge(tmp,sum.out,by='cellID',all.x=T)
    
    
    all.out = rbind(all.out,tmp2)
  }
  
}

all.out$totalReads = ifelse(is.na(all.out$totalReads),0,as.numeric(all.out$totalReads))

all.out$tumourType = ifelse(all.out$PDID %in% c("PD35918","PD37228","PD37104","PD36793"),'RCC','NB')



ggplot(all.out,aes(x=PDID,y=(nSNPs)))+
  geom_jitter(size=0.1,height = 0)+
  #geom_violin()+
  theme_bw() + ylim(c(0,10)) +scale_y_continuous(n.breaks = 10)+
  ylab('# SNVs with >=1 count')


write.csv(all.out,'/lustre/scratch117/casm/team274/mt22/CN_methods/snvCov_allout_2.csv')  





###------- Extract hSNPs coverage per cell per sample from AI output ------###
###------------------------------------------------------------------------###
all.out = data.frame()
for(sheet in c('NB_mani','RCC_mani')){
  # Import projMani
  projMani = read_excel(file.path(mainDir,"projectManifest.xlsx"),sheet = sheet)
  
  if(sheet == 'NB_mani'){
    srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')  
    setwd(file.path(mainDir,'CN_methods/alleleIntegrator_output/NB'))
  }else if (sheet == 'RCC_mani'){
    srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')  
    setwd(file.path(mainDir,'CN_methods/alleleIntegrator_output/RCC'))
  }
  
  #PDID=unique(projMani$PDID)[5]
  for(PDID in unique(projMani$PDID)){
    print(PDID)
    
    donorMani = projMani[projMani$PDID == PDID,]
    # Set Sample specific params
    patientDNA = unique(donorMani$patientDNA)
    ######################
    # Call and phase SNPs
    hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,PDID,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
    #Expectation is that we'll find ~ 3 million of them
    message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
    
    
    
    ### 2. Using AlleleCounter to look for these SNVs in scRNAseq data
    # Lift SNVs coordinate from refDNA(hg37) over to ref10X (hg38)
    #Liftover to GRCh38
    hSNPs38 = changeGenomeVersion(hSNPs,liftChain)
    # Get scRNAseq BAMS
    bams10X = unique(donorMani$bams10X)
    names(bams10X) = gsub('^4602','',donorMani$sangerSampleID)
    
    for(i in seq_along(bams10X)){
      if(!grepl(names(bams10X)[i],bams10X[i])){
        message(sprintf('[DonorID: %s ] Incorrect scRNAseq BAM file, please check: %s, %s',PDID,names(bams10X)[i],bams10X[i]))
      }
    }
    
    
    out = runAlleleCounter(loci=hSNPs38,refGenome = refGenome10X,bams = bams10X,
                           outputs=file.path(mainDir,'CN_methods/hSNPs_coverage2/',paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                           nParallel=nParallel)
    
    # Keep SNPs with >=5 read coverage only
    out.sub = out[out$Tot >=1,]
    # Aggregate: for each cell, number of altSNPs with coverage >=5
    sum.out = mcols(out.sub) %>% as.data.frame() %>% group_by(cellID) %>% summarise(nSNPs = n(), totalReads=sum(Tot))
    if(sheet == 'NB_mani'){
      srat@meta.data$PDID = srat@meta.data$PD_ID
    }
    tmp = data.frame(cellID = srat@meta.data[srat@meta.data$PDID == PDID,]$cellID,PDID=PDID,totalSNPs=nrow(mcols(hSNPs38)))
    tmp2 = merge(tmp,sum.out,by='cellID',all.x=T)
    
    all.out = rbind(all.out,tmp2)
  }
  
}

all.out$totalReads = ifelse(is.na(all.out$totalReads),0,as.numeric(all.out$totalReads))
all.out$tumourType = ifelse(all.out$PDID %in% c("PD35918","PD37228","PD37104","PD36793"),'RCC','NB')
write_csv(all.out,file = '../hSNPcov_allout_2.csv',col_names = T)


## Plot ####
all.out = read.csv('/lustre/scratch117/casm/team274/mt22/CN_methods/hSNPcov_allout.csv')
ggplot(all.out,aes(x=PDID,y=(nSNPs)))+
  geom_jitter(size=0.1,height = 0)+
  #geom_violin()+
  theme_bw() + ylim(c(0,10)) +scale_y_continuous(n.breaks = 10)+
  ylab('# hSNPs with >=1 count')








# Run AlleleCounter
runAlleleCounter = function (loci, refGenome, bams, outputs = NULL,labels = names(bams), 
                             minCounts = 0,verbose = TRUE, nParallel = 1){
  if (is.null(labels) || any(duplicated(labels))) {
    warning("No valid BAM labels found.  Setting to generic Sample1, Sample2, Sample3, etc.")
    labels = paste0("Sample", seq_along(bams))
  }
  
  params = list(bams = bams, refGenome = refGenome, tgtLoci = loci, 
                outputs = outputs, nParallel = nParallel)
  
  cnts = do.call(alleleCounter, params)
  out = list()
  for (i in seq_along(bams)) {
    lab = labels[i]
    dat = cnts[[i]]
    
    dat$cellID = paste0(lab, "_", dat$barcode)
    
    dat$scSource = bams[i]
    dat$sample = lab
    out[[i]] = dat
  }
  out = do.call(c, out)
  bases = c("A", "C", "G", "T")
  vars = c(bases, "Tot", "altCount", "refCount")
  tmp = as.matrix(mcols(out)[, bases])
  out$altCount = tmp[cbind(seq(length(out)), match(out$ALT, bases))]
  out$refCount = tmp[cbind(seq(length(out)), match(out$REF, bases))]
  message('DONE!')
  
  return(out)
}

