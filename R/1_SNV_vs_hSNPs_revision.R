# This script extract the SNVs and hSNPs coverage for all samples for April 2022 revision

setwd('~/lustre_mt22/CN_methods/')


#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(VariantAnnotation)
library(alleleIntegrator)

#########################
# Set Global parameters #
#########################
# Set global parameters
refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
skipIfExists = T
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204/'


###------- Extract SNVs coverage per cell per sample from CAVEMAN output ------###
###----------------------------------------------------------------------------###
# Run caveman analysis for all samples
all.out = data.frame()

for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('RCC','NB','Wilms','ATRT','Ewings')){
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    # 1. Set outdir
    outDir = file.path(mainDir,'caveman_output')
    
    if(!dir.exists(outDir)){
      print('Making new Dir')
      dir.create(outDir,recursive = T)
    }
    
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      donorMani = projMani[projMani$PDID == current_PDID,]
      # subset annnb.srat object to keep only cells of that sample
      srat.sub = subset(srat, subset = PDID == current_PDID)
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
      message(sprintf('[%s] Number of SNVs detected from Caveman output is %d',current_PDID,nrow(mcols(snvs))))
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
          message(sprintf('[DonorID: %s ] Incorrect scRNAseq BAM file, please check: %s, %s',current_PDID,names(bams10X)[i],bams10X[i]))
        }
      }
      
      
      out = runAlleleCounter(loci=snvs38,refGenome = refGenome10X,bams = bams10X,
                             outputs=file.path(outDir,paste0(current_PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                             nParallel=nParallel)
      
      # Keep SNPs with >=1 read coverage only
      out.sub = out[out$Tot >=1,]
      # Aggregate: for each cell, number of altSNPs with coverage >=5
      sum.out = mcols(out.sub) %>% as.data.frame() %>% group_by(cellID) %>% summarise(nSNPs = n(), totalReads=sum(Tot))
      
      if(all(grepl('^4602',srat.sub@meta.data$cellID))){
        cellID = gsub('^4602','',srat.sub@meta.data$cellID)
      }else if(all(grepl('^CG_SB_',sum.out$cellID))){
        cellID = gsub('\\.','_',srat.sub@meta.data$cellID)
      }else{
        cellID = srat.sub@meta.data$cellID
        message('What is going on?')
      }
      
      
      tmp = data.frame(cellID = cellID,
                       finalAnn = srat.sub@meta.data$annot,
                       PDID=current_PDID,totalSNPs=nrow(mcols(snvs38)))  
      tmp2 = merge(tmp,sum.out,by='cellID',all.x=T)
      tmp2$tumourType = tumourType
      
      all.out = rbind(all.out,tmp2)
    }
  }
}
    
all.out$totalReads = ifelse(is.na(all.out$totalReads),0,as.numeric(all.out$totalReads))
all.out$tumCat = ifelse(all.out$tumourType == 'RCC','adult','pediatric')
write.csv(all.out,'/lustre/scratch117/casm/team274/mt22/CN_methods/revision_2204/snvCov_allout.csv')  


# A few quick plots
ggplot(all.out,aes(x=PDID,y=(nSNPs)))+
  geom_jitter(size=0.1,height = 0)+
  #geom_violin()+
  theme_bw() + ylim(c(0,10)) +scale_y_continuous(n.breaks = 10)+
  ylab('# SNVs with >=1 count')


all.out$cat = ifelse(all.out$totalReads <5,all.out$totalReads,
                     ifelse(all.out$totalReads <=10,'5-10','11+'))
dd=all.out %>% group_by(tumCat,cat) %>% summarise(nCells = n())
dd$cat = factor(dd$cat,levels = c('0','1','2','3','4','5-10','11+'))

ggplot(dd,aes(x=cat,y=nCells,fill=tumCat))+
  geom_col(position = 'dodge')



## Metrics used in manuscript
tmp = snv %>% group_by(PDID) %>% summarise(nCells=n_distinct(cellID),nSNPs = sum(nSNPs,na.rm = T))
sum(all.out$totalReads)/((sum(tmp$nSNPs)/10000) * (sum(tmp$nCells)))

#####################################################################
#                       hSNPs coverage                          #####
#####################################################################

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







setwd('~/lustre_mt22/CN_methods/')


#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(VariantAnnotation)
library(alleleIntegrator)

#########################
# Set Global parameters #
#########################
# Set global parameters
refGenome = '/lustre/scratch119/realdata/mdt1/team78pipelines/reference/Human/GRCH37d5/genome.fa'
refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
skipIfExists = T
nParallel=24
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204/'


###------- Extract SNVs coverage per cell per sample from CAVEMAN output ------###
###----------------------------------------------------------------------------###
# Run caveman analysis for all samples
all.out = data.frame()

for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('RCC','NB','Wilms','ATRT','Ewings')){
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    # 1. Set outdir
    outDir = file.path(mainDir,'hSNP_coverage')
    
    if(!dir.exists(outDir)){
      print('Making new Dir')
      dir.create(outDir,recursive = T)
    }
    
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      donorMani = projMani[projMani$PDID == current_PDID,]
      # subset annnb.srat object to keep only cells of that sample
      srat.sub = subset(srat, subset = PDID == current_PDID)
      
      patientDNA = unique(donorMani$patientDNA)
      
      ######################
      ### 1. Call heterozygous SNPs
      hSNPs = findHetSNPs(patientDNA,refGenome,file.path(mainDir,'alleleIntegrator_output',tumourType,current_PDID,paste0(current_PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
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
      
      # Run alleleCounter through scRNAseq bam files
      out = runAlleleCounter(loci=hSNPs38,refGenome = refGenome10X,bams = bams10X,
                             outputs=file.path(outDir,paste0(current_PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                             nParallel=nParallel)
      
      #print(dim(mcols(out)))
      # Keep SNPs with >=1 read coverage only
      out.sub = out[out$Tot >=1,]
      #print(dim(mcols(out.sub)))
      # Aggregate: for each cell, number of altSNPs with coverage >=5
      sum.out = mcols(out.sub) %>% as.data.frame() %>% group_by(cellID) %>% summarise(nSNPs = n(), totalReads=sum(Tot))
      #sum.out$PDID = PDID
      #sum.out$totalSNPs = nrow(mcols(snvs38))
      if(all(grepl('^4602',srat.sub@meta.data$cellID))){
        cellID = gsub('^4602','',srat.sub@meta.data$cellID)
      }else if(all(grepl('^CG_SB_',sum.out$cellID))){
        cellID = gsub('\\.','_',srat.sub@meta.data$cellID)
      }else{
        cellID = srat.sub@meta.data$cellID
        message('What is going on?')
      }
      
      tmp = data.frame(tumourType=tumourType,cellID = cellID,PDID=current_PDID,totalSNPs=nrow(mcols(hSNPs38)))
      tmp2 = merge(tmp,sum.out,by='cellID',all.x=T)
      
      all.out = rbind(all.out,tmp2)
    }
  }
}
      
      
      
# histogram of how many reads covers any SNV in each cell
all.out$totalReads = ifelse(is.na(all.out$totalReads),0,as.numeric(all.out$totalReads))

write_csv(all.out,file = '/lustre/scratch117/casm/team274/mt22/CN_methods/revision_2204/hSNPcov_allout.csv',col_names = T)


## Plot ####
all.out = read.csv('/lustre/scratch117/casm/team274/mt22/CN_methods/hSNPcov_allout.csv')
ggplot(all.out,aes(x=PDID,y=(nSNPs)))+
  geom_jitter(size=0.1,height = 0)+
  #geom_violin()+
  theme_bw() + ylim(c(0,10)) +scale_y_continuous(n.breaks = 10)+
  ylab('# SNVs with >=1 count')

pdf('~/CN_methods/hSNPs_Hist_wN.pdf',width=3.5,height=3.5)
p=ggplot(all.out,aes(x=log10(totalReads)))+
  #geom_density()+
  geom_histogram(bins = 200,fill='white',col='black',size=0.3)+
  #facet_wrap(vars(PDID),scales = 'free')+
  xlab('Log10 (# Reads covering all SNVs per cell)') + ylab('# Cells') +
  theme_bw(base_size = 8)

print(p)
dev.off()
ylim(c(0,10)) +scale_y_continuous(n.breaks = 10)+
  ylab('# SNVs with >=1 count')




##### Calculate some metrics #######

# hSNPs
all.out = read.csv('/lustre/scratch117/casm/team274/mt22/CN_methods/revision_2204/hSNPcov_allout.csv')
mean(all.out$totalReads)
1e6*mean(all.out$totalReads)/3e9
10/(1e6*mean(all.out$totalReads)/3e9)
# SNVs
all.out = read.csv('/lustre/scratch117/casm/team274/mt22/CN_methods/revision_2204/snvCov_allout.csv')

summary_data = all.out %>% group_by(PDID) %>% summarise(totalReads = sum(totalReads,na.rm = T), totalCells =n_distinct(cellID),nSNVs = unique(totalSNPs))
summary_data$avg_snv_cov = summary_data$totalReads/(summary_data$totalCells*(summary_data$nSNVs/10000))
mean(summary_data$avg_snv_cov)
max(all.out$totalReads)
View(all.out[all.out$totalReads==556,]) # Ewings sample with whole genome duplication etc.

