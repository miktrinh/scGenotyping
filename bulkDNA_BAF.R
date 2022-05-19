# Plot Tumour bulkDNA BAF #### 
setwd('~/lustre_mt22/CN_methods/')
#############
# Libraries #
#############
library(alleleIntegrator)
library(readxl)
source('~/lustre_mt22/CN_methods/scripts/finalScripts/misc.R')
#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
minSegLen=1e6
subCl.minSegLen=2e7
skipIfExists = T
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204'


refGenome = '/lustre/scratch119/realdata/mdt1/team78pipelines/reference/Human/GRCH37d5/genome.fa'
nParallel=25

#----------- Run AlleleIntegrator on NB dataset (5 samples)
#----------------------------------------------------------- 
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('Ewings')){
    for(PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Running AlleleIntegrator for Sample %s - tumourType: %s',PDID,tumourType))
      # Set output directory
      outDir = file.path(mainDir,'alleleIntegrator_output',tumourType,PDID)
      if(!file.exists(outDir)){
        message(sprintf('[%s]: Cannot find output dir - Please check!...'))
        next()
      }
      
      # Set Sample specific params
      donorMani = projMani[projMani$PDID == PDID,]
      tumourDNA = unique(donorMani$tumourDNA[!is.na(donorMani$tumourDNA)])
      patientDNA = unique(donorMani$patientDNA[!is.na(donorMani$patientDNA)])
      
      ######################
      # Call heterozygous SNPs
      hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
      #Expectation is that we'll find ~ 3 million of them
      message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
      
      baf.out = generateCoverageAndBAF(BAM = tumourDNA,refGenome = refGenome,hSNPs=hSNPs,
                                       outPath = file.path(outDir,paste0(PDID,'_cov_BAF.RDS')),nParallel=nParallel)
    }
  }
}
      