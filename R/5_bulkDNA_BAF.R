# Plot Tumour bulkDNA BAF #### 
setwd('~/lustre_mt22/CN_methods/')
#############
# Libraries #
#############
library(alleleIntegrator)
library(readxl)
source('~/lustre_mt22/CN_methods/scripts/finalScripts/R/misc.R')
#########################
# Set Global parameters #
#########################
skipIfExists = T
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2'
plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots'


refGenome = '/lustre/scratch119/realdata/mdt1/team78pipelines/reference/Human/GRCH37d5/genome.fa'
nParallel=25

#----------- Run AlleleIntegrator on NB dataset (5 samples)
#----------------------------------------------------------- 
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c("ATRT")){
    for(PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Generating BAF plot from bulkDNA for Sample %s - tumourType: %s',PDID,tumourType))
      # Set output directory
      outDir = file.path(mainDir,'alleleIntegrator_output_completed',tumourType,PDID)
      if(!file.exists(outDir)){
        message(sprintf('[%s]: Cannot find output dir - Please check!...',PDID))
        dir.create(outDir,recursive = T)
        #next()
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
      if(!file.exists(file.path(outDir,paste0(PDID,'_BAF_granges.RDS')))){
        baf.out = generateCoverageAndBAF(BAM = tumourDNA,refGenome = refGenome,hSNPs=hSNPs,
                                         outPath = file.path(outDir,paste0(PDID,'_cov_BAF_july.tsv')),nParallel=nParallel)
        saveRDS(baf.out,file.path(outDir,paste0(PDID,'_BAF_granges.RDS')))  
      }else{
        baf.out = readRDS(file.path(outDir,paste0(PDID,'_BAF_granges.RDS')))
      }
      minCoverage=10
      #Filter to just the ones that we trust
      filt = baf.out[baf.out$coverage>=minCoverage,]
      # Subset randomly 50% of the points
      set.seed(2397)
      idx = sample(1:nrow(mcols(filt)), nrow(mcols(filt))/2, replace=FALSE)
      filt.sub = filt[idx,]
      dd = as.data.frame(mcols(filt.sub))
      dd$hSNP_pos = rownames(dd)
      dd = dd[,c('BAF','coverage', 'isHet','logR','A','C','G','T','Tot','refCnts','altCnts')]
      plotFun = function(noFrame=F,noPlot=FALSE,minCoverage=10){
        #Filter to just the ones that we trust
        filt = baf.out[baf.out$coverage>=minCoverage,]
        #Work out the chromosome boundaries
        chrsToPlot=c(1:22)
        chrs = chrsToPlot
        chrLens = seqlengths(filt)
        tmp = sapply(split(start(filt),as.character(seqnames(filt))),max)
        chrLens[is.na(chrLens)] = tmp[names(chrLens)[is.na(chrLens)]]
        chrLens = as.numeric(chrLens[chrs])
        
        
        x = start(filt) +cumsum(c(0,chrLens))[match(as.character(seqnames(filt)),chrs)]
        
        # Subset randomly 50% of the points
        set.seed(2397)
        idx = sample(1:nrow(mcols(filt)), nrow(mcols(filt))/2, replace=FALSE)
        filt.sub = filt[idx,]
        
        # BAF plot
        par(mfrow=c(1,1),mar=c(2.1,4.1,1.1,1.1))
        alpha = max(0.002,min(1,1e5/length(filt.sub)))
        plot(x[idx],filt.sub$BAF,
             col=rgb(0,0,0,alpha=alpha/2), 
             cex=0.01,
             las=2,
             xaxt='n',
             yaxt='n',
             xlab='',
             ylab='BAF',
             xaxs='i',
             yaxs='i',
             ylim=c(0,1),
             xlim=c(1,sum(chrLens)))
        
        axis(side=2, at=c(0,0.5,1),labels=c(0,0.5,1),las=1)
        abline(v=cumsum(chrLens),col='lightgrey')
        
      }
      
      saveFig(file.path(plotDir,paste0('SupFig_dnaBAF_',PDID,'_',tumourType)),plotFun,width = 5.8,height = 2.2,res=1000,rawData = dd)
    } 
  }
}

