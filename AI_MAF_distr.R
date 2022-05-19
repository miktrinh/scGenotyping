### BAF generated from alelleIntegrator distribution
### This is different to the BAF plot in supplementary figures because here, BAF are aggregated by genomic bin, not read coverage
### Using alelleIntegrator output when normal reference is provided
setwd('~/lustre_mt22/CN_methods/revision_2204_v2/')

#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(alleleIntegrator)
source('~/lustre_mt22/CN_methods/scripts/finalScripts/misc.R')


#########################
# Set Global parameters #
#########################
binSize=5e6
subCl.minSegLen = 2e7
chrToPlot = c(1:22)
skipIfExists = T
keepClonalSegs = F
# Import Manifest
projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]  
mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'


# Set outdir
outdir = file.path(mainDir,'sensitivity_analysis/AI_BAF_distribution')

if(!dir.exists(outdir)){
  print('Making new Dir')
  dir.create(outdir,recursive = T)
}

# For each tumour type
# loop through each PDID

for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('NB','RCC','ATRT','Wilms','Ewings')){
    #srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    #agg_gCnts_allSamples = tibble()
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Aggregating AI BAF by genomic bin Sample %s - tumourType: %s',current_PDID,tumourType))
      donorMani = projMani[projMani$PDID == current_PDID,]
      # subset annnb.srat object to keep only cells of that sample
      #srat.sub = subset(srat, subset = PDID == current_PDID)
      
      
      ####------------------ Generate Battenberg CN summary file ----------------####
      ### Work out which genomic bins to aggregate across all cell types, which will require aggregation across Tumour cells only
      
      btb.fp = unique(donorMani$battenbergFp)
      #----- Processing Battenberg data -------#
      if(current_PDID == "PD42184"){
        dna.data = processBTB(btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=chrToPlot,removeBalancedSegs=F,longFormat = F,keepClonalSegs = F,method = 'allelicRatio')  
      }else{
        dna.data = processBTB(btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=chrToPlot,removeBalancedSegs=F,longFormat = F,keepClonalSegs = F,method = 'allelicRatio')  
      }
      
      #dna.data = annotateBTB(btb.fp = btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = F,method = 'totalCN')  
      #dna.data = GRanges(dna.data$Chr,IRanges(dna.data$Start,dna.data$Stop),
      #                   matNum=dna.data$matNum,patNum=dna.data$patNum, 
      #                   totCN=dna.data$tumTot,tot2min=dna.data$tot2min, idx=dna.data$idx, Chr=dna.data$Chr, 
      #                   tumFrac = dna.data$tumFrac,clonalType=dna.data$type, posID=dna.data$posID)
      if(keepClonalSegs){
        all_dipSegs = GRanges()
        all_altSegs = GRanges()
        for(chr in chrToPlot){
          if(!chr %in% dna.data$chr){
            next
          }
          # Extract the relevant altered copy number segments - GROUND TRUTH
          
          # If there is subclonal segments within altSegs still
          segs = dna.data[(dna.data$chr == chr),]
          segs$toRemove = F
          if('sub' %in% segs$clonalType){
            print('Found sub!')
            # Find the major clone segment overallping with the subclonal segment
            for(k in 1:(sum(segs$clonalType == 'sub'))){
              seg.sub = segs[segs$clonalType == 'sub'][k]
              o = subsetByOverlaps(segs[segs$clonalType != 'sub' & segs$toRemove == F],seg.sub)
              segs$toRemove[segs$idx %in% o$idx] = T
              # merge and split this into 2 to remove the clonal segments
              if(length(o)>1){
                stop('Weird...need to check!')
              }else if(length(o) == 1){
                print('Here.........s')
                # Check if they are the same segment coordiantes
                if(start(seg.sub) < start(o)){
                  allSeg_tmp1 = o
                  end(allSeg_tmp1) = start(seg.sub) - 1 
                  allSeg_tmp1$toRemove = F
                  segs = append(segs,c(allSeg_tmp1))
                }
                
                if(end(seg.sub) < end(o)){
                  allSeg_tmp2 = o
                  start(allSeg_tmp2) = end(seg.sub) + 1
                  allSeg_tmp2$toRemove = F
                  segs = append(segs,c(allSeg_tmp2))
                }
              }else{
                message('Nothing!')
              }
            }
          }
          
          segs = segs[segs$clonalType != 'sub' & segs$toRemove == F]
          mcols(segs)$toRemove = NULL
          
          # "diploid segments" = segments with an expected BAF of 0.5
          # "altered segments" = segments with an expected BAF of NOT 0.5
          dipSegs = segs[(segs$chr) & segs$tumFrac == 0.5,]
          altSegs = segs[(segs$chr == chr) & (segs$tumFrac != 0.5)]
          
          if((length(altSegs) == 0) & (length(dipSegs) > 1)){ 
            # Merge multiple diploid segments to just 
            maxStop = max(end(dipSegs))
            dipSegs = dipSegs[start(dipSegs) == 1,]
            end(dipSegs) = maxStop
          }else if((length(altSegs) > 0) & (length(dipSegs) > 1)){
            print('Here ...')
            altSegs = altSegs[order(start(altSegs), decreasing = F),]
            dipSegs_new = GRanges()
            for(k in 1:(length(altSegs)+1)){
              current_altStop = ifelse(k<=length(altSegs),start(altSegs)[k],Inf)
              prev_altStop = ifelse(k>1,end(altSegs)[k-1],0)
              
              dipSegs_tmp = dipSegs[(start(dipSegs) > prev_altStop) & (end(dipSegs) < current_altStop)]
              if(length(dipSegs_tmp) == 0){
                next
              }
              maxStop = max(end(dipSegs_tmp))
              dipSegs_tmp = dipSegs_tmp[start(dipSegs_tmp) == min(start(dipSegs_tmp)),]
              end(dipSegs_tmp) = maxStop
              dipSegs_new = append(dipSegs_new,dipSegs_tmp)
            }
            dipSegs = dipSegs_new
          }
          all_dipSegs = append(all_dipSegs,dipSegs)
          all_altSegs = append(all_altSegs,altSegs)
        }
        
        all_dipSegs$segType = 'diploid'
        all_altSegs$segType = 'altered'
        allSegs = append(all_dipSegs,all_altSegs)
      }else{
        allSegs = dna.data
        allSegs$segType = ifelse(allSegs$tumFrac == 0.5,'diploid','altered')
      }
      
      # Remove segments which are shorter than binsize
      allSegs = allSegs[width(allSegs) >= binSize]
      allSegs$Start = start(allSegs)
      allSegs$Stop = end(allSegs)
      
      
      # Identify altered copy number segments withOUT subclones
      allSegs$posID = paste0(allSegs$chr,':',allSegs$Start,'_',allSegs$Stop)
      clonalSegs = unique(names(table(allSegs$posID)[table(allSegs$posID)>1])) # if a segment appears twice, its clonal (major and minor config)
      # There might still be other clonal segments if the subclone segment coordinate lies within the major clone segments (thus posID are different)
      allSegs = allSegs[!allSegs$posID %in% clonalSegs]
      
      
      
      ####------------------ Get MAF output ----------------####
      agg_gCnts_fp = file.path(outdir,paste0(tumourType,'_',current_PDID,'_aggregated_gCnts.txt'))
      if(skipIfExists & file.exists(agg_gCnts_fp)){
        agg_gCnts = read.delim(agg_gCnts_fp,sep = '\t')
      }else{
        gCnts_fp = file.path(mainDir,'alleleIntegrator_output',tumourType,current_PDID,paste0(current_PDID,'_gCnts_allhSNPs.RDS'))
        if(skipIfExists & file.exists(gCnts_fp)){
          gCnts = readRDS(gCnts_fp)
        }else{
          srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
          phCnts = readRDS(file.path(mainDir,'alleleIntegrator_output',tumourType,current_PDID,paste0(current_PDID,'_phCnts.RDS')))
          # Add...
          phCnts$altIsMum = ifelse(phCnts$altCountTum > phCnts$refCountTum,TRUE,FALSE)
          phCnts$matCount = ifelse(phCnts$altIsMum,phCnts$altCount,phCnts$refCount)
          phCnts$patCount = ifelse(phCnts$altIsMum,phCnts$refCount,phCnts$altCount)
          phCnts$informative = NA
          # Filter data
          #You don't **have** to provide cluster information here.  You could just tell filterCells which bacodes represent cells by specificying "passCellIDs".  But clustering information helps a lot with interpretation
          if(all(grepl('^4602',srat@meta.data$cellID[srat@meta.data$PDID==current_PDID]))){
            passCellIDs = gsub('^4602','',rownames(srat@meta.data[srat@meta.data$PDID==current_PDID,]))
            clusterIDs = setNames(srat@meta.data$annot[srat@meta.data$PDID==current_PDID],
                                  gsub('4602','',srat@meta.data$cellID[srat@meta.data$PDID==current_PDID]))
            
          }else if(all(grepl('^CG.SB.',srat@meta.data$cellID[srat@meta.data$PDID==current_PDID]))){
            passCellIDs = gsub('\\.','_',rownames(srat@meta.data[srat@meta.data$PDID==current_PDID,]))
            clusterIDs = setNames(srat@meta.data$annot[srat@meta.data$PDID==current_PDID],
                                  gsub('\\.','_',srat@meta.data$cellID[srat@meta.data$PDID==current_PDID]))
          }else{
            passCellIDs = rownames(srat@meta.data[srat@meta.data$PDID==current_PDID,])
            clusterIDs = setNames(srat@meta.data$annot[srat@meta.data$PDID==current_PDID],
                                  srat@meta.data$cellID[srat@meta.data$PDID==current_PDID])
          }
          gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=c('Leukocytes'),dropUninformative = F)
          saveRDS(gCnts,file.path(mainDir,'alleleIntegrator_output',tumourType,current_PDID,paste0(current_PDID,'_gCnts_allhSNPs.RDS')))  
        }
        
        
        
        # Perform segmentation: generate start and stop positions for each genomic bin
        gCnts$genomic_bin = NA
        gCnts$segType = NA
        gCnts$tumFrac = NA
        gCnts$matNum = NA
        gCnts$patNum = NA
        
        for(i in 1:length(allSegs)){
          seg = allSegs[i,]
          tumFrac = seg$tumFrac 
          if(tumFrac == 0.5 ){
            celltype = 'all'  
          }else{
            celltype = 'Tumour'
          }
          
          
          stops = seq(seg$Start+binSize,seg$Stop,binSize) # if the last segment is shorter than binSize --> will automatically be removed (ie, there wont be a stop position for the last segment. eg seq(2,10,3) gives 2,5,8)
          starts = seq(seg$Start,seg$Stop - binSize,binSize)  
          
          # Add genomic bin information to gCnts object
          for(j in 1:length(stops)){
            seg_start = starts[j]
            seg_stop = stops[j]
            if(length(gCnts[seqnames(gCnts) == seg$chr & start(gCnts) >= seg_start & end(gCnts) < seg_stop])>0){
              gCnts[seqnames(gCnts) == seg$chr & start(gCnts) >= seg_start & end(gCnts) < seg_stop]$genomic_bin = paste0(seg$chr,':',seg_start,'_',seg_stop)
              gCnts[seqnames(gCnts) == seg$chr & start(gCnts) >= seg_start & end(gCnts) < seg_stop]$segType = seg$segType
              gCnts[seqnames(gCnts) == seg$chr & start(gCnts) >= seg_start & end(gCnts) < seg_stop]$tumFrac = seg$tumFrac
              gCnts[seqnames(gCnts) == seg$chr & start(gCnts) >= seg_start & end(gCnts) < seg_stop]$matNum  = seg$matNum 
              gCnts[seqnames(gCnts) == seg$chr & start(gCnts) >= seg_start & end(gCnts) < seg_stop]$patNum  = seg$patNum 
            }
          }
        }
        
        
        gCnts_noNA = gCnts[!is.na(gCnts$genomic_bin)]
        if(length(gCnts_noNA[gCnts_noNA$segType != 'diploid' & gCnts_noNA$clusterID != 'Tumour']$genomic_bin)>0){
          gCnts_noNA[gCnts_noNA$segType != 'diploid' & gCnts_noNA$clusterID != 'Tumour']$genomic_bin = NA  
        }
        
        gCnts_noNA = gCnts_noNA[!is.na(gCnts_noNA$genomic_bin)]
        agg_gCnts = aggregateByLists(gCnts_noNA, assays = c("matCount", "patCount"), regionList = gCnts_noNA$genomic_bin, cellList = gCnts_noNA$segType)  
        
        colnames(agg_gCnts) = c('segType','genomicBin','matCount','patCount')
        agg_gCnts$totCount = agg_gCnts$patCount + agg_gCnts$matCount
        agg_gCnts$MAF = agg_gCnts$matCount/agg_gCnts$totCount
        agg_gCnts$segType = factor(agg_gCnts$segType,levels = c('diploid','altered'))
        agg_gCnts$PDID = current_PDID
        agg_gCnts$tumourType = tumourType
        m = match(agg_gCnts$genomicBin,gCnts_noNA$genomic_bin)
        sum(is.na(m))
        agg_gCnts$matNum = gCnts_noNA$matNum[m]
        agg_gCnts$patNum = gCnts_noNA$patNum[m]
        agg_gCnts$tumFrac = gCnts_noNA$tumFrac[m]
        
        write_delim(agg_gCnts,file.path(outdir,paste0(tumourType,'_',current_PDID,'_aggregated_gCnts.txt')),col_names = T,delim = '\t')
      }
      
      
      #agg_gCnts = read.delim(file.path(outdir,paste0(tumourType,'_',current_PDID,'_aggregated_gCnts.txt')),sep = '\t')
      plotFun = function(noFrame=FALSE,noPlot=FALSE){
        par(mar=c(1.4,1.2,0.5,0.5),xpd=T)
        # plot densities
        xmax = max(density(agg_gCnts$MAF)$x)+0.01
        xmin = min(density(agg_gCnts$MAF)$x)-0.01
        ymax = max(density(agg_gCnts$MAF)$y)
        # plot densities
        d_dipSegs <- density(agg_gCnts$MAF[agg_gCnts$segType == 'diploid' & agg_gCnts$totCount > 0])
        
        if(tumourType=='ATRT'){
          plot(d_dipSegs,main='',axes=F,xlim=c(xmin,xmax),ylim=c(0,121.5))
        }else{
          plot(d_dipSegs,main='',axes=F,xlim=c(xmin,xmax))
        }
        polygon(d_dipSegs, col="white", border="black")
        
        if((sum(agg_gCnts$tumFrac == 1)>2)){
          d_lossSegs <- density(agg_gCnts$MAF[agg_gCnts$tumFrac == 1 & agg_gCnts$totCount > 0])
          lines(d_lossSegs)
          polygon(d_lossSegs, col=adjustcolor("#E2E1E1",alpha.f=0.7), border="black")
        }
        
        if((sum(!agg_gCnts$tumFrac %in% c(1,0.5))>2)){
          d_gainSegs <- density(agg_gCnts$MAF[!agg_gCnts$tumFrac %in% c(1,0.5) & agg_gCnts$totCount > 0])
          lines(d_gainSegs)
          polygon(d_gainSegs, col=adjustcolor("#9C9C9B",alpha.f=0.7), border="black")
        }
        
        #n_cells = n_distinct(gCnts_noNA$cellID)
        
        axis(1,at=c(0,1/3,0.5,2/3,1),labels = c(0,'1/3','0.5','2/3',1), las=1,pos = 0.05,tck=-0.02,cex.axis=0.6,lwd.ticks = 0.8,hadj = 0.5,padj = -2.5,lwd = 0.8)
        axis(2,at=c(0,round(max(d_dipSegs$y))),labels = c('',''),las=1,tck=-0.00,cex.axis=0.7,lwd.ticks = 0.4,hadj = 0.45,padj = 0.5,lwd = 0.8)
        #axis(2,at=c(2000,2800),labels = c('',''),las=1,tck=-0.02,cex.axis=0.6,lwd.ticks = 0,hadj = 0.35,padj = 0.5,lwd = 0.8)
        title(sprintf('MAF distribution in %s - %s',tumourType,current_PDID),family='Helvetica',cex.main = 0.6)
        mtext(side=1,text = 'Maternal Alelle Frequency per 5Mb',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
        mtext(side=2,text = 'Density',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
      }
      saveFig(file.path(outdir,paste0(tumourType,'_',current_PDID,'_MAFdistr')),plotFun,width = 1.8,height = 1.8,rawData = agg_gCnts)    
    }
  }
}
   






### Plot aggregated distribution across all PDIDs and tumourTypes

# For each tumour type
# loop through each PDID
maf_allSamples = tibble()
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('NB','RCC','Wilms','ATRT','Ewings')){
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      
      if(file.exists(file.path(outdir,paste0(tumourType,'_',current_PDID,'_MAFdistr_rawData.tsv')))){
        agg_gCnts = read_delim(file.path(outdir,paste0(tumourType,'_',current_PDID,'_MAFdistr_rawData.tsv')),col_names = T,delim = '\t')
        maf_allSamples = rbind(maf_allSamples,agg_gCnts)  
      }
    }
  }
}


plotFun = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(1.4,1.2,0.5,0.5),xpd=T)
  # plot densities
  xmax = max(density(maf_allSamples$MAF)$x)+0.01
  xmin = min(density(maf_allSamples$MAF)$x)-0.01
  # plot densities
  d_dipSegs <- density(maf_allSamples$MAF[maf_allSamples$segType == 'diploid' & maf_allSamples$totCount > 0])
  plot(d_dipSegs,main='',axes=F,xlim=c(xmin,xmax))
  polygon(d_dipSegs, col="white", border="black")
  
  # plot densities
  #d_allSegs <- density(maf_allSamples$MAF[maf_allSamples$totCount > 0])
  #plot(d_allSegs,main='',axes=F,xlim=c(xmin,xmax))
  #polygon(d_allSegs, col="white", border="black")
  
  d_gainSegs <- density(maf_allSamples$MAF[maf_allSamples$segType != 'diploid' & maf_allSamples$tumFrac != 1 & maf_allSamples$totCount > 0])
  lines(d_gainSegs)
  polygon(d_gainSegs, col=adjustcolor("#9C9C9B",alpha.f=0.6), border="black")
  
  d_lossSegs <- density(maf_allSamples$MAF[maf_allSamples$tumFrac == 1 & maf_allSamples$totCount > 0])
  lines(d_lossSegs)
  polygon(d_lossSegs, col=adjustcolor("#E2E1E1",alpha.f=0.6), border="black")
  #n_cells = n_distinct(gCnts_noNA$cellID)
  
  axis(1,at=c(0,1/3,0.5,2/3,1),labels = c(0,'1/3',0.5,'2/3',1), las=1,pos = 0.05,tck=-0.02,cex.axis=0.6,lwd.ticks = 0.8,hadj = 0.5,padj = -2.5,lwd = 0.8)
  axis(2,at=c(0,round(max(d_dipSegs$y))),labels = c('',''),las=1,tck=-0.00,cex.axis=0.7,lwd.ticks = 0.4,hadj = 0.45,padj = 0.5,lwd = 0.8)
  #axis(2,at=c(2000,2800),labels = c('',''),las=1,tck=-0.02,cex.axis=0.6,lwd.ticks = 0,hadj = 0.35,padj = 0.5,lwd = 0.8)
  title(sprintf('MAF distribution in %s samples from %s tumour types',n_distinct(maf_allSamples$PDID),n_distinct(maf_allSamples$tumourType)),family='Helvetica',cex.main = 0.6)
  mtext(side=1,text = 'Maternal Alelle Frequency per 5Mb',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
  mtext(side=2,text = 'Density',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
}
saveFig(file.path(outdir,'allSamples_AI_MAFdistr'),plotFun,width = 2.4,height = 2.4,rawData = maf_allSamples) 







######################################################
## Setting arbitary threshold of 0.2
thresholds = seq(0,1.01,0.01)
roc_data_by_PDID = tibble()


for(threshold in thresholds){
  maf_allSamples$accuracy = NA
  maf_allSamples$accuracy = ifelse((maf_allSamples$segType == 'diploid') & (abs(maf_allSamples$MAF) <= threshold), 'TN',
                                      ifelse((maf_allSamples$segType == 'diploid') & (maf_allSamples$MAF < -threshold), 'FP',
                                             ifelse((maf_allSamples$segType == 'diploid') & (maf_allSamples$MAF > threshold),'FP',
                                                    
                                                    ifelse((maf_allSamples$segType != 'diploid') & (abs(maf_allSamples$MAF) < threshold),'FN',
                                                           
                                                           ifelse((maf_allSamples$segType != 'diploid') & (abs(maf_allSamples$MAF) >= threshold),'TP','Weird')))))
                                                                  
                                                                  #ifelse((maf_allSamples$segType != 'diploid') & (maf_allSamples$meanScore < -threshold),'FP', # wrong class classification
                                                                         
                                                                   #ifelse((maf_allSamples$segType < 2) & (maf_allSamples$meanScore < -threshold),'TP',
                                                                   #       ifelse((maf_allSamples$segType < 2) & (maf_allSamples$meanScore > threshold),'FP','Weird')))))))) # wrong class classification
  if(sum(maf_allSamples$accuracy == 'Weird')>0){
    stop('Please check weird accuracy assignment...')
  }
  tmp = maf_allSamples %>% group_by(accuracy,PDID) %>% summarise(n=n())
  tmp$threshold = threshold
  
  roc_data_by_PDID = rbind(roc_data_by_PDID,tmp)
}

roc_data_by_PDID$method = 'alleleIntegrator'
write_delim(roc_data_by_PDID,file.path(mainDir,'sensitivity_analysis/ROC/AI_ROC_byPDID_rawdata.txt'),delim = '\t',col_names = T)

roc_data_allSamples = roc_data_by_PDID %>% group_by(accuracy,threshold,method) %>% summarise(total = sum(n))
write_delim(roc_data_allSamples,file.path(mainDir,'sensitivity_analysis/ROC/AI_ROC_rawdata.txt'),delim = '\t',col_names = T)

roc_data_by_PDID = pivot_wider(roc_data_by_PDID,names_from = 'accuracy',values_from = 'n')
roc_data_by_PDID$FP = ifelse(is.na(roc_data_by_PDID$FP),0,roc_data_by_PDID$FP)
roc_data_by_PDID$FN = ifelse(is.na(roc_data_by_PDID$FN),0,roc_data_by_PDID$FN)
roc_data_by_PDID$TP = ifelse(is.na(roc_data_by_PDID$TP),0,roc_data_by_PDID$TP)
roc_data_by_PDID$TN = ifelse(is.na(roc_data_by_PDID$TN),0,roc_data_by_PDID$TN)
# calculate y-axis True Positive Rate (TPR) = TP/(TP+FN)
roc_data_by_PDID$TPR = roc_data_by_PDID$TP/(roc_data_by_PDID$TP + roc_data_by_PDID$FN)
# calculate x-axis False Positive Rate (FPR) = FP/(FP+TN)
roc_data_by_PDID$FPR = roc_data_by_PDID$FP/(roc_data_by_PDID$FP + roc_data_by_PDID$TN)
write_delim(roc_data_by_PDID,file.path(mainDir,'sensitivity_analysis/ROC/AI_ROC_byPDID_dataForPlot.txt'),delim = '\t',col_names = T)











### Read in the ROC data from expression methods
roc_data_by_PDID = read.delim(file.path(mainDir,'sensitivity_analysis/ROC/AI_ROC_byPDID_dataForPlot.txt'),sep = '\t')
roc_data_by_PDID$method = 'alleleIntegrator'

roc_data_by_PDID_expression = read_delim(file.path(mainDir,'sensitivity_analysis/ROC/expressionMethod_ROC_byPDID_dataForPlot_2.txt'),delim = '\t')

roc_by_PDID = rbind(roc_data_by_PDID,roc_data_by_PDID_expression)

##=================== plot ROC for each PDID =================####
m=match(roc_by_PDID$PDID,projMani$PDID)
sum(is.na(m))
roc_by_PDID$tumourType = projMani$TumourType[m]
roc_by_PDID$AUC = NA


lty = c(alleleIntegrator=1,
        CK = 5,
        inferCNV = 3)


for(tumourType in unique(roc_by_PDID$tumourType)){
  for(PDID in unique(roc_by_PDID$PDID[roc_by_PDID$tumourType == tumourType])){
    
    ### Calculate AUC ####
    for(method in unique(roc_by_PDID$method)){
      dat = roc_by_PDID[roc_by_PDID$PDID == PDID & roc_by_PDID$method == method,]
      height = (dat$TPR[-1]+dat$TPR[-nrow(dat)])/2
      width = -diff(dat$FPR)
      AUC = sum(height*width)
      roc_by_PDID$AUC[roc_by_PDID$method==method & roc_by_PDID$PDID == PDID] = AUC
    }
    
    
    data = roc_by_PDID[roc_by_PDID$PDID == PDID,]
    plotFun = function(noFrame=FALSE,noPlot=FALSE){
      par(mar=c(2,2,0.5,0.5),xpd=T)
      line = -8
      
      plot(1, type="n", axes=F,frame.plot = T,xlim=c(0,1),ylim=c(0,1))
      
      for(method in c('alleleIntegrator','CK','inferCNV')){
        
        lines(y= data$TPR[data$method == method],x=data$FPR[data$method == method],lty=lty[method],
              type = "l",las=1,frame.plot = T,axes = F,xlim=c(0,1),ylim=c(0,1),lwd=1.5)  
        
        mtext(paste('AUC =',format(unique(data$AUC[data$method == method]), digits=3),
                    '  Method: ',method),
              side=3, adj=0.9,line=line,cex=0.6,font = 1)
        line = line - 1
      }
      
      axis(side = 2,at =seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,las=1,pos = -0.04,lwd.ticks = 0.5,hadj = 0.1)
      axis(side = 1,at = seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,pos=-0.04,las=0,padj=-2.5)
      
      abline(a=0,b=1,lty=1,col='grey')
      
      mtext(side=1,text = 'FPR',family='Helvetica',font = 1,cex = 0.7,line = 0.9)
      mtext(side=2,text = 'TPR',family='Helvetica',font = 1,cex = 0.7,line = 0.9)
    }
    saveFig(file.path(mainDir,paste0('sensitivity_analysis/ROC/byPDID/',tumourType,'_',PDID,'_inhouseROC')),plotFun,width = 3,height = 3,rawData = data) 
    
    
  }  
}


### Get some metrics for the manuscript
summary_data = roc_by_PDID %>% group_by(tumourType,method) %>% summarise(mean_AUC = mean(AUC,na.rm=T))


#####For all samples ####

roc_data_AI_allSamples = read.delim(file.path(mainDir,'sensitivity_analysis/ROC/AI_ROC_rawdata.txt'),sep='\t')
roc_data_AI_allSamples$method = 'alleleIntegrator'
roc_data_expr_allSamples = read.delim(file.path(mainDir,'sensitivity_analysis/ROC/expressionMethod_ROC_rawdata_2.txt'),sep='\t')
roc_data_allSamples = rbind(roc_data_AI_allSamples,roc_data_expr_allSamples)

roc_data_allSamples = pivot_wider(roc_data_allSamples,names_from = 'accuracy',values_from = 'total')
roc_data_allSamples$FP = ifelse(is.na(roc_data_allSamples$FP),0,roc_data_allSamples$FP)
roc_data_allSamples$FN = ifelse(is.na(roc_data_allSamples$FN),0,roc_data_allSamples$FN)
roc_data_allSamples$TP = ifelse(is.na(roc_data_allSamples$TP),0,roc_data_allSamples$TP)
roc_data_allSamples$TN = ifelse(is.na(roc_data_allSamples$TN),0,roc_data_allSamples$TN)
roc_data_allSamples$TPR = roc_data_allSamples$TP/(roc_data_allSamples$TP + roc_data_allSamples$FN)
# calculate x-axis False Positive Rate (FPR) = FP/(FP+TN)
roc_data_allSamples$FPR = roc_data_allSamples$FP/(roc_data_allSamples$FP + roc_data_allSamples$TN)


#######
roc_data_allSamples$AUC = NA
### Calculate AUC ####
for(method in unique(roc_data_allSamples$method)){
  dat = roc_data_allSamples[roc_data_allSamples$method == method,]
  dat = dat[order(dat$TPR,decreasing = T),]
  height = (dat$TPR[-1]+dat$TPR[-nrow(dat)])/2
  width = -diff(dat$FPR)
  AUC = sum(height*width)
  roc_data_allSamples$AUC[roc_data_allSamples$method==method] = AUC
}


data = roc_data_allSamples
lty = c(alleleIntegrator=1,
        CK = 5,
        inferCNV = 3)

plotFun = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(2,2,0.5,0.5),xpd=T)
  line = -8
  
  plot(1, type="n", axes=F,frame.plot = T,xlim=c(0,1),ylim=c(0,1))
  
  for(method in c('alleleIntegrator','CK','inferCNV')){
    
    lines(y= data$TPR[data$method == method],x=data$FPR[data$method == method],lty=lty[method],
          type = "l",las=1,frame.plot = T,axes = F,xlim=c(0,1),ylim=c(0,1),add=T,lwd=1.5)  
    
    mtext(paste('AUC =',format(unique(data$AUC[data$method == method]), digits=3),
                '  Method: ',method),
          side=3, adj=0.9,line=line,cex=0.6,font = 1)
    line = line - 1
  }
  
  axis(side = 2,at =seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,las=1,pos = -0.04,lwd.ticks = 0.5,hadj = 0.1)
  axis(side = 1,at = seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,pos=-0.04,las=0,padj=-2.5)
  
  abline(a=0,b=1,lty=1,col='grey')
  
  mtext(side=1,text = 'False Positive Rate',family='sans',font = 1,cex = 0.7,line = 0.9)
  mtext(side=2,text = 'True Positive Rate',family='sans',font = 1,cex = 0.7,line = 0.9)
}
saveFig(file.path(mainDir,'sensitivity_analysis/ROC/inhouseROC_allSamples_v2'),plotFun,width = 3,height = 3,rawData = data) 









#### A (slightly) different way to compute ROC and AUC ####
# Data
maf_allSamples$segType = ifelse(maf_allSamples$tumFrac == 0.5,'no_changed','altered')
testres = maf_allSamples$MAF
truestat = maf_allSamples$segType
## OR for expression based methods
all_out_meanScore
for(method in c('CK',"inferCNV")){
  out_meanScore_sub = all_out_meanScore[all_out_meanScore$method == method,]
  
  out_meanScore_sub$segType = ifelse(out_meanScore_sub$nTot == 2,'no_changed','altered')
  out_meanScore_sub$abs_meanScore = abs(out_meanScore_sub$meanScore)
  out_meanScore_sub = out_meanScore_sub[order(out_meanScore_sub$abs_meanScore,decreasing = F),]
  testres = out_meanScore_sub$abs_meanScore  
  truestat = out_meanScore_sub$segType
  
}
# Summary table (Table I in the paper)
( tab=as.matrix(table(truestat, testres)) )
( tot=colSums(tab) )                            # Number of patients w/ each test result
( truepos=unname(rev(cumsum(rev(tab[1,])))) )   # Number of true positives
( falsepos=unname(rev(cumsum(rev(tab[2,])))) )  # Number of false positives
( totpos=sum(tab[1,]) )                         # The total number of positives (one number)
( totneg=sum(tab[2,]) )                         # The total number of negatives (one number)
(sens=truepos/totpos)                           # Sensitivity (fraction true positives)
(omspec=falsepos/totneg)                        # 1 − specificity (false positives)
sens=c(sens,0); omspec=c(omspec,0)


plot(omspec, sens, type="b", xlim=c(0,1), ylim=c(0,1), lwd=1,pch=19,
     xlab="1 − specificity", ylab="Sensitivity") # perhaps with xaxs="i"
grid()
abline(0,1, col="red", lty=2)


height = (sens[-1]+sens[-length(sens)])/2
width = -diff(omspec) # = diff(rev(omspec))
sum(height*width)












### A different way to compute ROC and AUC ####
library(pROC)
pROC_obj = list()
i=1
for(method in c('AI','CK','inferCNV')){
  if(method == 'CK'){
    data = read.delim(file.path(mainDir,paste0('sensitivity_analysis/expressionFC_distribution/allSamples_',method,'_exprDistr_rawData.tsv')))
    data = data[,c("meanScore","celltype")]
    data$celltype = factor(data$celltype,levels = c('all','Tumour'))
  }else if(method == 'inferCNV'){
    data = read.delim(file.path(mainDir,paste0('sensitivity_analysis/expressionFC_distribution/allSamples_',method,'_exprDistr_rawData.tsv')))
    data = data[,c("meanScore","celltype")]
    data$celltype = factor(data$celltype,levels = c('all','Tumour'))
  }else if(method == 'AI'){
    data = read.delim(file.path(mainDir,paste0('sensitivity_analysis/AI_BAF_distribution/allSamples_AI_MAFdistr_v2_rawData.tsv')))
    data = data[,c("MAF","segType")]
    data$segType = factor(data$segType,levels = c('diploid','altered'))
  }
  
  colnames(data) = c('prediction','truth')
  
  pROC_obj[[i]] = roc(data$truth ~ data$prediction,
                  smoothed = T,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=F, auc.polygon=F, max.auc.polygon=F, grid=F,
                  print.auc=F, show.thres=F)
  
  names(pROC_obj)[i] = method
  i = i + 1
}

cols = c("#cb6751",
          "#7aa457",
          "#9e6ebd")

plotFun = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(2,2,0.5,0.5),xpd=T)
  line = -8
  
  plot(1, type="n", axes=F,frame.plot = T,xlim=c(0,1),ylim=c(0,1))
  
  for(i in 1:length(pROC_obj)){
    lines(y= pROC_obj[[i]]$sensitivities,x=1-pROC_obj[[i]]$specificities,col=cols[i],
          type = "l",las=1,frame.plot = T,axes = F,xlim=c(0,1),ylim=c(0,1),add=T,lwd=1.5)  
    
    mtext(paste('AUC =',format((pROC_obj[[i]]$auc), digits=3),
                '  Method: ',names(pROC_obj)[i]),
          side=3, adj=0.9,line=line,cex=0.6,font = 1, col=cols[i])
    line = line - 1
  }
  
  axis(side = 2,at =seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,las=1,pos = -0.04,lwd.ticks = 0.5,hadj = 0.1)
  axis(side = 1,at = seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,pos=-0.04,las=0,padj=-2.5)
  
  abline(a=0,b=1,lty=1,col='grey')
  
  mtext(side=1,text = 'FPR',family='Helvetica',font = 1,cex = 0.7,line = 0.9)
  mtext(side=2,text = 'TPR',family='Helvetica',font = 1,cex = 0.7,line = 0.9)
}
saveFig(file.path(mainDir,'pROC'),plotFun,width = 3,height = 3) 


avg="none",spread.estimate="stderror",xaxs='i',yaxs='i',
xlab="False Positive Rate",ylab="True Positive Rate",add=T
mtext(paste('AUC =',format(mean(perf$auc), digits=3), '±', format(sd(perf$auc), digits=3),
            '  Model: ',model.names[i]),
      side=3, adj=1, line=lines[i], cex=cex, col=fg.col)
abline(a=0,b=1,lty=1,col=bg.col)

title(main=paste(c('10-fold cross validation', paste(selected.features, collapse = ' + ')),
                 collapse = '\n'), cex.main = 1.5)


sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
## definition shape.
plot(sens.ci, type="bars")

