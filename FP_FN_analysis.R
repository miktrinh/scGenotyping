# Sensitivity and Specificity Analysis - Comparing CopyKat inferred CNprofile to ground truth
#############
# Libraries #
#############
library(tidyverse)

# Get ChromInfo
chromInfo = read.delim('~/CN_method_tmp/chrom_abspos_kb.txt',sep = '\t')

# Import NB and RCC sratObj
if(tumour == 'RCC'){
  sheet = 'RCC_mani'
  srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')
  copykat.results = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v2_rcc.CK.normREF.default.2397/v2_CKresults_default_normREF_80perc_2397.rds')
  # Extract CK output for RCC
  CKdata = srat@meta.data %>% group_by(PDID,finalAnn,CKpred.normREF.default.80perc.2397) %>% summarise(nCells = n())
  colnames(CKdata) = c('PD_ID','cell_type','CopyKat_output','nCells')
  CKdata$PD_ID = factor(CKdata$PD_ID,levels = c("PD35918","PD36793","PD37104","PD37228"))
}else if(tumour == 'NB'){
  sheet = 'NB_mani'
  srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')
  copykat.results = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_nb.CK.normREF.default.2397/v5_CKresults_default_normREF_80perc_2397.rds')
  # Extract CK output for NB
  CKdata = srat@meta.data %>% group_by(PD_ID,cell_type,CKpred.normREF.default.80perc.2397) %>% summarise(nCells = n())
  colnames(CKdata) = c('PD_ID','cell_type','CopyKat_output','nCells')
  CKdata$PD_ID = factor(CKdata$PD_ID,levels = c("PD42184","PD42752-1","PD42752-2","PD46693","PD43255"))
}

# Import Manifest
projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = sheet)
if(sheet == 'RCC_mani'){
  projMani = projMani[!grepl('^n_',projMani$SampleID),]  
}



DimPlot(srat,group.by = 'finalAnn')

CKdata$CopyKat_output = ifelse(is.na(CKdata$CopyKat_output),'Uncalled',
                               ifelse(CKdata$CopyKat_output == 'aneuploid','Aneuploid','Diploid'))
######################################################
## Setting arbitary threshold of 0.2
binSize=5e7
threshold = 0.2
main.srat = srat
#----- Processing CopyKat results -------#
out = data.frame()

for(k in 1:length(copykat.results)){
  
  # Get copyKat CNA matrix and prediction
  CNA_summary_byCellType = data.frame()
  CNA_mat = copykat.results[[k]]$CNAmat
  colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
  colnames(CNA_mat) = gsub('\\.','-',colnames(CNA_mat))
  pred = as.data.frame(copykat.results[[k]]$prediction)
  
  sample = names(copykat.results)[k]
  
  if(length(sample) > 1){
    message(paste0('More than 1 sample detected: k=',k,', samples are ',sample))
  }else{
    message(paste0('Checking sample ',sample))
  }
  
  # subset srat object to keep only cells of that sample
  srat = subset(main.srat, subset = PDID == sample)
  
  
  ####------------------ Generate Battenberg CN summary file ----------------####
  donorMani = projMani[projMani$PDID == sample,]
  btb.fp = unique(donorMani$battenbergFp)
  PDID=sample
  #----- Processing Battenberg data -------#
  dna.data = annotateBTB(btb.fp,subCl.minSegLen = 2e7,PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat=F,method = 'totalCN')  
  
  dna.data$cna = ifelse(dna.data$tot2min == '2:1',FALSE,TRUE)
  segs = dna.data[dna.data$type=='maj',]
  # subset by annotated cell type
  for(celltype in unique(srat$finalAnn)){
    
    CNA_mat_sub = CNA_mat[,c(1:3,which(colnames(CNA_mat) %in% rownames(srat@meta.data[srat@meta.data$finalAnn == celltype,])))]
    
    if(ncol(CNA_mat_sub) == 4){
      chrom_tmp=data.frame(celltype = celltype,CNA_mat_sub)
      colnames(chrom_tmp)[5] = 'mean_logCN'
    }else if (ncol(CNA_mat_sub) > 4){
      chrom_tmp = data.frame(celltype = celltype,CNA_mat_sub[,c(1:3)],mean_logCN = apply(CNA_mat_sub[,-c(1:3)],MARGIN = 1,FUN = mean))
    }else if (ncol(CNA_mat_sub) < 4){
      chrom_tmp = data.frame(celltype = celltype,CNA_mat_sub[,c(1:3)],mean_logCN = NA)
    }
    
    CNA_summary_byCellType = rbind(CNA_summary_byCellType,chrom_tmp)
  }
  
  
  #----- Plotting copykat results! -------#
  
  for(celltype in unique(CNA_summary_byCellType$celltype)){
    for(chr in unique(chromInfo$chrom)){
      if(!chr %in% unique(segs$Chr)){
        next
      }
      chrsegs = segs[segs$Chr == chr,]
      chrsegs = chrsegs[order(chrsegs$Start,decreasing = F),]
      chr.cna = CNA_summary_byCellType[CNA_summary_byCellType$chrom == chr &
                                         CNA_summary_byCellType$celltype == celltype,]
      
      for(j in 1:nrow(chrsegs)){
        
        seg = chrsegs[j,]
        stops = seq(seg$Start,seg$Stop,binSize)
        starts = seq(seg$Start,seg$Stop,binSize)
        cna = seg$cna
        if(length(stops) > 0){
          stops = c(stops[-1],seg$Stop)
        }
        for(i in 1:length(stops)){
          start = starts[i]
          stop = stops[i]
          tmp = chr.cna[chr.cna$chrompos > start & chr.cna$chrompos<=stop,]
          if(nrow(tmp) == 0){
            tmp2 = data.frame(chr,start,stop,meanScore=0,cna=cna,celltype=celltype,PDID=sample)
          }else if(nrow(tmp) > 0){
            tmp2 = data.frame(chr,start,stop,meanScore = mean(tmp$mean_logCN,na.rm=T),cna=cna,celltype=celltype,PDID=sample)
          }
          out = rbind(out,tmp2)
        }
        
      }
    }
    
  }
  
}


out$absMeanScore = abs(out$meanScore)
out$predCall = ifelse(out$absMeanScore > 0.2,'yes','no')

d = out %>% group_by(PDID,celltype,cna,predCall) %>% tally()

d$truth = ifelse(d$celltype == 'Tumour',as.character(d$cna),FALSE)

d$result = ifelse(d$truth == TRUE & d$predCall == 'yes','TP',
                  ifelse(d$truth == TRUE & d$predCall == 'no','FN',
                         ifelse(d$truth == FALSE & d$predCall == 'yes','FP',
                                ifelse(d$truth == FALSE & d$predCall == 'no','TN','?'))))



saveRDS(list(out,d),'/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_nb.CK.Quant_v2.RDS')


