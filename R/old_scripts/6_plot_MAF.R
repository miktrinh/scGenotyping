# This script contains code to plot the output of copyKat and inferCNV (overlapping shifts in average expression level)

setwd('~/lustre_mt22/CN_methods/')

#############
# Libraries #
#############
library(tidyverse)
library(Seurat)
library(readxl)
library(alleleIntegrator)
source('scripts/finalScripts/R/misc.R')

#########################
# Set Global parameters #
#########################
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]
mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
ai_outDir = file.path(mainDir,'alleleIntegrator_output_completed')
# 1. Set plotDir - where the figures are saved
plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots'

# Figure 2C and S2 - Example of AlleleIntegrator output for normal and tumour RCC cells ####
plotAndSave_MAF_AIoutput = function(cov=500,subCl.minSegLen = 2e7,minSegLen=1e6,chrToPlot = c(1:22),skipIfExists = T,projMani,mainDir,ai_outDir,plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots',chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
                                    tumourType_toPlot=NULL,PDID_toPlot = NULL,plotFile_prefix=NULL){
  if(is.null(tumourType_toPlot)){
    tumourType_toPlot = unique(projMani$TumourType)
  }
  if(is.null(PDID_toPlot)){
    PDID_toPlot = unique(projMani$PDID[projMani$TumourType %in% tumourType_toPlot])
  }
  
  
  if(!dir.exists(plotDir)){
    print('Making new Dir')
    dir.create(plotDir,recursive = T)
  }
  
  # Import chromInfo
  chromInfo = read.delim(chromInfo_path,sep = '\t')
  for(tumourType in tumourType_toPlot){
    if(!tumourType %in% c('NB','Wilms','ATRT','Ewings')){
      next
    }
      
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
      
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
        
        if(!current_PDID %in% PDID_toPlot){
          next
        }
        
        message(sprintf('Plotting alleleIntegrator MAF output for Sample %s - tumourType: %s',current_PDID,tumourType))
        # subset annnb.srat object to keep only cells of that sample
        srat.sub = subset(srat, subset = PDID == current_PDID)
        
        
        ####------------------ Get MAF output ----------------####
        gCnts_fp = file.path(ai_outDir,tumourType,current_PDID,paste0(current_PDID,'_gCnts_allhSNPs.RDS'))
        if(skipIfExists & file.exists(gCnts_fp)){
          gCnts = readRDS(gCnts_fp)
        }else{
          phCnts = readRDS(file.path(ai_outDir,tumourType,current_PDID,paste0(current_PDID,'_phCnts.RDS')))
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
          saveRDS(gCnts,file.path(ai_outDir,tumourType,current_PDID,paste0(current_PDID,'_gCnts_allhSNPs.RDS')))  
        }
        
        
        
        ctmp = aggregateByLists(gCnts, assays = c("matCount", "patCount"), gCnts$clusterID)  
        
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
        out2$clusterID = factor(out2$clusterID,levels = unique(c(unique(out2$clusterID),'Tumour')))
        
        #----- Plotting AlleleIntegrator MAF results! -------#
        out2 = out2[out2$chr != 23,]
        chromInfo2 = chromInfo[chromInfo$chrom != 23,]
        
        ####------------------ Generate Battenberg CN summary file ----------------####
        donorMani = projMani[projMani$PDID == current_PDID,]
        btb.fp = unique(donorMani$battenbergFp)
        #----- Processing Battenberg data -------#
        #dna.data = annotateBTB(btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'allelicRatio')  
        dna.data = processBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=chrToPlot,removeBalancedSegs=F,longFormat = T,keepClonalSegs = T,method = 'allelicRatio')  
        
        # Remove X chromosome 
        dna.data = dna.data[dna.data$chr != 23,]
        dna.data$clusterID = 'Tumour'
        for(ct in unique(out2$clusterID)){
          if(ct != 'Tumour'){
            tmp = rbind(data.frame(chr=1,matNum=1,patNum=1,clonalType='maj',totCN=2,tot2min='2:1',tumFrac=0.5,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,clusterID=ct),
                        data.frame(chr=1,matNum=1,patNum=1,clonalType='maj',totCN=2,tot2min='2:1',tumFrac=0.5,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,clusterID=ct))
            dna.data = rbind(dna.data,tmp)
          }
        }
        
        
        subCl.dna = dna.data[dna.data$clonalType == 'sub',]
        majCl.dna = dna.data[dna.data$clonalType != 'sub',]
        
        
        plotFun = function(noFrame=T,noPlot=FALSE){
          layout(mat=matrix(c(1:nlevels(out2$clusterID)),ncol=1),
                 heights = rep(2,nlevels(out2$clusterID)))
          
          
          mar = c(0.1,0.6,0.8,0.6)
          text.cex2 = 0.7
          text.cex1 = 0.7
          text.y = 1.37
          
          for(celltype in levels(out2$clusterID)){
            
            par(mar=mar,xpd=TRUE)
            tmp = out2[out2$clusterID == celltype,]
            dna = majCl.dna[majCl.dna$clusterID == celltype,]
            dna = dna[order(dna$abspos_kb,decreasing = F),]
            
            ncells = nrow(srat.sub@meta.data[srat.sub@meta.data$PDID == current_PDID & srat.sub@meta.data$annot == celltype,])
            
            
            # Plot main frame
            plot(out2$midPos*1000, out2$MAF.readCovBin,
                 las=1,
                 type='n',xaxt='n',yaxt='n',
                 ylim=c(-0.1,1.3),
                 frame.plot=F)
            
            text(x=1e3,y=text.y,paste0(celltype,'_',current_PDID),cex=text.cex1,family = 'Helvetica',font=2,adj = 0)
            text(x=2.86e9,y=text.y,paste0('n=',ncells),cex=text.cex2,family = 'Helvetica',font=1,adj = 1)
            axis(2,at=c(0,0.5,1),labels = c(0,'1/2',1),las=1,pos = 0,tck = -.02,lwd = 0.3,cex.axis=0.7,hadj = 0.2,padj = 0.5)
            axis(4,at=c(0,0.5,1),labels = c(0,'1/2',1),las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.3,cex.axis=0.7,hadj = 1.0,padj = 0.5,col.axis='#b02c46')
            
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
            
            # Plot chromosome number
            #text(x=(xleft+xright)/2,y = 1.095,labels = c(1:22),cex = c(rep(0.7,10),rep(0.62,4),rep(0.53,4),rep(0.31,4)),font = 1)
            
            # Plot AI MAF
            a = ifelse(celltype %in% c('Tumour','Leukocytes'),0.7,0.85)
            points(tmp$midPos*1000,tmp$MAF.readCovBin,
                   pch=19,
                   cex=0.03,col=colAlpha('black',a))
            
            # Plot ground truth
            # Subclone
            if(celltype == 'Tumour' & (nrow(subCl.dna) > 0)){
              for(chr in unique(subCl.dna$chr)){
                for(posID  in unique(subCl.dna[subCl.dna$chr == chr,]$posID )){
                  lines(x=subCl.dna[subCl.dna$chr == chr & subCl.dna$posID  == posID ,]$abspos_kb*1000,subCl.dna[subCl.dna$chr == chr & subCl.dna$posID  == posID ,]$tumFrac,col='#4169E1',lwd=1.1)      
                }
              }
            }
            
            # Major clone CN profile
            lines(x=dna$abspos_kb*1000,dna$tumFrac,col='#b02c46',lwd=1.0)
          }
        }
        
        
        #saveFig(file.path(outdir,paste0('FigS2_',tumourType,'_MAF_',current_PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500,rawData = out2)  
        saveFig(file.path(plotDir,paste0(plotFile_prefix,tumourType,'_MAF_',current_PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500,rawData = out2)  
        
        
      }
    
  }
}








      
      








if(sampleList != 'all' & sheet == 'RCC_mani'){
  saveFig(file.path(res,paste0('v2_RCC_AI_',PDID)),plotFun,width = 2.65,height = 2.15,res=500)
}else if(sampleList == 'all' & sheet == 'RCC_mani'){
  saveFig(file.path(res,paste0('FigS2_AI_',PDID,'_',subCl.minSegLen)),plotFun,width = 2.3,height = 3.5,res=500)  
}else if(sampleList != 'all' & sheet == 'NB_mani'){
  saveFig(file.path(res,paste0('v2_NB_AI_',PDID)),plotFun,width = 2.5,height = 2.265,res=500)  
}else if(sampleList == 'all' & sheet == 'NB_mani'){
  saveFig(file.path(res,paste0('FigS2_AI_',PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500)  
}
plotMAF = function(tumour,sampleList='all'){
  if(tumour == 'RCC'){
    sheet = 'RCC_mani'
  }else if(tumour == 'NB'){
    sheet = 'NB_mani'
  }
  
  # Import Manifest
  projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = sheet)
  if(sheet == 'RCC_mani'){
    projMani = projMani[!grepl('^n_',projMani$SampleID),]  
  }
  outdir = '/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/'
  
  #PDID = unique(projMani$PDID)[1]
  for(PDID in unique(projMani$PDID)){
    donorMani = projMani[projMani$PDID == PDID,]
    
    ####------------------ Get MAF output ----------------####
    if(sheet == 'RCC_mani'){
      f = list.files(outdir,pattern = paste0(PDID,"_gCnts_allhSNPs.RDS"))  
    }else if(sheet == 'NB_mani'){
      f = list.files(outdir,pattern = paste0(PDID,"tmp_gCnts_allhSNPs.RDS"))
    }
    
    gCnts = readRDS(file.path(outdir,f))
    ctmp = aggregateByLists(gCnts, assays = c("matCount", "patCount"), gCnts$clusterID)  
    
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
    out2 = out2[out2$chr != 23,]
    chromInfo2 = chromInfo[chromInfo$chrom != 23,]
    
    if(sheet == 'RCC_mani'){
      out2$clusterID[out2$clusterID == 'Renal_cell_carcinoma'] = 'Tumour'
      out2$clusterID[out2$clusterID == 'Proximal_tubuluar_cells'] = 'PTC'
      out2$clusterID = factor(out2$clusterID,levels = c('Tumour','PTC','Leukocytes'))  
      if(sampleList != 'all'){
        out2=out2[out2$clusterID != 'Leukocytes',]
        out2$clusterID = factor(out2$clusterID,levels = c('Tumour','PTC'))  
      }
    }else if(sheet == 'NB_mani'){
      out2$clusterID = factor(out2$clusterID,levels = c('Tumour','Mesenchyme','Endothelium','Leukocytes')) 
    }
    
    
    ####------------------ Generate Battenberg CN summary file ----------------####
    btb.fp = unique(donorMani$battenbergFp)[!is.na(unique(donorMani$battenbergFp))]
    #----- Processing Battenberg data -------#
    dna.data = annotateBTB(btb.fp,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'allelicRatio')  
    
    # Remove X chromosome 
    dna.data = dna.data[dna.data$Chr != 23,]
    if(sheet=='NB_mani'){
      dna.data$clusterID = 'Tumour'
      tmp = rbind(data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,clusterID='Endothelium'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,clusterID='Endothelium'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,clusterID='Leukocytes'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,clusterID='Leukocytes'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,clusterID='Mesenchyme'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,clusterID='Mesenchyme'))
    }else if(sheet=='RCC_mani'){
      dna.data$clusterID = 'Tumour'
      tmp = rbind(data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,clusterID='PTC'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,clusterID='PTC'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,clusterID='Leukocytes'),
                  data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,clusterID='Leukocytes'))
    }
    
    dna.data = rbind(dna.data,tmp)
    subCl.dna = dna.data[dna.data$type == 'sub',]
    majCl.dna = dna.data[dna.data$type == 'maj',]
    
    
    plotFun = function(noFrame=T,noPlot=FALSE){
      layout(mat=matrix(c(1:nlevels(CNA_summary_byCellType$celltype)),ncol=1),
             heights = rep(2,nlevels(CNA_summary_byCellType$celltype)))
      
      # Set plotting params
      params = data.frame(sheet = c('NB_mani','RCC_mani','NB_mani','RCC_mani'),
                          sampleList=c('all','all','PD42184','PD37228'),
                          mar = c('0.1,0.6,0.8,0.6','0.1,0.6,0.8,0.6','0.5,0.8,0.2,0.6','0.4,0.8,0.2,0.6'),
                          text.cex2=c(0.7,0.7,0.7,0.7),
                          text.cex1=c(0.7,0.7,0.7,0.7),
                          text.y = c(1.37,1.37,1.32,1.32))
      
      for(celltype in levels(out2$clusterID)){
        par(mar=as.vector(as.numeric(strsplit(params[params$sheet == sheet & params$sampleList == sampleList,]$mar,split = ',')[[1]])),xpd=TRUE)
        tmp = out2[out2$clusterID == celltype,]
        dna = majCl.dna[majCl.dna$clusterID == celltype,]
        dna = dna[order(dna$abspos_kb,decreasing = F),]
        
        if(sheet == 'RCC_mani'){
          ncells = nrow(rcc.srat@meta.data[rcc.srat@meta.data$PDID == PDID & rcc.srat@meta.data$finalAnn == celltype,])  
        }else if (sheet == 'NB_mani'){
          ncells = nrow(nb.srat@meta.data[nb.srat@meta.data$PD_ID == PDID & nb.srat@meta.data$finalAnn == celltype,])  
        }
        
        # Plot main frame
        plot(out2$midPos*1000, out2$MAF.readCovBin,
             las=1,
             type='n',xaxt='n',yaxt='n',
             ylim=c(-0.1,1.3),
             frame.plot=F)
        
        text(x=1e3,y=ifelse('all'%in% sampleList,1.37,1.32),paste0(celltype,'_',PDID),cex=params[params$sheet == sheet & params$sampleList == sampleList,]$text.cex1,family = 'Helvetica',font=2,adj = 0)
        text(x=2.86e9,y=ifelse('all'%in% sampleList,1.37,1.32),paste0('n=',ncells),cex=params[params$sheet == sheet & params$sampleList == sampleList,]$text.cex2,family = 'Helvetica',font=1,adj = 1)
        axis(2,at=c(0,0.5,1),labels = c(0,'1/2',1),las=1,pos = 0,tck = -.02,lwd = 0.3,cex.axis=0.7,hadj = 0.2,padj = 0.5)
        axis(4,at=c(0,0.5,1),labels = c(0,'1/2',1),las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.3,cex.axis=0.7,hadj = 1.0,padj = 0.5,col.axis='#b02c46')
        
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
        
        # Plot chromosome number
        #text(x=(xleft+xright)/2,y = 1.095,labels = c(1:22),cex = c(rep(0.7,10),rep(0.62,4),rep(0.53,4),rep(0.31,4)),font = 1)
        
        # Plot AI MAF
        a = ifelse(celltype %in% c('Tumour','Leukocyte'),0.7,0.85)
        points(tmp$midPos*1000,tmp$MAF.readCovBin,
               pch=19,
               cex=0.03,col=colAlpha('black',a))
        
        # Plot ground truth
        # Subclone
        if(celltype == 'Tumour' & (nrow(subCl.dna) > 0)){
          for(chr in unique(subCl.dna$Chr)){
            for(Idx in unique(subCl.dna[subCl.dna$Chr == chr,]$Idx)){
              lines(x=subCl.dna[subCl.dna$Chr == chr & subCl.dna$Idx == Idx,]$abspos_kb*1000,subCl.dna[subCl.dna$Chr == chr & subCl.dna$Idx == Idx,]$tumFrac,col='#4169E1',lwd=1.1)      
            }
          }
        }
        
        # Major clone CN profile
        lines(x=dna$abspos_kb*1000,dna$tumFrac,col='#b02c46',lwd=1.0)
      }
    }
    
    if(sampleList != 'all' & sheet == 'RCC_mani'){
      saveFig(file.path(res,paste0('v2_RCC_AI_',PDID)),plotFun,width = 2.65,height = 2.15,res=500)
    }else if(sampleList == 'all' & sheet == 'RCC_mani'){
      saveFig(file.path(res,paste0('FigS2_AI_',PDID,'_',subCl.minSegLen)),plotFun,width = 2.3,height = 3.5,res=500)  
    }else if(sampleList != 'all' & sheet == 'NB_mani'){
      saveFig(file.path(res,paste0('v2_NB_AI_',PDID)),plotFun,width = 2.5,height = 2.265,res=500)  
    }else if(sampleList == 'all' & sheet == 'NB_mani'){
      saveFig(file.path(res,paste0('FigS2_AI_',PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500)  
    }
  }
}
