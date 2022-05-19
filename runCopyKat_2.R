# Run Copy Kat
setwd('~/lustre_mt22/CN_methods/')

#############
# Libraries #
#############
# Load libraries
library(copykat)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap)
library(readxl)
source('scripts/finalScripts/misc.R')



#########################
# Set Global parameters #
#########################
seed = 2397
normREF.fracToUsed = 1
n.cores=20
skipIfExists=T
normREF=T
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]  
mainDir = '~/lustre_mt22/CN_methods/revision_2204/'

for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('Wilms','ATRT','Ewings')){
    # Import annotated Seurat objects
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    
    # 1. define outdir
    if(normREF){
      outDir = file.path(mainDir,'CopyKAT_output/normREF',tumourType)  
    }else{
      outDir = file.path(mainDir,'CopyKAT_output/noNormREF',tumourType)  
    }
    
    if(!dir.exists(outDir)){
      print('Making new Dir')
      dir.create(outDir,recursive = T)
    }else{
      message('Existing output directory found!')
    }
    # if normREF is specified, add 'normREF' tag to the output file / column title
    if(normREF){
      col_to_add = paste0('CKpred.normREF.default.',normREF.fracToUsed*100,'perc.',seed) 
      srat[[paste0('CK.normREF.',seed,'.',normREF.fracToUsed*100,'perc')]] = NA
    }else{
      col_to_add = paste0('CKpred.default.',seed)  
    }
    
    # Add output column to seurat object
    srat[[col_to_add]] = NA
    
    # Run copyKat for each sample
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Running copyKat for Sample %s - tumourType: %s',current_PDID,tumourType))
      
      # Subset for relevant info
      donorMani = projMani[projMani$PDID == current_PDID,]
      sub.srat = subset(srat,subset = PDID == current_PDID)
      
      # if normREF is specified, add 'normREF' tag to the output file / column title
      if(normREF){
        normREF.cellIDs = rownames(sub.srat@meta.data[sub.srat@meta.data$annot == 'Leukocytes',])  
        # Subsampling normal reference cellIDs (if normREF.fracToUsed < 1, otherwise just reshuffle)
        set.seed(seed)
        norm.cell.names = sample(normREF.cellIDs,size = length(normREF.cellIDs)*normREF.fracToUsed)
        # Mark on the big seurat object which cells were used as REF
        m = match(norm.cell.names,rownames(srat@meta.data))
        sum(is.na(m))
        srat[[paste0('CK.normREF.',seed,'.',normREF.fracToUsed*100,'perc')]][m,1] = TRUE  
        
      }else{
        normREF.cellIDs = NULL  
        norm.cell.names = ''
      }
      
      
      # 3. Run CopyKat - WITH NORMAL REFERENCE
      # Get matrix of raw counts
      mtx = sub.srat@assays$RNA@counts  
      
      set.seed(seed)
      
      pred1_path = file.path(outDir,paste0(tumourType,'_',current_PDID,'_CK.',seed,'_prediction.csv'))
      if(file.exists(pred1_path) & skipIfExists){
        pred = read.csv(file.path(outDir,paste0(tumourType,'_',current_PDID,'_CK.',seed,'_prediction.csv')))
      }else{
        copykat.results <- copyKat_inhouse_v2(rawmat=mtx, id.type="S", plot.genes = T,
                                sam.name=paste0(outDir,'/',current_PDID),
                                norm.cell.names = norm.cell.names,n.cores=n.cores)
        #copykat.results = copykat_inhouse(rawmat=mtx, id.type="S", plot.genes = "TRUE",
        #                                  sam.name=paste0(outDir,'/',current_PDID), 
        #                                  n.cores=n.cores, norm.cell.names = norm.cell.names, seed = seed)  
        #saveRDS(copykat.results,outPath)  
        pred = as.data.frame(copykat.results$prediction)
        #m = match(rownames(pred),gsub('-','.',rownames(srat@meta.data)))
        pred1 = pred[pred$cell.names %in% pred$cell.names[pred$copykat.pred !='not.defined'],]
        pred2 = pred[!pred$cell.names %in% pred1$cell.names,]
        print(nrow(pred) == (nrow(pred1) + nrow(pred2)))
        pred = rbind(pred1,pred2)
        
        write_csv(pred,file.path(outDir,paste0(tumourType,'_',current_PDID,'_CK.',seed,'_prediction.csv')))  
      }
      
      
      m = match(pred$cell.names,rownames(srat@meta.data))
      message(sprintf('Sample %s: found %d cellID mismatches!',current_PDID,sum(is.na(m))))
      srat[[col_to_add]][m,1] = pred[,'copykat.pred']  
    }
    
    DimPlot(srat,group.by = col_to_add)
    
    # Save the seurat object?
    saveRDS(srat,file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    srat@meta.data$cellID = rownames(srat@meta.data)
    if(is.null(norm.cell.names)){
      write_csv(as.data.frame(srat@meta.data),file.path(outDir,paste0(tumourType,'_CK.',seed,'_metadata.csv')))  
    }else{
      write_csv(as.data.frame(srat@meta.data),file.path(outDir,paste0(tumourType,'_CK.normREF.',seed,'.',normREF.fracToUsed*100,'perc_metadata.csv')))
    }
    
  }
}


subCl.minSegLen = 2e7
# Import chromInfo
chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('Wilms','ATRT','Ewings')){
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    
    # 1. Set outdir
    # Make output dir
    if(normREF){
      outDir = file.path(mainDir,'CopyKAT_output/normREF',tumourType)  
    }else{
      outDir = file.path(mainDir,'CopyKAT_output/noNormREF',tumourType)  
    }
    
    if(!dir.exists(outDir)){
      warning(sprintf('No existing output directory found for PDID %s - Tumourtype: %s ... Please check!',PDID,tumourType))
    }
    
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Plotting copyKat output for Sample %s - tumourType: %s',current_PDID,tumourType))
      # Get copyKat CNA matrix and prediction
      CNA_mat = read.delim(file.path(outDir,paste0(current_PDID,'_copykat_CNA_results.txt')))
      colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
      colnames(CNA_mat) = gsub('\\.','-',colnames(CNA_mat))
      
      pred = read.delim(file.path(outDir,paste0(tumourType,'_',current_PDID,'_CK.',seed,'_prediction.csv')))
      #pred = pred1
      CNA_summary_byCellType = data.frame()
      
      
      # subset annnb.srat object to keep only cells of that sample
      srat.sub = subset(srat, subset = PDID == current_PDID)
      
      # subset by annotated cell type
      for(celltype in unique(srat.sub$annot)){
        CNA_mat_sub = CNA_mat[,c(1:3,which(colnames(CNA_mat) %in% gsub('\\.','-',rownames(srat.sub@meta.data[srat.sub@meta.data$annot == celltype,]))))]
        
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
      
      # Remove X chromosome 
      CNA_summary_byCellType = CNA_summary_byCellType[CNA_summary_byCellType$chrom != 23,]
      CNA_summary_byCellType$celltype[CNA_summary_byCellType$celltype == 'Tumour Cells'] = 'Tumour'
      if(tumourType == 'Wilms'){
        CNA_summary_byCellType$celltype = factor(CNA_summary_byCellType$celltype,levels=c("Tumour", 'Leukocytes','Endothelium',"Muscle"))  
      }else if(tumourType == 'Ewings'){
        CNA_summary_byCellType$celltype = factor(CNA_summary_byCellType$celltype,levels=c("Tumour", 'Leukocytes',"Endothelium"))
      }else if(tumourType == 'ATRT'){
        CNA_summary_byCellType$celltype = factor(CNA_summary_byCellType$celltype,levels=c("Tumour", 'Leukocytes',"Endothelium","Fibroblast","Epithelium"))
      }
      
      
      
      ####------------------ Generate Battenberg CN summary file ----------------####
      donorMani = projMani[projMani$PDID == current_PDID,]
      btb.fp = unique(donorMani$battenbergFp)
      #----- Processing Battenberg data -------#
      dna.data = annotateBTB(btb.fp = btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'totalCN')  
      # Remove X chromosome 
      dna.data = dna.data[dna.data$Chr != 23,]
      
      dna.data$celltype = 'Tumour'
      for(ct in unique(CNA_summary_byCellType$celltype)){
        if(ct != 'Tumour'){
          tmp = rbind(data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype=ct),
                      data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype=ct))
          dna.data = rbind(dna.data,tmp)
        }
      }
      
      
      
      dna.data$log_CNratio = log(as.numeric(dna.data$tumTot)/2)
      subCl.dna = dna.data[dna.data$type == 'sub',]
      majCl.dna = dna.data[dna.data$type == 'maj',]
      
      
      #----- Plotting copykat results! -------#
      chromInfo2 = chromInfo[chromInfo$chrom != 23,]
      
      plotFun = function(noFrame=TRUE,noPlot=FALSE){
        noFrame=T
        # Set layout
        layout(mat=matrix(c(1:length(unique(CNA_summary_byCellType$celltype))),ncol=1),
               heights = rep(2,length(unique(CNA_summary_byCellType$celltype))))

        for(celltype in unique(CNA_summary_byCellType$celltype)[order(unique(CNA_summary_byCellType$celltype))]){
          
          if(length(unique(CNA_summary_byCellType$celltype)) == 2){
            par(mar=c(2.7,0.6,2.7,0.6),xpd=TRUE)  
          }else{
            par(mar=c(0.5,0.6,0.5,0.6),xpd=TRUE) 
          }
          
          tmp = CNA_summary_byCellType[CNA_summary_byCellType$celltype == celltype,]
          dna = majCl.dna[majCl.dna$celltype == celltype,]
          ncells = nrow(srat.sub@meta.data[srat.sub@meta.data$PDID == current_PDID & srat.sub@meta.data$annot == celltype,])
          
          # Set params for plotting
          ylim = c(round(min(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)-0.1,digits = 1),round(max(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)+0.2,digits = 1))
          ybottom=min(round(min(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)-0.1,digits = 1),round(log(0.5)/2,2))
          ytop=max(round(max(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)+0.1,digits = 1),round(max(dna.data$log_CNratio)/2,2))
          ytext = ytop + 0.1
          
          
          # Plot main frame
          plot(CNA_summary_byCellType$abspos, CNA_summary_byCellType$mean_logCN,
               las=1,
               type='n',
               ylim=ylim,
               xlab=ifelse(noFrame,'','Genomic Position'),
               ylab=ifelse(noFrame,'',''),
               xaxt=ifelse(noFrame,'n','s'),
               yaxt=ifelse(noFrame,'n','s'),
               frame.plot=F)
          
          
          #Plot background chromosome
          xleft = c(0,chromInfo2[chromInfo2$arm == 'q' & chromInfo2$chrom!=22,]$abspos*1000)
          xright = c(chromInfo2[chromInfo2$arm == 'q',]$abspos*1000)
          
          
          # Plot y axes
          if(ytop >0 & ytop < 0.35){ # Wilms
            scale.factor = 2
            axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.7,cex.axis=0.6,hadj = 1.5,col='black',
                 at=c(round(log(0.5)/2,2),0,round(log(3/2)/2,2)),col.axis = '#b02c46',
                 labels = c(1,2,3))  
          }else if(ytop >= 0.63){ # PD47706 (Ewings)
            scale.factor = 4
            axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.7,cex.axis=0.6,hadj = 1.5,col='black',
                 at=c(round(log(0.5)/4,2),0,round(log(4/2)/4,2),round(log(6/2)/4,2),round(log(7/2)/4,2)),col.axis = '#b02c46',
                 labels = c(1,2,4,6,7))  
            ytop = ytop/2
            ytext = ytop + 0.1/4

          }else if(ytop >= 0.35){ # ATRT
            scale.factor = 2
            axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.7,cex.axis=0.6,hadj = 1.5,col='black',
                 at=c(round(log(0.5)/2,2),0,round(log(3/2)/2,2),round(log(4/2)/2,2)),col.axis = '#b02c46',
                 labels = c(1,2,3,4))  
          }
          ybottom=min(round(min(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)-0.1,digits = 1),round(log(0.5)/scale.factor,2))
          axis(2,at=c(ybottom+0.05,0,ytop-0.05),labels = c('Low',0,'High'),las=1,pos = 0,tck = -.00,lwd = 0.3,cex.axis=0.6,hadj = 0.3,padj = 0.5)
          
          text(x=9e8,y=ytext,paste0(celltype,'_',current_PDID,' (n=',ncells,')'),cex=0.6,family = 'Helvetica',font=2,adj = 0)
          
          col = replicate(c('white','lightgrey'),n = 22/2)
          rect(xleft=xleft,
               xright=xright,
               ybottom=ybottom,
               ytop=ytop,
               col = col,
               lty = 'blank')
          #Black surrounding border
          rect(xleft=min(xleft),
               xright=max(xright),
               ybottom=ybottom,
               ytop=ytop,
               col = colAlpha('white',0.0001),
               border = 'black',lwd = 0.4)
          
          #text(x=(xleft+xright)/2,y = 0.72,labels = c(1:22),cex = c(rep(0.7,10),rep(0.62,4),rep(0.53,4),rep(0.31,4)),font = 1)
          
          # Plot ground truth
          rect(xleft = dna[dna$posType=='Start',]$abspos_kb*1000,
               xright = dna[dna$posType=='Stop',]$abspos_kb*1000,
               ybottom = 0,
               ytop=dna[dna$posType=='Start',]$log_CNratio/scale.factor,
               col=colAlpha('#b02c46',1),
               border=colAlpha('#b02c46',1),lwd = 0.8)
          
          # Plot subclone total CN
          if(celltype == 'Tumour cells' & nrow(subCl.dna) > 0){
            for(chr in unique(subCl.dna$Chr)){
              for(logR in unique(subCl.dna[subCl.dna$Chr == chr,]$log_CNratio)){
                d = subCl.dna[subCl.dna$Chr == chr & subCl.dna$log_CNratio == logR,]
                if(logR == 0){
                  lines(d$abspos_kb*1000,d$log_CNratio/2,col='#4169E1',lwd=1.0)
                }else{
                  rect(xleft = d[d$posType=='Start',]$abspos_kb*1000,
                       xright = d[d$posType=='Stop',]$abspos_kb*1000,
                       ybottom = 0,
                       ytop=d[d$posType=='Start',]$log_CNratio/2,
                       col=colAlpha('#4169E1',0.7),
                       border=colAlpha('#4169E1',0.7),lwd = 0)
                }
              }
            }
          }
          
          # Plot CopyKat output
          lines(x=tmp$abspos,tmp$mean_logCN,col='black',lwd=0.7)
          segments(x0=min(xleft),x1 = max(xright), 
                   y0=0.2, y1=0.2,
                   col = 'darkgrey',lty = 'dashed',lwd = 0.4)
          segments(x0=min(xleft),x1 = max(xright), 
                   y0=-0.2, y1=-0.2,
                   col = colAlpha('darkgrey',0.9),lty = 'dashed',lwd = 0.4)
        }
      }
      
      saveFig(file.path(outDir,paste0('../FigS2_',tumourType,'_CK_',current_PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500)  
    }
  }
}
  

















