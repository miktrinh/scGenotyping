#' Generates the Figures for paper

#################
# Materials map #
#################

## Main Figures ##

#F1A - Summary of methods.
#F1B - Histogram of SNV coverage.
#F1C - Histogram of hSNP coverage.
#F1D - hSNPs alelle frequency in normal vs tumour samples.



#F2A - RCC UMAP coloured by PDID + small panels for expr of CA9.
#F2B - Barplots of AI and CK tumour/normal classification for RCC.
#F2C - Example of alleleIntegrator output (MAF) and CopyKat output (average expression) for PTC and tumour RCC cells (PD35918)

#F2D - NB UMAP coloured by PDID.
#F2E - Barplots of AI and CK tumour/normal classification for NB
#F2F - Example of alleleIntegrator output (MAF) and CopyKat output for Mesenchyme and tumour NB cells (PD42752-2 and PD42184)



#F3A - Barplots of AI and CK tumour/normal classification for Wilms, Ewings, and ATRT.
#F3B - ROC for AI, CK and inferCNV across all samples.
#F3C - Distribution of MAF for AI, and of log_Expression for CK and inferCNV.
#F3D - Subclonal analysis for PD46693 (DNA-derived BAF at hSNPs, AI-derived MAF for major and minor clones)


## Supplementary Figures ##

#FS1 - CopyKat and AI output for RCC dataset
#FS2 - CopyKat and AI output for NB dataset + MA plot for PD46693a majorCl/minorCl

## Supplementary Table ##
#TS1-2 - List of DEGs between PD46693 majorCl vs minorCl 
#TS3:11 - Battenberg output for all samples from RCC + NB datasets



###################
# Set Project Dir #
###################
setwd('/lustre/scratch117/casm/team274/mt22/CN_methods/')



#############
# Libraries #
#############
library(tidyverse)
library(plotrix)
library(Seurat)
library(readxl)
library(alleleIntegrator)
#Plotting libraries
#library(ggplot2)
#library(plotrix)
#library(reshape2)
#library(beeswarm)
library(RColorBrewer)
#library(survminer)
#source('scripts/plotParts.R')
source('scripts/finalScripts/R/misc.R')

####################
# Useful functions #
####################
plotAndSave_MAF_AIoutput = function(cov=500,subCl.minSegLen = 2e7,minSegLen=1e6,chrToPlot = c(1:22),skipIfExists = T,
                                    projMani,mainDir,ai_outDir,
                                    plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots',
                                    chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
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
    if(!tumourType %in% c('RCC','NB','Wilms','ATRT','Ewings')){
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


plot_CK_inferCNV_avgExpr = function(projMani,mainDir,inferCNV_dir = file.path(mainDir,'inferCNV_output'),
                                    CK_dir = file.path(mainDir,'CopyKAT_output'),subCl.minSegLen = 2e7,minSegLen = 1e6,
                                    plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots',
                                    chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
                                    skipIfExists = T,
                                    normREF=T,chrToPlot = c(1:22),
                                    tumourType_toPlot=NULL,PDID_toPlot = NULL,plotFile_prefix=NULL){
  # Import chromInfo
  chromInfo = read.delim(chromInfo_path,sep = '\t')
  
  
  for(tumourType in tumourType_toPlot){
    if(!tumourType %in% c('NB','RCC','Wilms','ATRT','Ewings')){next}
    
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    
    # 1. Set outdir
    # Make output dir
    if(normREF){
      outDir = file.path(plotDir,'normREF')  
    }else{
      outDir = file.path(plotDir,'noNormREF')  
    }
    
    if(!dir.exists(outDir)){
      print('Making new Dir')
      dir.create(outDir,recursive = T)
    }
    
    
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      if(!current_PDID %in% PDID_toPlot){
        next
      }
      
      message(sprintf('Plotting inferCNV and copyKat output for Sample %s - tumourType: %s',current_PDID,tumourType))
      # subset annnb.srat object to keep only cells of that sample
      srat.sub = subset(srat, subset = PDID == current_PDID)
      
      
      #### 1. Import inferCNV output, ie. final heatmap matrix ####
      if(normREF){
        expr.mtx.fp = file.path(inferCNV_dir,'normREF',tumourType,current_PDID,'expr.infercnv.dat')
        gene_order_fp = file.path(inferCNV_dir,'normREF',tumourType,current_PDID,'run.final.infercnv_obj')
      }else{
        expr.mtx.fp = file.path(inferCNV_dir,'noNormREF',tumourType,current_PDID,'expr.infercnv.dat')
        gene_order_fp = file.path(inferCNV_dir,'noNormREF',tumourType,current_PDID,'run.final.infercnv_obj')
      }
      expr.mtx = read.delim(expr.mtx.fp,sep = '\t')
      
      # Get absolute genomic position for each gene
      chromInfo.sub = chromInfo[chromInfo$arm == 'q',]
      chromInfo.sub$chrom = chromInfo.sub$chrom + 1  
      # Import gene_order
      gene_order = readRDS(gene_order_fp)
      gene_order = gene_order@gene_order
      gene_order$chrom = gsub('chr','',gene_order$chr)
      gene_order = gene_order[gene_order$chrom %in% chrToPlot,]
      # Calculate absolute genomic position for each gene
      gene_order$prev_end = 0
      m = match(gene_order$chrom[gene_order$chrom != 1],chromInfo.sub$chrom)
      sum(is.na(m))
      gene_order$prev_end[gene_order$chrom != 1] = chromInfo.sub$abspos_kb[m]
      gene_order$abspos = ((gene_order$start+gene_order$stop)/2) + 1000*gene_order$prev_end
      
      #### Calculate the mean expression shifts for each gene across all cells belonging to a given celltype
      inferCNVsummary_byCellType = data.frame()
      for(celltype in unique(srat.sub$annot)){
        CNA_mat_sub = expr.mtx[,c(which(gsub('^X','',colnames(expr.mtx)) %in% gsub('-','.',rownames(srat.sub@meta.data[srat.sub@meta.data$annot == celltype,]))))]
        chrom_tmp = data.frame(celltype = celltype,gene = rownames(CNA_mat_sub),mean_logCN = apply(CNA_mat_sub,MARGIN = 1,FUN = function(x){mean(log(x))}))
        inferCNVsummary_byCellType = rbind(inferCNVsummary_byCellType,chrom_tmp)
      }
      
      # Add gene order and genomic position info
      inferCNVsummary_byCellType = inferCNVsummary_byCellType[inferCNVsummary_byCellType$gene %in% rownames(gene_order),]
      m = match(inferCNVsummary_byCellType$gene,rownames(gene_order))
      sum(is.na(m))
      inferCNVsummary_byCellType$chrom = gene_order$chrom[m]
      inferCNVsummary_byCellType$abspos = gene_order$abspos[m]
      inferCNVsummary_byCellType = inferCNVsummary_byCellType[order(inferCNVsummary_byCellType$abspos,decreasing = F),]
      # Remove X chromosome 
      inferCNVsummary_byCellType = inferCNVsummary_byCellType[inferCNVsummary_byCellType$chrom != 23,]
      inferCNVsummary_byCellType$celltype[inferCNVsummary_byCellType$celltype == 'Tumour Cells'] = 'Tumour'
      if(tumourType == 'Wilms'){
        inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c("Tumour", 'Leukocytes','Endothelium',"Muscle"))  
      }else if(tumourType == 'Ewings'){
        inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c("Tumour", 'Leukocytes',"Endothelium"))
      }else if(tumourType == 'ATRT'){
        inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c("Tumour", 'Leukocytes',"Endothelium","Fibroblast","Epithelium"))
      }else if(tumourType == 'RCC'){
        inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c('Leukocytes',"PTC",'Tumour'))
      }
      

      
      
      #### 2. Import copyKat CNA matrix and prediction ####
      if(normREF){
        cna.mtx.fp = file.path(CK_dir,'normREF',tumourType,paste0(current_PDID,'_copykat_CNA_results.txt'))
      }else{
        cna.mtx.fp = file.path(CK_dir,'noNormREF',tumourType,paste0(current_PDID,'_copykat_CNA_results.txt'))
      }
      
      CNA_mat = read.delim(cna.mtx.fp)
      colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
      
      CNA_summary_byCellType = data.frame()
      # subset by annotated cell type
      for(celltype in unique(srat.sub$annot)){
        CNA_mat_sub = CNA_mat[,c(1:3,which(colnames(CNA_mat) %in% gsub('-','.',rownames(srat.sub@meta.data[srat.sub@meta.data$annot == celltype,]))))]
        
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
      }else if(tumourType == 'RCC'){
        inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c('Leukocytes',"PTC",'Tumour'))
      }
      
      
      
      
      
      
      
      ####------------------ Generate Battenberg CN summary file ----------------####
      donorMani = projMani[projMani$PDID == current_PDID,]
      btb.fp = unique(donorMani$battenbergFp)
      #----- Processing Battenberg data -------#
      #dna.data = annotateBTB(btb.fp = btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'totalCN')  
      # Remove X chromosome 
      #dna.data = dna.data[dna.data$Chr != 23,]
      
      dna.data = processBTB(btb.fp,PDID=current_PDID,minSegLen=minSegLen,subCl.minSegLen = subCl.minSegLen,tgtChrs=chrToPlot,removeBalancedSegs=F,longFormat = T,keepClonalSegs = T,method = 'totalCN')  
      
      dna.data$celltype = 'Tumour'
      for(ct in unique(CNA_summary_byCellType$celltype)){
        if(ct != 'Tumour'){
          tmp = rbind(data.frame(chr=1,matNum=1,patNum=1,clonalType='maj',totCN=2,tumFrac=0.5,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype=ct),
                      data.frame(chr=1,matNum=1,patNum=1,clonalType='maj',totCN=2,tumFrac=0.5,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype=ct))
          dna.data = rbind(dna.data,tmp)
        }
      }
      
      
      
      dna.data$log_CNratio = log(as.numeric(dna.data$totCN)/2)
      subCl.dna = dna.data[dna.data$clonalType == 'sub',]
      majCl.dna = dna.data[dna.data$clonalType != 'sub',]
      
      
      #----- Plotting copykat results! -------#
      chromInfo2 = chromInfo[chromInfo$chrom != 23,]
      
      plotFun = function(noFrame=FALSE,noPlot=FALSE){
        
        # Set layout
        layout(mat=matrix(c(1:length(unique(CNA_summary_byCellType$celltype))),ncol=1),
               heights = rep(2,length(unique(CNA_summary_byCellType$celltype))))
        
        for(celltype in unique(CNA_summary_byCellType$celltype)[order(unique(CNA_summary_byCellType$celltype))]){
          
          if(length(unique(CNA_summary_byCellType$celltype)) == 2){
            par(mar=c(2.7,0.6,2.7,0.6),xpd=TRUE)  
          }else{
            par(mar=c(0.6,0.6,0.8,0.6),xpd=TRUE) 
          }
          
          copykat.data = CNA_summary_byCellType[CNA_summary_byCellType$celltype == celltype,]
          infercnv.data = inferCNVsummary_byCellType[inferCNVsummary_byCellType$celltype == celltype,]
          dna = majCl.dna[majCl.dna$celltype == celltype,]
          ncells = nrow(srat.sub@meta.data[srat.sub@meta.data$PDID == current_PDID & srat.sub@meta.data$annot == celltype,])
          
          # Set params for plotting
          ylim.ck = c(round(min(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)-0.1,digits = 1),round(max(CNA_summary_byCellType[!is.na(CNA_summary_byCellType$mean_logCN),]$mean_logCN)+0.1,digits = 1))
          ylim.inf = c(round(min(inferCNVsummary_byCellType[!is.na(inferCNVsummary_byCellType$mean_logCN),]$mean_logCN)-0.1,digits = 1),round(max(inferCNVsummary_byCellType[!is.na(inferCNVsummary_byCellType$mean_logCN),]$mean_logCN)+0.1,digits = 1))
          ylim = c(min(ylim.ck[1],ylim.inf[1]),max(ylim.ck[2],ylim.inf[2]))
          ybottom=min(ylim[1],round(log(0.5)/2,2))
          ytop=max(ylim[2],round(max(dna.data$log_CNratio)/2,2))
          
          
          
          # Plot main frame
          plot(CNA_summary_byCellType$abspos, CNA_summary_byCellType$mean_logCN,
               las=1,
               type='n',
               ylim=ylim,
               #xlab=ifelse(noFrame,'','Genomic Position'),
               #ylab=ifelse(noFrame,'',''),
               #xaxt=ifelse(noFrame,'n','s'),
               #yaxt=ifelse(noFrame,'n','s'),
               xlab = '',ylab='',xaxt='n',yaxt='n',
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
            ybottom = ybottom/2
            #ytext = ytop + 0.1/4
            
          }else if(ytop >= 0.35){ # ATRT
            scale.factor = 4
            axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.7,cex.axis=0.6,hadj = 1.5,col='black',
                 at=c(round(log(0.5)/4,2),0,round(log(3/2)/4,2),round(log(4/2)/4,2)),col.axis = '#b02c46',
                 labels = c(1,2,3,4))  
            ytop = ytop/2
            ybottom = ybottom/2
          }
          
          #ybottom=min(ylim[1],round(log(0.5)/scale.factor,2))
          axis(2,at=c(ybottom+0.05,0,ytop-0.05),labels = c(ybottom+0.05,0,ytop-0.05),las=1,pos = 0,tck = -.00,lwd = 0.3,cex.axis=0.6,hadj = 0.3,padj = 0.5)# labels = c('Low',0,'High'),
          ytext = ytop + 0.1/scale.factor
          
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
          if(celltype == 'Tumour' & nrow(subCl.dna) > 0){
            for(chr in unique(subCl.dna$chr)){
              for(logR in unique(subCl.dna[subCl.dna$chr == chr,]$log_CNratio)){
                d = subCl.dna[subCl.dna$chr == chr & subCl.dna$log_CNratio == logR,]
                if(logR == 0){
                  lines(d$abspos_kb*1000,d$log_CNratio/2,col='#4169E1',lwd=1.0)
                }else{
                  rect(xleft = d[d$posType=='Start',]$abspos_kb*1000,
                       xright = d[d$posType=='Stop',]$abspos_kb*1000,
                       ybottom = 0,
                       ytop=d[d$posType=='Start',]$log_CNratio/scale.factor,
                       col=colAlpha('#4169E1',0.7),
                       border=colAlpha('#4169E1',0.7),lwd = 0)
                }
              }
            }
          }
          
          
          #### 4. Plot average expression shifts output ####
          ### CopyKat output
          lines(x=copykat.data$abspos,copykat.data$mean_logCN,col='black',lwd=0.4)
          ### inferCNV output
          lines(x=infercnv.data$abspos,infercnv.data$mean_logCN,col='blue',lwd=0.4)
          
        }
      }
      copykat.data$method = 'CK'
      infercnv.data$method = 'inferCNV' 
      dd = rbind(copykat.data[,colnames(copykat.data)[-3]],infercnv.data[,colnames(copykat.data)[-3]])
      
      saveFig(file.path(outDir,paste0(plotFile_prefix,tumourType,'_CK.inferCNV_',current_PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500,rawData = dd)  
    }
  }  
}




#########################
# Set Global parameters #
#########################
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]
mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
refGenome = '/lustre/scratch119/realdata/mdt1/team78pipelines/reference/Human/GRCH37d5/genome.fa'
nParallel=25

plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots'
# first version of all plots are at revision_2204/plots
#resDir = '~/lustre_mt22/CN_methods/revision_2204_v2/Plots/'

if(!dir.exists(plotDir)){dir.create(plotDir)}


# Figure 1B - SNVs Coverage ####
fig1b_snvHist = function(){
  all.out = read.csv('/lustre/scratch117/casm/team274/mt22/CN_methods/revision_2204/snvCov_allout.csv')
  
  all.out$cat = ifelse(all.out$totalReads <5,all.out$totalReads,
                       ifelse(all.out$totalReads <=10,'5-10','>10'))
  dd=all.out %>% group_by(tumCat,cat) %>% summarise(nCells = n())
  dd$cat = factor(dd$cat,levels = c('0','1','2','3','4','5-10','>10'))
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    ylim_max = max(dd$nCells)
    ylab = round(seq(0,ylim_max+10,length.out = 5),digits = -2)
    
    par(mar=c(1.0,2.3,0.2,0.5),xpd=T)
    plot(seq(0,9,1),rep(c(0),10),ylim=c(-3,ylim_max),xlim=c(0,10),
         las=1,
         type='n',frame.plot = F,axes = F)
    if(!noFrame){
      #axis(2,at=c(1,0.5,0),las=1,pos = 0.1,tck=-0.02,cex.axis=0.75,lwd.ticks = 0,hadj = 0.5)
      axis(side = 2,at =ylab,labels = ylab,
           tck=-0.015,lwd = 0.8,cex.axis=0.4,las=1,pos = -0.2,lwd.ticks = 0.5,hadj = 0.05)
      axis(side = 1,at = c(0,10),labels = c('',''),
           tck=-0.001,lwd = 0.8,cex.axis=0.3,las=1,pos = -10,lwd.ticks = 0,hadj = 0.4)
      
      text(x=seq(0.55,9.55,1.5),y=-350,as.character(levels(dd$cat)),cex = 0.4)
      
      mtext(side=1,text = '# Reads',family='Helvetica',font = 1,cex = 0.55,line = 0.1)
      mtext(side=2,text = '# Cells',family='Helvetica',font = 1,cex = 0.55,line = 1)
      
    }
    
    
    
    if(!noPlot){
      ### Bar Plot for pediatric cancer
      xleft = seq(0,9,1.5)
      xright = xleft + 0.5
      ybottom = 0
      ytop= dd[dd$tumCat=='pediatric',]$nCells[order(dd[dd$tumCat=='pediatric',]$cat)]
      
      rect(xleft=xleft,
           xright=xright,
           ybottom=ybottom,
           ytop=ytop,
           col='white',
           lwd = 0.7,
           border = 'black')
      
      
      ### Bar Plot for Adult cancer (RCC)
      xleft = seq(0.6,9.6,1.5)
      xright = xleft + 0.5
      ybottom = 0
      ytop= dd[dd$tumCat=='adult',]$nCells[order(dd[dd$tumCat=='adult',]$cat)]
      
      rect(xleft=xleft,
           xright=xright,
           ybottom=ybottom,
           ytop=ytop,
           col='black',
           lwd = 0.7,
           border = 'black')
      
      cols = c('Pediatric cancer' = 'white','Adult cancer'='black')
      legend(y=9000, x=5.5,legend=names(cols),fill = cols,lwd = 0,cex = 0.45,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 1.0,bty = 'n')
      legend(y=7000, x=5.5,legend=sprintf('Total cells = %d',sum(dd$nCells)),lwd = 0,cex = 0.4,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 1.0,bty = 'n')
    }
    
  }
  saveFig(file.path(plotDir,paste0('Fig1b_SNVcov_distr_allCells')),plotFun,width = 2.1,height = 2.2,rawData = dd)   
}


# Figure 1C - hSNPs Coverage ####
fig1c_hsnpHist = function(){
  all.out = read.csv('/lustre/scratch117/casm/team274/mt22/CN_methods/revision_2204/hSNPcov_allout.csv')
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(1.2,2.5,0.5,1.5),xpd=T)
    n_cells = n_distinct(all.out$cellID)
    dat = hist(log10(all.out$totalReads),plot = F,breaks = 50)
    
    
    
    if(!noPlot){
      hist(log10(all.out$totalReads),breaks = 50,freq = T,main = '',axes = F)  
    }
    
    if(!noFrame){
      
      if(noPlot){
        plot(0,0,xlim=c(1,5),
             ylim=c(0,3000),
             las=1,
             type='n',frame.plot = F,axes = F)
      }
      
      axis(1,at=seq(round(min(dat$breaks)),round(max(dat$breaks))),las=1,pos = 0.05,tck=-0.02,cex.axis=0.7,lwd.ticks = 0.8,hadj = 0.5,padj = -2.5,lwd = 0.8)
      axis(2,at=seq(0,round(max(dat$counts)),1000),las=1,tck=-0.02,cex.axis=0.7,lwd.ticks = 0.4,hadj = 0.45,padj = 0.5,lwd = 0.8)
      title(sprintf('Heterozygous SNP coverage in %s cells',prettyNum(n_cells,big.mark = ',')),family='Helvetica',cex.main = 0.6)
      mtext(side=1,text = '# Reads (log 10)',family='Helvetica',font = 1,cex = 0.7,line = 0.2)
      mtext(side=2,text = '# Cells',family='Helvetica',font = 1,cex = 0.7,line = 1.6)
      
    }
    
  }
  
  saveFig(file.path(plotDir,paste0('Fig1c_hSNPcov_distr')),plotFun,width = 2.7,height = 2.6,rawData = all.out)  
}


# Figure 1D - hSNP BAF for RCC sample PD35918 ####
# code extracted from bulkDNA_BAF.R
fig1d_pd35918_hSNP_baf = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2'
  projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
  PDID = 'PD35918'
  tumourType = 'RCC'
  # Generate BAF plot for bulkDNA PD35918
  message(sprintf('Running AlleleIntegrator for Sample %s - tumourType: %s',PDID,tumourType))
  # Set output directory
  outDir = file.path(mainDir,'alleleIntegrator_output_completed',tumourType,PDID)
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
  
  baf.out = generateCoverageAndBAF(BAM = tumourDNA,refGenome = refGenome,#hSNPs=hSNPs,
                                   outPath = file.path(outDir,paste0(PDID,'_cov_BAF_SNPs1k_v2.tsv')),nParallel=nParallel)
  minCoverage=10
  #Filter to just the ones that we trust
  filt = baf.out[baf.out$coverage>=minCoverage,]
  # Subset randomly 50% of the points
  set.seed(2397)
  idx = sample(1:nrow(mcols(filt)), nrow(mcols(filt))/10, replace=FALSE)
  idx = idx[order(idx,decreasing = F)]
  filt.sub = filt[idx,]
  filt.sub$hSNP_pos = paste0(seqnames(filt.sub),':',start(filt.sub))
  dd = as.data.frame(mcols(filt.sub))
  dd = dd[,c('hSNP_pos','BAF','coverage', 'isHet','logR','A','C','G','T','Tot','refCnts','altCnts')]
  
  plotFun = function(noFrame=F,noPlot=FALSE,minCoverage=10){
    #Filter to just the ones that we trust
    #filt = baf.out[baf.out$coverage>=minCoverage,]
    #Work out the chromosome boundaries
    chrsToPlot=c(1:22)
    chrs = chrsToPlot
    chrLens = seqlengths(filt)
    tmp = sapply(split(start(filt),as.character(seqnames(filt))),max)
    chrLens[is.na(chrLens)] = tmp[names(chrLens)[is.na(chrLens)]]
    chrLens = as.numeric(chrLens[chrs])
    
    
    x = start(filt) +cumsum(c(0,chrLens))[match(as.character(seqnames(filt)),chrs)]
    
    # Subset randomly 50% of the points
    #set.seed(2397)
    #idx = sample(1:nrow(mcols(filt)), nrow(mcols(filt))/4, replace=FALSE)
    #filt.sub = filt[idx,]
    
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
  
  saveFig(file.path(plotDir,paste0('Fig1d_dnaBAF_',PDID,'_',tumourType)),plotFun,width = 5.8,height = 2.2,res=1000,rawData = dd)
} 


########################################################################################################################

# Figure 2A - RCC UMAP ####
fig2a_rccUMAP = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
  # Import the seurat object
  rcc.srat = readRDS(file.path(mainDir,'sc_seuratObjects','RCC','RCC_ann.RDS'))
  rcc.srat@meta.data$finalAnn = rcc.srat@meta.data$annot
  dd = cbind(rcc.srat@meta.data[,c("cellID","PDID","finalAnn")],rcc.srat@reductions$umap@cell.embeddings)
  
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    colfunc <- colorRampPalette(c("#383838", "#E2E2E2"))
    ccs = colfunc(length(unique(dd$PDID)))
    ccs = c(PD35918=ccs[1],
            PD37228 = ccs[2],
            PD37104 = ccs[3],
            PD36793 = ccs[4]
    )
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         xlim=c(-13,17),
         ylim=c(-13,17), cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Renal Cell Carcinoma'),
         frame.plot=F)
    
    #axis(2,hadj=0.2,padj=0.6,cex.axis=0.6,tck = -.02,las=1,lwd = 0.6,at=c(-10,0,10,15),labels = c(-10,'','',15))
    #axis(1,hadj=0.5,padj=-2.2,cex.axis=0.6,tck = -.02,las=1,lwd = 0.6,at=c(-10,0,10,15),labels = c(-10,'','',15))
    #mtext('UMAP 1',side=1,line = -0.1,family = 'Helvetica',cex = 0.65)
    #mtext('UMAP 2',side=2,line = 0.1,family = 'Helvetica',cex = 0.65)
    
    if(!noPlot){
      #Add density contours
      #addDensityContours(dd$UMAP_1,dd$UMAP_2,dd$finalAnn,col=colAlpha('black',0.4),nGrid = 2000)
      points(dd$UMAP_1,dd$UMAP_2,
             col = ccs[dd$PDID],
             pch = 19,
             cex=0.01)
      
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ finalAnn,data=dd,FUN=mean)
      
      #Position tweaks
      mids[mids$finalAnn=='Leukocyte','UMAP_2'] = mids[mids$finalAnn=='Leukocyte','UMAP_2'] + 4.9
      mids[mids$finalAnn=='Leukocyte','UMAP_1'] = mids[mids$finalAnn=='Leukocyte','UMAP_1'] - 0.3
      mids[mids$finalAnn=='PTC','UMAP_2'] = mids[mids$finalAnn=='PTC','UMAP_2'] - 4.9
      mids[mids$finalAnn=='PTC','UMAP_1'] = mids[mids$finalAnn=='PTC','UMAP_1'] -0.8
      mids[mids$finalAnn=='Tumour','UMAP_2'] = mids[mids$finalAnn=='Tumour','UMAP_2'] + 3
      mids[mids$finalAnn=='Tumour','UMAP_1'] = mids[mids$finalAnn=='Tumour','UMAP_1'] + 0.9
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$finalAnn,cex = 0.7,xpad = 1.2,ypad = 1.9,border = T,
                   #bg=ccs[mids$finalAnn],
                   col=c('black','black','#b02c46'),font=c(1,1,2))
      
      #Add PDID barplot
      pdid = table(dd$PDID,dd$finalAnn)
      pdid = sweep(pdid, 2, colSums(pdid), "/")
      p.mids = aggregate(cbind(UMAP_1,UMAP_2) ~ finalAnn,data=dd,FUN=mean)
      for(i in unique(p.mids$finalAnn)){
        if(i == 'Leukocytes'){
          xshift = -5.6
          yshift = -5.3
        }else if(i=='PTC'){
          xshift = 4.3
          yshift = -2.8
        }else if(i=='Tumour'){
          xshift = 4.9
          yshift = 1.0
        }
        dat = pdid[c('PD35918','PD36793','PD37104','PD37228'),colnames(pdid) == i]*4.5
        xleft = p.mids[p.mids$finalAnn == i,]$UMAP_1 + xshift
        xright = xleft + 1
        ybottom = p.mids[p.mids$finalAnn == i,]$UMAP_2 + yshift + c(0,cumsum(dat)[-length(dat)])
        ytop= ybottom + dat
        
        #rect(xleft=xleft,
        #     xright=xright,
        #     ybottom=ybottom,
        #     ytop=ytop,
        #     col = ccs[names(dat)],lwd = 0.5,
        #     border = 'black')
        
      }
      legend(x=9.5, y=-6.8,legend=names(ccs),fill = ccs,lwd = 0,cex = 0.65,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    }
  }
  
  saveFig(file.path(plotDir,'Fig2a_RCC_PDIDumap'),plotFun,rawData=dd,width = 2.9,height = 2.65,res = 500)
  
  
  
  
  
  ## Plot expression level of key marker genes on UMAP
  nbreak=100
  for(tgt in c('CA9')){
    # Extract marker genes expression level (scale.data)
    dd[[tgt]] = rcc.srat@assays$RNA@scale.data[rownames(rcc.srat@assays$RNA@scale.data) == tgt,]
    #This adds a column of color values based on the expression values
    dd[[paste0(tgt,'.col')]] = ifelse(dd[[tgt]] == 0, 'grey','#b02c46')
    dd[[paste0(tgt,'.col')]] = colAlpha(dd[[paste0(tgt,'.col')]],seq(0.05,1,(1-0.05)/nbreak)[as.numeric(cut(dd[[tgt]],breaks = nbreak))])
  }
  
  
  # Figure 2a.bonus - RCC marker genes expression
  for(tgt in c('CA9')){
    plotFun2 = function(noFrame=FALSE,noPlot=FALSE){
      par(mar=c(0.1,0.1,0.1,0.1))
      plot(dd$UMAP_1,dd$UMAP_2,
           las=1,
           type='n',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           frame.plot=T)
      text(x=13.0,y=15.6,labels = tgt,col='#b02c46',family="Helvetica", font=4,cex=0.7)
      
      if(!noPlot){
        #Add density contours
        #pl$addDensityContours(dd$UMAP_1,dd$UMAP_2,dd$annot,col=pl$colAlpha('black',0.4),nSplits=10)
        points(dd$UMAP_1,dd$UMAP_2,
               col = dd[[paste0(tgt,'.col')]],
               pch = 19,
               cex=0.01)
      }
    }
    
    saveFig(file.path(plotDir,paste0('Fig2a_bonus_RCC_',tgt,'umap')),plotFun2,rawData=dd,width = 2.2,height = 2.2,res = 500)
  }
}


# Figure 2B - Barplots of RCC AI/CK classification ####
fig2b_rccBarplot = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
  # Import the seurat object
  rcc.srat = readRDS(file.path(mainDir,'sc_seuratObjects','RCC','RCC_ann.RDS'))
  rcc.srat@meta.data$finalAnn = rcc.srat@meta.data$annot
  # Extract CK output
  dd = rcc.srat@meta.data
  if('AIcall_normREF_v2' %in% colnames(dd)){
    dd$AIcall_normREF = as.character(dd$AIcall_normREF_v2)
  }
  
  dd$AIcall_normREF = ifelse(is.na(dd$AIcall_normREF),'Uncalled',
                             ifelse(dd$AIcall_normREF == 'abbFrac','Tumour',
                                    ifelse(dd$AIcall_normREF == 'normFrac','Normal','Uncalled')))
  dd$CKpred.normREF.default.100perc.2397 = ifelse(is.na(dd$CKpred.normREF.default.100perc.2397),'Uncalled',
                                                  ifelse(dd$CKpred.normREF.default.100perc.2397 == 'aneuploid','Aneuploid',
                                                         ifelse(dd$CKpred.normREF.default.100perc.2397 == 'diploid','Diploid','Uncalled')))
  dd$PDID = as.factor(dd$PDID)
  dd$finalAnn = factor(dd$finalAnn, levels = c('PTC','Leukocytes','Tumour'))
  dd = dd[,c("cellID","PDID","finalAnn","AIcall_normREF","CKpred.normREF.default.100perc.2397")]
  
  
  for(method in c('AI','CK')){
    if(method == 'AI'){
      col = 'AIcall_normREF'
      cols = c(Normal = '#c8c8c8',
               Tumour = '#b02c46',
               Uncalled = '#474646')  
      cols_umap = c(Normal = 'black',
                    Tumour = '#b02c46',
                    Uncalled = '#dedede')
      
      dd.sub = dd[dd$AIcall_normREF != 'Uncalled',]
      dd.sub$AIcall_normREF = factor(dd.sub$AIcall_normREF,levels = c('Tumour','Normal'))
      
    }else if(method == 'CK'){
      col = 'CKpred.normREF.default.100perc.2397'
      cols = c(Diploid = '#c8c8c8',
               Aneuploid = '#b02c46',
               Uncalled = '#474646')
      
      cols_umap = c(Diploid = 'black',
                    Aneuploid = '#b02c46',
                    Uncalled = '#dedede')
      
      dd.sub = dd[dd$CKpred.normREF.default.100perc.2397 != 'Uncalled',]
      dd.sub$CKpred.normREF.default.100perc.2397 = factor(dd.sub$CKpred.normREF.default.100perc.2397,levels = c('Aneuploid','Diploid'))
    }
  
    
    # Define the layout
    plotFun = function(noFrame=FALSE,noPlot=FALSE){
      
      layout(matrix(c(1,2,3,4),ncol=1),heights = c(0.2,1,1,1.8))
      
      par(mar=c(0,0.6,0.8,0.6))
      plot(0, 0,
           las=1,
           type='n',frame.plot = F,axes = F)
      title('RCC',cex.main=1,family = 'Helvetica',font=2)

      for(celltype in levels(dd.sub$finalAnn)){
        par(mar=c(0.05,1,0.2,0.1),xpd=TRUE)
        
        tmp=as.matrix(table(dd.sub[dd.sub$finalAnn == celltype,col],dd.sub$PDID[dd.sub$finalAnn == celltype]))
        tmp = sweep(tmp,2,colSums(tmp),'/')
        
        if(celltype == levels(dd.sub$finalAnn)[length(levels(dd.sub$finalAnn))]){
          par(mar=c(2.5,1,0.5,0.1),xpd=TRUE)
          barplot(tmp,
                  col=cols[rownames(tmp)],
                  space=0.12,axes = FALSE,
                  las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
          text(x = seq(0.62,5.7,by = 1.13),y = -0.06,colnames(tmp),cex=0.6,family = 'Helvetica',font=1,srt=90,adj = 1)
          text(x = -0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
          text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$finalAnn == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
        }else{
          barplot(tmp,
                  col=cols[rownames(tmp)],
                  space=0.12,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                  las = 1,main = '',border = F)
          text(x =-0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
          text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$finalAnn == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
        }
      }
    }
    saveFig(file.path(plotDir,paste0('Fig2b_RCC_barplot_',method,'_noUncalled')),plotFun,width = (0.17+0.15*n_distinct(dd.sub$PDID)),height = 2.75,rawData = dd.sub,res=500)
  }
}




# Figure 2C - MAF and average expression PD37228 ####
# PlotandSaveMAF function from 6_plot_MAF.R script

fig2c_MAF_CKexpr_PD37228 = function(){
  ai_outDir = file.path(mainDir,'alleleIntegrator_output_completed')
  plotAndSave_MAF_AIoutput(cov=500,subCl.minSegLen = 2e7,minSegLen=1e6,chrToPlot = c(1:22),skipIfExists = T,
                           projMani=projMani,mainDir=mainDir,ai_outDir=ai_outDir,
                           plotDir=plotDir,
                           chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
                           tumourType_toPlot='RCC',PDID_toPlot = 'PD37228',plotFile_prefix='Fig2c_')
  
  
  plot_CK_inferCNV_avgExpr(normREF=T,subCl.minSegLen = 2e7,minSegLen = 1e6,chrToPlot = c(1:22),skipIfExists = T,
                           projMani=projMani,mainDir=mainDir,
                           inferCNV_dir = file.path(mainDir,'inferCNV_output'),
                           CK_dir = file.path(mainDir,'CopyKAT_output'),
                           plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots',
                           chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
                           tumourType_toPlot='RCC',PDID_toPlot = 'PD37228',plotFile_prefix='Fig2c_')
}



# Figure 2D - NB UMAP ####
fig2d_nbUMAP = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
  # Import the seurat object
  nb.srat = readRDS(file.path(mainDir,'sc_seuratObjects','NB','NB_ann.RDS'))
  dd = cbind(nb.srat@meta.data[,c("cellID","PDID","finalAnn")],nb.srat@reductions$umap@cell.embeddings)
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    colfunc <- colorRampPalette(c("#E2E2E2","#303030"))
    ccs = colfunc(length(unique(dd$PDID)))
    ccs = c('PD42752-1'=ccs[1],
            'PD42752-2' = ccs[2],
            PD46693 = ccs[3],
            PD42184 = ccs[4],
            PD43255 = ccs[5])
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         xlim=c(-15.5,15),
         ylim=c(-15,20), cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',
         ylab='',
         main=ifelse(noFrame,'','Neuroblastoma'),
         xaxt='n',
         yaxt='n',
         frame.plot=F)
    
    
    #axis(2,hadj=0.6,padj=0.6,cex.axis=0.8,tck = -.02,las=1,lwd = 0.7,at=c(-15,-10,0,10,15),labels = c('',-10,'',10,''))
    #axis(1,hadj=0.5,padj=-1.3,cex.axis=0.8,tck = -.02,las=1,lwd = 0.7,at=c(-15,-10,0,10,15),labels = c('',-10,'',10,''))
    #axis(2,hadj=0.3,padj=0.5,cex.axis=0.6,tck = -.02,las=1,lwd = 0.7)
    #axis(1,hadj=0.5,padj=-2.5,cex.axis=0.6,tck = -.02,las=1,lwd = 0.7)
    #mtext('UMAP 1',side=1,line = 1,family = 'Helvetica',cex = 0.6)
    #mtext('UMAP 2',side=2,line = 1.2,family = 'Helvetica',cex = 0.6)
    #mtext('UMAP 1',side=1,line = -0.1,family = 'Helvetica',cex = 0.65)
    #mtext('UMAP 2',side=2,line = 0.1,family = 'Helvetica',cex = 0.65)
    
    if(!noPlot){
      #Add density contours
      #addDensityContours(dd$UMAP_1,dd$UMAP_2,dd$finalAnn,col=colAlpha('black',0.4),nGrid = 2000)
      points(dd$UMAP_1,dd$UMAP_2,
             col = ccs[dd$PDID],
             pch = 19,
             cex=0.01)
      
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ finalAnn,data=dd,FUN=mean)
      
      #Position tweaks
      mids[mids$finalAnn=='Leukocyte','UMAP_2'] = mids[mids$finalAnn=='Leukocytes','UMAP_2'] + 4.9
      mids[mids$finalAnn=='Leukocyte','UMAP_1'] = mids[mids$finalAnn=='Leukocytes','UMAP_1'] + 3.7
      mids[mids$finalAnn=='Endothelium','UMAP_1'] = mids[mids$finalAnn=='Endothelium','UMAP_1'] - 0.4
      mids[mids$finalAnn=='Endothelium','UMAP_2'] = mids[mids$finalAnn=='Endothelium','UMAP_2'] - 0.6
      mids[mids$finalAnn=='Mesenchyme','UMAP_1'] = mids[mids$finalAnn=='Mesenchyme','UMAP_1'] + 8
      mids[mids$finalAnn=='Mesenchyme','UMAP_2'] = mids[mids$finalAnn=='Mesenchyme','UMAP_2'] + 3.9
      mids[mids$finalAnn=='Tumour','UMAP_1'] = mids[mids$finalAnn=='Tumour','UMAP_1'] - 0.4
      mids[mids$finalAnn=='Tumour','UMAP_2'] = mids[mids$finalAnn=='Tumour','UMAP_2'] - 4.1
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$finalAnn,cex = 0.7,xpad = 1.2,ypad = 1.9,border = T,
                   bg=ccs[mids$finalAnn],
                   col=c('black','black','black','#b02c46'),font=c(1,1,1,2))
      
      #Add PDID barplot
      pdid = table(dd$PDID,dd$finalAnn)
      pdid = sweep(pdid, 2, colSums(pdid), "/")
      p.mids = aggregate(cbind(UMAP_1,UMAP_2) ~ finalAnn,data=dd,FUN=mean)
      for(i in unique(p.mids$finalAnn)){
        if(i == 'Leukocytes'){
          xshift = 7.3
          yshift = -3
        }else if(i=='Endothelium'){
          xshift = -2.5
          yshift = 1.2
        }else if(i=='Mesenchyme'){
          xshift = -2.2
          yshift = 2.9
        }else if(i=='Tumour'){
          xshift = -7.3
          yshift = -7.5
        }
        dat = pdid[c('PD42184','PD42752-2','PD43255','PD42752-1','PD46693'),colnames(pdid) == i]*5
        xleft = p.mids[p.mids$finalAnn == i,]$UMAP_1 + xshift
        xright = xleft + 1
        ybottom = p.mids[p.mids$finalAnn == i,]$UMAP_2 + yshift + c(0,cumsum(dat)[-length(dat)])
        ytop= ybottom + dat
        #rect(xleft=xleft,
        #     xright=xright,
        #     ybottom=ybottom,
        #     ytop=ytop,
        #     col = ccs[names(dat)],lwd = 0.7,
        #     border = 'black')
        
      }
      legend(x=6, y=-7,legend=names(ccs)[5:1],fill = ccs[5:1],lwd = 0,cex = 0.65,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    }
  }
  
  
  saveFig(file.path(plotDir,'Fig2d_NB_PDIDumap'),plotFun,rawData=dd,width = 3,height = 3,res=500)
}



# Figure 2E - Barplots of NB AI/CK classification ####
fig2e_nbBarplot = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
  # Import the seurat object
  nb.srat = readRDS(file.path(mainDir,'sc_seuratObjects','NB','NB_ann.RDS'))
  
  # Extract CK output
  dd = nb.srat@meta.data
  if('AIcall_normREF_v2' %in% colnames(dd)){
    dd$AIcall_normREF = as.character(dd$AIcall_normREF_v2)
  }
  
  dd$AIcall_normREF = ifelse(is.na(dd$AIcall_normREF),'Uncalled',
                             ifelse(dd$AIcall_normREF == 'abbFrac','Tumour',
                                    ifelse(dd$AIcall_normREF == 'normFrac','Normal','Uncalled')))
  dd$CKpred.normREF.default.100perc.2397 = ifelse(is.na(dd$CKpred.normREF.default.100perc.2397),'Uncalled',
                                                  ifelse(dd$CKpred.normREF.default.100perc.2397 == 'aneuploid','Aneuploid',
                                                         ifelse(dd$CKpred.normREF.default.100perc.2397 == 'diploid','Diploid','Uncalled')))
  dd$PDID = as.factor(dd$PDID)
  dd$finalAnn = factor(dd$finalAnn, levels = c('Endothelium','Mesenchyme','Leukocytes','Tumour'))
  dd = dd[,c("cellID","PDID","finalAnn","AIcall_normREF","CKpred.normREF.default.100perc.2397")]
  
  
  for(method in c('AI','CK')){
    if(method == 'AI'){
      col = 'AIcall_normREF'
      cols = c(Normal = '#c8c8c8',
               Tumour = '#b02c46',
               Uncalled = '#474646')  
      cols_umap = c(Normal = 'black',
                    Tumour = '#b02c46',
                    Uncalled = '#dedede')
      
      dd.sub = dd[dd$AIcall_normREF != 'Uncalled',]
      dd.sub$AIcall_normREF = factor(dd.sub$AIcall_normREF,levels = c('Tumour','Normal'))
      
    }else if(method == 'CK'){
      col = 'CKpred.normREF.default.100perc.2397'
      cols = c(Diploid = '#c8c8c8',
               Aneuploid = '#b02c46',
               Uncalled = '#474646')
      
      cols_umap = c(Diploid = 'black',
                    Aneuploid = '#b02c46',
                    Uncalled = '#dedede')
      
      dd.sub = dd[dd$CKpred.normREF.default.100perc.2397 != 'Uncalled',]
      dd.sub$CKpred.normREF.default.100perc.2397 = factor(dd.sub$CKpred.normREF.default.100perc.2397,levels = c('Aneuploid','Diploid'))
    }
    
    
    # Define the layout
    plotFun = function(noFrame=FALSE,noPlot=FALSE){
      
      layout(matrix(c(1,2,3,4,5),ncol=1),heights = c(0.2,1,1,1,1.8))
      
      par(mar=c(0,0.6,0.8,0.6))
      plot(0, 0,
           las=1,
           type='n',frame.plot = F,axes = F)
      title('NB',cex.main=1,family = 'Helvetica',font=2)
      
      for(celltype in levels(dd.sub$finalAnn)){
        par(mar=c(0.05,1,0.2,0.1),xpd=TRUE)
        
        tmp=as.matrix(table(dd.sub[dd.sub$finalAnn == celltype,col],dd.sub$PDID[dd.sub$finalAnn == celltype]))
        tmp = sweep(tmp,2,colSums(tmp),'/')
        
        if(celltype == levels(dd.sub$finalAnn)[length(levels(dd.sub$finalAnn))]){
          par(mar=c(2.5,1,0.5,0.1),xpd=TRUE)
          barplot(tmp,
                  col=cols[rownames(tmp)],
                  space=0.12,axes = FALSE,
                  las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
          text(x = seq(0.62,5.7,by = 1.13),y = -0.06,colnames(tmp),cex=0.6,family = 'Helvetica',font=1,srt=90,adj = 1)
          text(x = -0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
          text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$finalAnn == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
        }else{
          barplot(tmp,
                  col=cols[rownames(tmp)],
                  space=0.12,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                  las = 1,main = '',border = F)
          text(x =-0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
          text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$finalAnn == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
        }
      }
    }
    saveFig(file.path(plotDir,paste0('Fig2e_NB_barplot_',method,'_noUncalled')),plotFun,width = (0.17+0.15*n_distinct(dd.sub$PDID)),height = 2.75,rawData = dd.sub,res=500)
  }
}




# Figure 2F - MAF and average expression PD42752-2 PD42184 ####
# PlotandSaveMAF function from 6_plot_MAF.R script
fig2f_MAF_CKexpr_NB = function(){
  
  for(PDID in c("PD42752-2",'PD42184')){
    ai_outDir = file.path(mainDir,'alleleIntegrator_output_completed')
    plotAndSave_MAF_AIoutput(cov=500,subCl.minSegLen = 2e7,minSegLen=1e6,chrToPlot = c(1:22),skipIfExists = T,
                             projMani=projMani,mainDir=mainDir,ai_outDir=ai_outDir,
                             plotDir=plotDir,
                             chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
                             tumourType_toPlot='NB',PDID_toPlot = PDID,plotFile_prefix='Fig2f_')
    
    
    plot_CK_inferCNV_avgExpr(normREF=T,subCl.minSegLen = 2e7,minSegLen = 1e6,chrToPlot = c(1:22),skipIfExists = T,
                             projMani=projMani,mainDir=mainDir,
                             inferCNV_dir = file.path(mainDir,'inferCNV_output'),
                             CK_dir = file.path(mainDir,'CopyKAT_output'),
                             plotDir='~/lustre_mt22/CN_methods/revision_2204_v2/Plots',
                             chromInfo_path = '/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',
                             tumourType_toPlot='NB',PDID_toPlot = PDID,plotFile_prefix='Fig2f_')
  }
  
}







########################################################################################################################

# Figure 3A - Barplots of AI/CK classification for Ewings/Wilms/ATRT ####
fig3a_othersBarplot = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
  
  #if(tumourType %in% c('RCC','NB','Wilms','ATRT','Ewings')){
  for(tumourType in c('Wilms','ATRT','Ewings')){
    # Import the seurat object
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    # Extract AI/CK call
    dd = srat@meta.data
    dd$finalAnn = dd$annot
    if('AIcall_normREF_v2' %in% colnames(dd)){
      dd$AIcall_normREF = as.character(dd$AIcall_normREF_v2)
    }
    # Define relevant cell types to keep
    if(tumourType %in% c('ATRT','Ewings','Wilms')){
      celltype_order = c('Endothelium','Leukocytes','Tumour')
    }else if(tumourType == 'RCC'){
      dd$finalAnn = dd$annot
      celltype_order = c('PTC','Leukocytes','Tumour')
    }else if(tumourType == 'NB'){
      celltype_order = c('Endothelium','Mesenchyme','Leukocytes','Tumour')
    }
    
  
  
    dd$AIcall_normREF = ifelse(is.na(dd$AIcall_normREF),'Uncalled',
                               ifelse(dd$AIcall_normREF == 'abbFrac','Tumour',
                                      ifelse(dd$AIcall_normREF == 'normFrac','Normal','Uncalled')))
    dd$CKpred.normREF.default.100perc.2397 = ifelse(is.na(dd$CKpred.normREF.default.100perc.2397),'Uncalled',
                                                    ifelse(dd$CKpred.normREF.default.100perc.2397 == 'aneuploid','Aneuploid',
                                                           ifelse(dd$CKpred.normREF.default.100perc.2397 == 'diploid','Diploid','Uncalled')))
    dd$PDID = as.factor(dd$PDID)
    dd = dd[dd$finalAnn %in% celltype_order,]
    dd$finalAnn = factor(dd$finalAnn, levels = celltype_order)
    
    dd = dd[,c("cellID","PDID","finalAnn","AIcall_normREF","CKpred.normREF.default.100perc.2397")]
    
    
    for(method in c('AI','CK')){
      if(method == 'AI'){
        col = 'AIcall_normREF'
        cols = c(Normal = '#c8c8c8',
                 Tumour = '#b02c46',
                 Uncalled = '#474646')  
        cols_umap = c(Normal = 'black',
                      Tumour = '#b02c46',
                      Uncalled = '#dedede')
        
        dd.sub = dd[dd$AIcall_normREF != 'Uncalled',]
        dd.sub$AIcall_normREF = factor(dd.sub$AIcall_normREF,levels = c('Tumour','Normal'))
        
      }else if(method == 'CK'){
        col = 'CKpred.normREF.default.100perc.2397'
        cols = c(Diploid = '#c8c8c8',
                 Aneuploid = '#b02c46',
                 Uncalled = '#474646')
        
        cols_umap = c(Diploid = 'black',
                      Aneuploid = '#b02c46',
                      Uncalled = '#dedede')
        
        dd.sub = dd[dd$CKpred.normREF.default.100perc.2397 != 'Uncalled',]
        dd.sub$CKpred.normREF.default.100perc.2397 = factor(dd.sub$CKpred.normREF.default.100perc.2397,levels = c('Aneuploid','Diploid'))
      }
      
      
      # Define the layout
      plotFun = function(noFrame=FALSE,noPlot=FALSE){
        
        if(tumourType == 'NB'){
          layout(matrix(c(1,2,3,4,5),ncol=1),heights = c(0.2,1,1,1,1.8))  
        }else{
          layout(matrix(c(1,2,3,4),ncol=1),heights = c(0.2,1,1,1.8))
        }
        
        par(mar=c(0,0.6,0.8,0.6))
        plot(0, 0,
             las=1,
             type='n',frame.plot = F,axes = F)
        title(tumourType,cex.main=1,family = 'Helvetica',font=2)
        
        for(celltype in levels(dd.sub$finalAnn)){
          par(mar=c(0.05,1,0.2,0.1),xpd=TRUE)
          
          tmp=as.matrix(table(dd.sub[dd.sub$finalAnn == celltype,col],dd.sub$PDID[dd.sub$finalAnn == celltype]))
          tmp = sweep(tmp,2,colSums(tmp),'/')
          
          if(celltype == levels(dd.sub$finalAnn)[length(levels(dd.sub$finalAnn))]){
            par(mar=c(2.5,1,0.5,0.1),xpd=TRUE)
            barplot(tmp,
                    col=cols[rownames(tmp)],
                    space=0.12,axes = FALSE,
                    las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
            text(x = seq(0.62,5.7,by = 1.13),y = -0.06,colnames(tmp),cex=0.6,family = 'Helvetica',font=1,srt=90,adj = 1)
            text(x = -0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
            text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$finalAnn == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
          }else{
            barplot(tmp,
                    col=cols[rownames(tmp)],
                    space=0.12,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                    las = 1,main = '',border = F)
            text(x =-0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
            text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$finalAnn == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
          }
        }
      }
      
      saveFig(file.path(plotDir,paste0('Fig3a_',tumourType,'_barplot_',method,'_noUncalled')),plotFun,width = (0.17+0.15*n_distinct(dd.sub$PDID)),height = 2.75,rawData = dd.sub,res=500)
    }
  }
}



# Figure 3B - ROC for all samples ####
fig3b_ROC_inferCNV_CK_AI = function(){
  mainDir = "~/lustre_mt22/CN_methods/revision_2204_v2/"
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
    
    plot(1, type="n", axes=F,
         frame.plot = ifelse(noFrame,F,T),
         xlim=c(0,1),ylim=c(0,1))
    
    if(!noPlot){
      for(method in c('alleleIntegrator','CK','inferCNV')){
        
        lines(y= data$TPR[data$method == method],x=data$FPR[data$method == method],lty=lty[method],
              type = "l",las=1,frame.plot = T,axes = F,xlim=c(0,1),ylim=c(0,1),lwd=1.5)  
        
        mtext(paste('AUC =',format(unique(data$AUC[data$method == method]), digits=3),
                    '  Method: ',method),
              side=3, adj=0.9,line=line,cex=0.6,font = 1)
        line = line - 1
      }
      abline(a=0,b=1,lty=1,col='grey')
    }
    
    
    if(!noFrame){
      axis(side = 2,at =seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,las=1,pos = -0.04,lwd.ticks = 0.5,hadj = 0.1)
      axis(side = 1,at = seq(0,1,0.2),tck=-0.015,lwd = 0.8,cex.axis=0.6,pos=-0.04,las=0,padj=-2.5)  
      
      mtext(side=1,text = 'False Positive Rate',family='sans',font = 1,cex = 0.7,line = 0.9)
      mtext(side=2,text = 'True Positive Rate',family='sans',font = 1,cex = 0.7,line = 0.9)
    }
    
  }
  #saveFig(file.path(plotDir,'sensitivity_analysis/ROC/inhouseROC_allSamples_v2'),plotFun,width = 3,height = 3,rawData = data) 
  saveFig(file.path(plotDir,'Fig3b_inhouseROC_allSamples'),plotFun,width = 3,height = 3,rawData = data) 
}



# Figure 3C.1 - MAF Distribution for all samples ####
fig3c_MAF_distr_AI = function(){
  # Import Manifest
  projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
  # Remove RCC normal samples from projMani
  projMani = projMani[!grepl('^n_',projMani$SampleID),]  
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
  
  ### Plot aggregated distribution across all PDIDs and tumourTypes
  
  # For each tumour type
  # loop through each PDID
  maf_allSamples = tibble()
  for(tumourType in unique(projMani$TumourType)){
    if(tumourType %in% c('NB','RCC','Wilms','ATRT','Ewings')){
      for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
        
        if(file.exists(file.path(mainDir,'sensitivity_analysis/AI_BAF_distribution',paste0(tumourType,'_',current_PDID,'_MAFdistr_rawData.tsv')))){
          agg_gCnts = read_delim(file.path(mainDir,'sensitivity_analysis/AI_BAF_distribution',paste0(tumourType,'_',current_PDID,'_MAFdistr_rawData.tsv')),col_names = T,delim = '\t')
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
  
  maf_allSamples = maf_allSamples[]
  saveFig(file.path(plotDir,'Fig3c_allSamples_AI_MAFdistr'),plotFun,width = 2.4,height = 2.4,rawData = maf_allSamples) 
  
}


# Figure 3C.2 -AvgExpr Distribution for all samples ####
fig3c_AvgExpr_distr_inferCNV_CK = function(){
  
  # Import Manifest
  projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
  # Remove RCC normal samples from projMani
  projMani = projMani[!grepl('^n_',projMani$SampleID),]  
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'

  ### Plot aggregated distribution across all PDIDs and tumourTypes
  
  # For each tumour type
  # loop through each PDID
  
  binSize=5e6
  all_out_meanScore = tibble()
  for(tumourType in unique(projMani$TumourType)){
    if(tumourType %in% c('RCC','Wilms','ATRT','Ewings','NB')){
      for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
        out_meanScore = read_delim(file.path(mainDir,'sensitivity_analysis/expressionFC_distribution',paste0(tumourType,'_',current_PDID,'_meanScore_inferCNV_CK_',binSize,'.txt')),col_names = T,delim = '\t')
        all_out_meanScore = rbind(all_out_meanScore,out_meanScore)
      }
    }
  }
  
  
  all_out_meanScore = all_out_meanScore[!is.na(all_out_meanScore$meanScore),]
  ## Plot distribution of log_expression_fold_change for diploid segments and aneuploid segments
  
  for(method in c('CK',"inferCNV")){
    out_meanScore_sub = all_out_meanScore[all_out_meanScore$method == method,]
    plotFun = function(noFrame=FALSE,noPlot=FALSE){
      
      par(mar=c(1.4,1.2,0.5,1),xpd=T)
      # plot densities
      xmax = max(density(out_meanScore_sub$meanScore)$x)+0.01
      xmin = min(density(out_meanScore_sub$meanScore)$x)-0.01
      d_dipSegs <- density(out_meanScore_sub$meanScore[out_meanScore_sub$celltype == 'all'])
      plot(d_dipSegs,main='',axes=F,xlim=c(xmin,xmax))
      polygon(d_dipSegs, col="white", border="black")
      
      #d_allSegs <- density(out_meanScore_sub$meanScore)
      #plot(d_allSegs,main='',axes=F,xlim=c(xmin,xmax))
      #polygon(d_allSegs, col="white", border="black")
      
      if(sum(out_meanScore_sub$nTot >2) >2){
        d_up <- density(out_meanScore_sub$meanScore[out_meanScore_sub$celltype != 'all' & out_meanScore_sub$nTot > 2])
        lines(d_up)
        polygon(d_up, col=adjustcolor("#9C9C9B",alpha.f=0.7), border="black")
      }
      
      if(sum(out_meanScore_sub$nTot <2) >2){
        d_down <- density(out_meanScore_sub$meanScore[out_meanScore_sub$celltype != 'all' & out_meanScore_sub$nTot < 2])
        lines(d_down)
        polygon(d_down, col=adjustcolor("#E2E1E1",alpha.f=0.7), border="black")  
      }
      
      x_text = ifelse((xmin > log(1/2)) & (xmax < log(3/2)),
                      list(c(round(xmin,digits = 3),0,round(xmax,digits = 3))),
                      ifelse((xmin > log(1/2)) & (xmax > log(3/2)),
                             list(c(round(xmin,digits = 3),0,round(xmax,digits = 3),log(3/2))),
                             ifelse((xmin < log(1/2)) & (xmax > log(3/2)),
                                    list(c(round(xmin,digits = 3),log(1/2)),0,round(xmax,digits = 3),log(3/2))),
                             ifelse((xmin < log(1/2)) & (xmax < log(3/2)),
                                    list(c(round(xmin,digits = 3),log(1/2),0,round(xmax,digits = 3))),print('Weird'))))
      
      if(x_text == 'Weird'){
        stop('Whoops!')
      }
      axis(1,at=x_text[[1]],
           las=1,pos = 0.05,tck=-0.02,cex.axis=0.6,lwd.ticks = 0.8,hadj = 0.5,padj = -2.5,lwd = 0.8)
      axis(2,at=c(0,round(max(d_dipSegs$y))),labels = c('',''),las=1,tck=-0.00,cex.axis=0.7,lwd.ticks = 0.4,hadj = 0.45,padj = 0.5,lwd = 0.8)
      title(sprintf('%s logFC distribution in %s samples from %s tumour types ',method,n_distinct(out_meanScore_sub$PDID),n_distinct(out_meanScore_sub$tumourType)),family='Helvetica',cex.main = 0.6)
      mtext(side=1,text = 'log expression FC per 5Mb',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
      mtext(side=2,text = 'Density',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
    }
    saveFig(file.path(plotDir,paste0('Fig3c_allSamples_',method,'_exprDistr')),plotFun,width = 2.4,height = 2.4,rawData = out_meanScore_sub) 
  }
  
}


# Figure 3D.1 - hSNP DNA-derived BAF for PD46693
# code extracted from bulkDNA_BAF.R
fig3d_pd46693_hSNP_baf = function(){
  mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2'
  projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
  PDID = 'PD46693'
  tumourType = 'NB'
  nParallel=60
  refGenome
  # Generate BAF plot for bulkDNA PD35918
  message(sprintf('Running AlleleIntegrator for Sample %s - tumourType: %s',PDID,tumourType))
  # Set output directory
  outDir = file.path(mainDir,'alleleIntegrator_output_completed',tumourType,PDID)
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
  
  baf.out = generateCoverageAndBAF(BAM = tumourDNA,refGenome = refGenome, hSNPs=hSNPs,
                                   outPath = file.path(outDir,paste0(PDID,'_cov_BAF_july.tsv')),nParallel=nParallel)
  minCoverage=10
  #Filter to just the ones that we trust
  filt = baf.out[baf.out$coverage>=minCoverage,]
  # Filter for chromosome
  filt = filt[seqnames(filt) %in% c(1:22)]
  # Subset randomly 50% of the points
  set.seed(2397)
  idx = sample(1:nrow(mcols(filt)), nrow(mcols(filt))/4, replace=FALSE)
  idx = idx[order(idx,decreasing = F)]
  filt.sub = filt[idx,]
  filt.sub$hSNP_pos = paste0(seqnames(filt.sub),':',start(filt.sub))
  dd = as.data.frame(mcols(filt.sub))
  dd = dd[,c('hSNP_pos','BAF','coverage', 'isHet','logR','A','C','G','T','Tot','refCnts','altCnts')]
  
  plotFun = function(noFrame=F,noPlot=FALSE,minCoverage=10){
    #Filter to just the ones that we trust
    #filt = baf.out[baf.out$coverage>=minCoverage,]
    #Work out the chromosome boundaries
    chrsToPlot=c(1:22)
    chrs = chrsToPlot
    chrLens = seqlengths(filt)
    tmp = sapply(split(start(filt),as.character(seqnames(filt))),max)
    chrLens[is.na(chrLens)] = tmp[names(chrLens)[is.na(chrLens)]]
    chrLens = as.numeric(chrLens[chrs])
    
    
    x = start(filt) +cumsum(c(0,chrLens))[match(as.character(seqnames(filt)),chrs)]
    
    # Subset randomly 50% of the points
    #set.seed(2397)
    #idx = sample(1:nrow(mcols(filt)), nrow(mcols(filt))/4, replace=FALSE)
    #filt.sub = filt[idx,]
    
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
  
  saveFig(file.path(plotDir,paste0('Fig3d_dnaBAF_',PDID,'_',tumourType)),plotFun,width = 5.8,height = 2.2,res=1000,rawData = dd)
} 



# Figure 3D.2 - Subclonal analysis for PD46693 (DNA-derived BAF at hSNPs, AI-derived MAF for major and minor clones)####
fig3d_PD46693_subclones_MAF = function(){
  # Load Libraries
  library(tidyverse)
  library(Seurat)
  
  # Import PD46693 sratObj
  pd46693 = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/NB/PD46693_probAbb.rds')
  DimPlot(pd46693,group.by = 'cell_type')
  FeaturePlot(pd46693,features = 'chr4.probAbberant')
  
  
  # Subset for Tumour population only ####
  tumour = subset(pd46693,subset = cell_type == 'Tumour')
  
  # Clustering
  tumour = NormalizeData(tumour)
  tumour = FindVariableFeatures(tumour)
  tumour = ScaleData(tumour, features = rownames(tumour))
  tumour = RunPCA(tumour, npcs = 75)
  ElbowPlot(tumour, ndims = 75)
  tumour = FindNeighbors(tumour, dims=1:60)
  tumour = FindClusters(tumour,resolution = 3)
  tumour = RunUMAP(tumour, dims=1:55,n.neighbors = 30,min.dist = 0.6)
  
  sum(is.na(tumour$chr4.probAbberant))
  tumour$chr4.probAbberant = as.numeric(tumour$chr4.probAbberant)
  FeaturePlot(tumour,features = 'chr4.probAbberant',cols = brewer.pal(5,'RdBu')[4:1],pt.size = 0.4)+ggtitle('PD46693 - Tumour')
  
  table(tumour$clonalType,tumour$clonalType)
  tumour$clonalType = ifelse(tumour$chr4.probAbberant > 0.99,'Minor',ifelse(tumour$chr4.probAbberant <0.01,'Major','uncalled'))
  
  
  # Allele Integrator Dot plots for the major and minor clones ####
  
  # Import Manifest
  projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = 'NB_mani')
  outdir = '/lustre/scratch117/casm/team274/mt22/CN_methods/alleleIntegrator_output/'
  
  PDID = 'PD46693'
  donorMani = projMani[projMani$PDID == PDID,]
  
  ####------------------ Generate Battenberg CN summary file ----------------####
  # Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
  # read BTB 
  btb.fp = unique(donorMani$battenbergFp)
  #----- Processing Battenberg data -------#
  dna.data = annotateBTB(btb.fp,subCl.minSegLen = 1e7,PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'allelicRatio')  
  
  # Remove X chromosome 
  dna.data = dna.data[dna.data$Chr != 23,]
  
  dna.data$clusterID = 'Tumour'
  
  subCl.dna = dna.data[dna.data$type == 'sub',]
  majCl.dna = dna.data[dna.data$type == 'maj',]
  
  ####------------------ Get MAF output ----------------####
  f = list.files(outdir,pattern = paste0(PDID,"tmp_gCnts_allhSNPs.RDS"))
  
  gCnts = readRDS(file.path(outdir,f))
  m = match(gCnts$cellID,clonal@meta.data$cellID)
  sum(is.na(m))
  gCnts$clonalType = NA
  gCnts$clonalType[!is.na(m)] = clonal@meta.data$clonalType[m[!is.na(m)]]
  
  for(i in c('Major','Minor')){
    message(sprintf('Analyzing %s clone', i))
    gCnts.tmp = gCnts[!is.na(gCnts$clonalType)]
    gCnts.tmp = gCnts.tmp[gCnts.tmp$clonalType == i]
    
    ctmp = aggregateByLists(gCnts.tmp, assays = c("matCount", "patCount"), gCnts.tmp$clusterID)  
    
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
    cov = 500
    
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
    # Remove X chromosome 
    dna.data = dna.data[dna.data$Chr != 23,]
    if(i == 'Major'){
      dna.data[dna.data$Chr == 4,]$patNum = 1
      dna.data[dna.data$Chr == 4,]$tumFrac = 1/2
      dna.data[dna.data$Chr == 4,]$config = '2:1'
    }else if(i == 'Minor'){
      dna.data[dna.data$Chr == 4,]$patNum = 0
      dna.data[dna.data$Chr == 4,]$tumFrac = 1
      dna.data[dna.data$Chr == 4,]$config = '2:0'
    }
    
    dna.data$clusterID = 'Tumour'
    
    
    
    out2 = out2[out2$chr != 23,]
    chromInfo2 = chromInfo[chromInfo$chrom != 23,]
    ncells_total = nrow(pd46693@meta.data[pd46693@meta.data$cell_type == 'Tumour',])#nrow(nb.srat@meta.data[nb.srat@meta.data$PD_ID == PDID & nb.srat@meta.data$finalAnn == celltype,])  
    ncells = length(unique(gCnts.tmp$cellID))
    
    plotFun = function(noFrame=T,noPlot=FALSE){
      celltype='Tumour'
      
      par(mar=c(0.5,0.8,1.5,0.1),xpd=TRUE)
      tmp = out2[out2$clusterID == celltype,]
      
      
      
      
      # Plot main frame
      plot(out2$midPos*1000, out2$MAF.readCovBin,
           las=1,
           type='n',xaxt='n',yaxt='n',
           main = PDID, cex.main = 0.8,
           ylim=c(-0.1,1.3),
           frame.plot=F)
      ncells=824
      ncells_total=1283
      text(x=1e6,y=1.32,paste0(i,' clone (',ncells,'/',ncells_total,')'),cex=0.7,family = 'Helvetica',font=1,adj = 0)
      axis(2,at=c(0,1/2,1),labels = c(0,'1/2',1),las=1,pos = 0,tck = -.02,lwd = 0.3,cex.axis=0.7,hadj = 0.3,padj = 0.5)
      
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
      
      
      # Plot AI MAF
      points(tmp$midPos*1000,tmp$MAF.readCovBin,
             pch=19,
             cex=0.1,col=colAlpha('black',0.8))
      
      # Plot ground truth
      if(i=='Major'){
        lines(x=majCl.dna$abspos_kb*1000,majCl.dna$tumFrac,col='#b02c46',lwd=1.5)  
      }else if(i=='Minor'){
        lines(x=subCl.dna$abspos_kb*1000,subCl.dna$tumFrac,col='#b02c46',lwd=1.5)  
      }
      
    }
    
    saveFig(file.path(plotDir,paste0('Fig3d_PD46693_',i,'clone_',cov,'_new')),plotFun(ncells_total=ncells_total,ncells = ncells),width = 6,height = 1.4,res=500,rawData = out2)
  }
  
}



########

## Supplementary Figure 7 - PD46693 major vs minor MA plot
# Data generated by the script 'subcloneAnalysis_PD46693.R'
figS7_PD46693_MAplot = function(){
  dd = read.delim('/lustre/scratch117/casm/team274/mt22/CN_methods/PD46693_dataForMAplot.txt',header = T,sep = '\t')
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    sigcol = c(yes = 'black',no='lightgrey',red = 'red')
    
    par(mfrow=c(1,1),mar=c(4,4,0.2,0.5))
    plot(dd$A,dd$M,
         las=1,
         type='n',cex.axis=0.6,tck=-0.03,cex.lab=0.7,
         xlab= 'Mean log2 normalized count',ylab='log2 FC',
         frame.plot=T)
    
    segments(x0=-30,x1 = 0, 
             y0=0, y1=0,
             col = 'black')
    points(dd$A,dd$M,pch=19,cex=0.02,col=sigcol[dd$sig])
    points(dd[dd$sig == 'red',]$A,dd[dd$sig == 'red',]$M,pch=19,cex=0.02,col=sigcol[dd[dd$sig == 'red',]$sig])
    text(dd[dd$sig == 'red',]$M~dd[dd$sig == 'red',]$A,
         labels=dd[dd$sig == 'red',]$gene,
         data=dd, cex=0.4, font=1,pos=1)
  }
  
  saveFig(file.path(plotDir,paste0('FigS7_')),plotFun,rawData = dd,width = 4,height = 2.6,res=500,rawData = dd)
  
}










