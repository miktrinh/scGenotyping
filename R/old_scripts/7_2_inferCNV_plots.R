# This script contains code to plot the output of copyKat and inferCNV (overlapping shifts in average expression level)

setwd('~/lustre_mt22/CN_methods/')

#############
# Libraries #
#############
library(tidyverse)
library(Seurat)
library(readxl)
# sudo apt install jags
library(infercnv)

source('scripts/finalScripts/R/misc.R')

#########################
# Set Global parameters #
#########################
normREF=T



# Import chromInfo
chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')

# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]  
mainDir = '~/lustre_mt22/CN_methods/revision_2204/'
tumourType_toPlot = unique(projMani$TumourType)
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
      
      # Scale around zero
      #inferCNVsummary_byCellType$mean_logCN = log(inferCNVsummary_byCellType$mean_logCN)# - 1
      #inferCNVsummary_byCellType$mean_logCN = inferCNVsummary_byCellType$mean_logCN - 1
      
      
      
      
      
      
      
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











run.final.obj = readRDS('./revision_2204/inferCNV_output/normREF/Wilms/PD48777/run.final.infercnv_obj')
plot_cnv(run.final.obj, k_obs_groups = 1, 
         cluster_by_groups = T, cluster_references = T, 
         out_dir = './revision_2204/inferCNV_output/test', x.center = 1, 
         x.range = 'auto', title = "inferCNV_test", 
         output_filename = "inferCNV_test", output_format = 'png', 
         write_expr_matrix = TRUE, png_res = 300, useRaster = T)
expr = read.delim('./revision_2204/inferCNV_output/test/expr.inferCNV_test.dat',sep = '\t')

# plot heatmap, where row = cells, col = genes, split by annotation
mtx = t(as.matrix(expr))
toKeep = sample(1:nrow(mtx),size = nrow(mtx)*0.1)
mtx = mtx[toKeep,]
m = match(rownames(mtx), gsub('-','.',rownames(srat@meta.data)))
sum(is.na(m))
type = srat@meta.data$annot[m]


library(ComplexHeatmap)
Heatmap(mtx,split = type,cluster_rows = T,cluster_columns = F,show_row_dend = F,show_column_dend = F,
        show_row_names = F,show_column_names = F)


# Convert this heatmap to wigle line
subCl.minSegLen = 2e7




























