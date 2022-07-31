# Sensitivity and Specificity Analysis - Comparing CopyKat inferred CNprofile to ground truth

setwd('~/lustre_mt22/CN_methods/revision_2204_v2/')
#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(GenomicRanges)
source('~/lustre_mt22/CN_methods/scripts/finalScripts/R/misc.R')

#########################
# Set Global parameters #
#########################
normREF = T
binSize=5e6
subCl.minSegLen = 2e7
minSegLen=1e6
chrToPlot = c(1:22)
keepClonalSegs = F

# Get ChromInfo
chromInfo = read.delim('~/lustre_mt22/chrom_abspos_kb.txt',sep = '\t')
# Import Manifest
projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]  
mainDir = '~/lustre_mt22/CN_methods/revision_2204/'
# Set output dir for plots
outdir = file.path('~/lustre_mt22/CN_methods/revision_2204_v2/sensitivity_analysis/expressionFC_distribution')

if(!dir.exists(outdir)){
  print('Making new Dir')
  dir.create(outdir,recursive = T)
}

# For each tumour type
# loop through each PDID

for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('RCC','NB','Wilms','ATRT','Ewings')){
    if(tumourType == 'Wilms'){
      celltypes = c('Tumour','Muscle','Leukocytes')
    }else if(tumourType == 'ATRT'){
      celltypes = c('Tumour','Endothelium','Leukocytes')
    }else if(tumourType == 'Ewings'){
      celltypes = c('Tumour','Endothelium','Leukocytes')
    }else if(tumourType == 'RCC'){
      celltypes = c('Tumour','PTC','Leukocytes')
    }else if(tumourType == 'NB'){
      celltypes = c('Tumour','Mesenchyme','Endothelium','Leukocytes')
    }
    
    #### 1. Set outdir
    #if(normREF){
    #  outDir = file.path(mainDir,'plots/normREF')  
    #}else{
    #  outDir = file.path(mainDir,'plots/noNormREF')  
    #}
    #if(!dir.exists(outDir)){
    #  print('Making new Dir')
    #  dir.create(outDir,recursive = T)
    #}
    
    #### 2. Import seuratObj 
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    
    #### 3. Processing CopyKat and inferCNV results 
    #-----------------------------------------------#
    
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){

      message(paste0('Checking sample ',current_PDID))
      donorMani = projMani[projMani$PDID == current_PDID,]
      # subset srat object to keep only cells of that sample
      srat.sub = subset(srat, subset = PDID == current_PDID)
      
      #### 3. Import copyKat CNA matrix and prediction ####
      #------------- copyKat ---------------#####
      if(normREF){
        cna.mtx.fp = file.path(mainDir,'CopyKAT_output/normREF',tumourType,paste0(current_PDID,'_copykat_CNA_results.txt'))
        expr.mtx.fp = file.path(mainDir,'inferCNV_output/normREF',tumourType,current_PDID,'expr.infercnv.dat')
        gene_order_fp = file.path(mainDir,'inferCNV_output/normREF',tumourType,current_PDID,'run.final.infercnv_obj')
      }else{
        
        cna.mtx.fp = file.path(mainDir,'CopyKAT_output/noNormREF',tumourType,paste0(current_PDID,'_copykat_CNA_results.txt'))
        expr.mtx.fp = file.path(mainDir,'inferCNV_output/noNormREF',tumourType,current_PDID,'expr.infercnv.dat')
        gene_order_fp = file.path(mainDir,'inferCNV_output/noNormREF',tumourType,current_PDID,'run.final.infercnv_obj')
      }
      CNA_mat = read.delim(cna.mtx.fp)
      colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
      #colnames(CNA_mat) = gsub('\\.','-',colnames(CNA_mat))
      
      # Calculate the average CK values (expression level) for each 200kb bin across all cells within each cell type
      CNA_summary_byCellType = data.frame()
      for(celltype in celltypes){
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
      
      
      
      
      
      
      
      
      #------------- inferCNV ---------------#####
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
      for(celltype in celltypes){
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
      inferCNVsummary_byCellType$start = gene_order$start[m]
      inferCNVsummary_byCellType$stop = gene_order$stop[m]
      inferCNVsummary_byCellType = inferCNVsummary_byCellType[order(inferCNVsummary_byCellType$abspos,decreasing = F),]
      # Remove X chromosome 
      inferCNVsummary_byCellType = inferCNVsummary_byCellType[inferCNVsummary_byCellType$chrom != 23,]
      inferCNVsummary_byCellType$celltype[inferCNVsummary_byCellType$celltype == 'Tumour Cells'] = 'Tumour'
      #if(tumourType == 'Wilms'){
      #  inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c("Tumour", 'Leukocytes','Endothelium',"Muscle"))  
      #}else if(tumourType == 'Ewings'){
      #  inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c("Tumour", 'Leukocytes',"Endothelium"))
      #}else if(tumourType == 'ATRT'){
      #  inferCNVsummary_byCellType$celltype = factor(inferCNVsummary_byCellType$celltype,levels=c("Tumour", 'Leukocytes',"Endothelium","Fibroblast","Epithelium"))
      #}
      
      # Scale around zero
      #inferCNVsummary_byCellType$mean_logCN = log(inferCNVsummary_byCellType$mean_logCN)
      
      
      
      
      
      
      
      ####------------------ Generate Battenberg CN summary file ----------------####
      btb.fp = unique(donorMani$battenbergFp)
      #----- Processing Battenberg data -------#
      #dna.data = annotateBTB(btb.fp = btb.fp,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = F,method = 'totalCN')  
      #dna.data = GRanges(dna.data$Chr,IRanges(dna.data$Start,dna.data$Stop),
      #                   matNum=dna.data$matNum,patNum=dna.data$patNum, 
      #                   totCN=dna.data$tumTot,tot2min=dna.data$tot2min, idx=dna.data$idx, Chr=dna.data$Chr, 
      #                   tumFrac = dna.data$tumFrac,clonalType=dna.data$type, posID=dna.data$posID)
      dna.data = processBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID=current_PDID,tgtChrs=chrToPlot,removeBalancedSegs=F,longFormat = F,keepClonalSegs = F,method = 'totalCN')  
      
      if(keepClonalSegs){
        all_dipSegs = GRanges()
        all_altSegs = GRanges()
        for(chr in chrToPlot){
          if(! chr %in% dna.data$chr){
            next
          }
          
  
          # If there is subclonal segments within altSegs still
          segs = dna.data[(dna.data$chr == chr),]
          segs$toRemove = F
          if('sub' %in% segs$clonalType){
            segs_new = GRanges()
            # Find the major clone segment overallping with the subclonal segment
            for(k in 1:(sum(segs$clonalType == 'sub'))){
              seg.sub = segs[segs$clonalType == 'sub'][k]
              o = subsetByOverlaps(segs[segs$clonalType != 'sub'  & segs$toRemove == F],seg.sub)
              segs$toRemove[segs$idx %in% o$idx] = T
              # merge and split this into 2 to remove the clonal segments
              if(length(o)>1){
                stop('Weird...need to check!')
              }else if(length(o) == 1){
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
          
          segs = segs[segs$toRemove == F]
          mcols(segs)$toRemove = NULL
          
          
          # Extract the relevant altered copy number segments - GROUND TRUTH
          # "diploid segments" = segments with total CN == 2
          # "altered segments" = segments with total CN != 2
          dipSegs = segs[(segs$chr) & segs$totCN == 2,]
          altSegs = segs[(segs$chr == chr) & (segs$totCN != 2)]
          
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
      
      
      
      # For each cell type, calculate the average CK values per genomic bin
      # For diploid genomic segments (including copy number neutral), aggregate the values across all cells from all celltypes
      # For altered genomic segments, only aggregate across tumour cells, split by gains and losses
      # Exclude genomic bins spanning the boundary across a normal and altered regions
      # Exclude altered genomic bins which are clonal (ie. containing subclones)
      
      out_meanScore = tibble()
      # Perform segmentation: generate start and stop positions for each genomic bin
      for(i in 1:length(allSegs)){
        seg = allSegs[i,]
        chr=seg$chr
        # Matrix of copykat/inferCNV values across the chromosome - DATA
        chr.expr.CK = CNA_summary_byCellType[CNA_summary_byCellType$chrom == seg$chr,]
        chr.expr.infCNV = inferCNVsummary_byCellType[inferCNVsummary_byCellType$chrom == seg$chr,]
        
        
        nTot = seg$totCN 
        if(nTot == 2 ){
          celltype = 'all'  
        }else{
          celltype = 'Tumour'
        }
        
        stops = seq(seg$Start+binSize,seg$Stop,binSize) # if the last segment is shorter than binSize --> will automatically be removed (ie, there wont be a stop position for the last segment. eg seq(2,10,3) gives 2,5,8)
        starts = seq(seg$Start,seg$Stop - binSize,binSize)  
        for(j in 1:length(stops)){
          start = starts[j]
          stop = stops[j]
          if(celltype == 'Tumour'){
            tmp.ck = chr.expr.CK[chr.expr.CK$chrompos >= start & chr.expr.CK$chrompos<=stop & chr.expr.CK$celltype == 'Tumour',]  
            tmp.inf = chr.expr.infCNV[chr.expr.infCNV$start >= start & chr.expr.infCNV$stop <= stop & chr.expr.infCNV$celltype == 'Tumour',] 
          }else{
            tmp.ck = chr.expr.CK[chr.expr.CK$chrompos >= start & chr.expr.CK$chrompos<=stop,]  
            tmp.inf = chr.expr.infCNV[chr.expr.infCNV$start >= start & chr.expr.infCNV$stop <= stop,] 
          }
          
          if(nrow(tmp.ck) == 0){ # if no CK values in any cells are available for this segment (which is very unlikely?...)
            #message(sprintf('%s Strange situation...! %s %s %s',current_PDID,i,j,chr))
            tmp2.ck = data.frame(chr,start,stop,meanScore=NA,nTot=nTot,celltype=celltype,PDID=current_PDID,tumourType=tumourType,method='CK')
          }else{
            tmp2.ck = data.frame(chr,start,stop,meanScore = mean(tmp.ck$mean_logCN,na.rm=T),nTot=nTot,celltype=celltype,PDID=current_PDID,tumourType=tumourType,method='CK')
          }
          
          if(nrow(tmp.inf) == 0){
            #message(sprintf('%s Strange situation...!',current_PDID))
            tmp2.inf = data.frame(chr,start,stop,meanScore=NA,nTot=nTot,celltype=celltype,PDID=current_PDID,tumourType=tumourType,method='inferCNV')
          }else{
            tmp2.inf = data.frame(chr,start,stop,meanScore = mean(tmp.inf$mean_logCN,na.rm=T),nTot=nTot,celltype=celltype,PDID=current_PDID,tumourType=tumourType,method='inferCNV')
          }
          
          out_meanScore = rbind(out_meanScore,rbind(tmp2.ck,tmp2.inf))
        }
      }
      
      #write_delim(out_meanScore,file.path(mainDir,paste0('sensitivity_analysis/expressionFC_distribution/',tumourType,'_',current_PDID,'_meanScore_inferCNV_CK_',binSize,'.txt')),col_names = T,delim = '\t')
      write_delim(out_meanScore,file.path(outdir,paste0(tumourType,'_',current_PDID,'_meanScore_inferCNV_CK_',binSize,'.txt')),col_names = T,delim = '\t')
      
      out_meanScore = out_meanScore[!is.na(out_meanScore$meanScore),]
      ## Plot distribution of log_expression_fold_change for diploid segments and aneuploid segments
      #out_meanScore_ck = out_meanScore[out_meanScore$method == 'CK',]
      for(method in c('CK',"inferCNV")){
        out_meanScore_sub = out_meanScore[out_meanScore$method == method,]
        plotFun = function(noFrame=FALSE,noPlot=FALSE){
          
          par(mar=c(1.4,1.2,0.5,1),xpd=T)
          # plot densities
          xmax = max(density(out_meanScore_sub$meanScore)$x)+0.01
          xmin = min(density(out_meanScore_sub$meanScore)$x)-0.01
          ymax = max(density(out_meanScore_sub$meanScore)$y)
          d_dipSegs <- density(out_meanScore_sub$meanScore[out_meanScore_sub$celltype == 'all'])
          if(method == 'inferCNV'){
            plot(d_dipSegs,main='',axes=F,xlim=c(xmin,xmax),ylim=c(0,ymax))  
          }else if(method == 'CK'){
            plot(d_dipSegs,main='',axes=F,xlim=c(xmin,xmax))
          }
          
          polygon(d_dipSegs, col="white", border="black")
          
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
          
          #n_cells = n_distinct(gCnts_noNA$cellID)
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
          #axis(2,at=c(2000,2800),labels = c('',''),las=1,tck=-0.02,cex.axis=0.6,lwd.ticks = 0,hadj = 0.35,padj = 0.5,lwd = 0.8)
          title(sprintf('%s logFC distribution in %s - %s',method,tumourType,current_PDID),family='Helvetica',cex.main = 0.6)
          mtext(side=1,text = 'log expression FC per 5Mb',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
          mtext(side=2,text = 'Density',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
        }
        saveFig(file.path(outdir,paste0(tumourType,'_',current_PDID,'_',method,'_exprDistr')),plotFun,width = 2.4,height = 2.4,rawData = out_meanScore_sub) 
        
      }
    }
  }
}
      


### Plot aggregated distribution across all PDIDs and tumourTypes
    
# For each tumour type
# loop through each PDID

all_out_meanScore = tibble()
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('RCC','Wilms','ATRT','Ewings','NB')){
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      out_meanScore = read_delim(file.path(outdir,paste0(tumourType,'_',current_PDID,'_meanScore_inferCNV_CK_',binSize,'.txt')),col_names = T,delim = '\t')
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
    
    #n_cells = n_distinct(gCnts_noNA$cellID)
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
    #axis(2,at=c(2000,2800),labels = c('',''),las=1,tck=-0.02,cex.axis=0.6,lwd.ticks = 0,hadj = 0.35,padj = 0.5,lwd = 0.8)
    title(sprintf('%s logFC distribution in %s samples from %s tumour types ',method,n_distinct(out_meanScore_sub$PDID),n_distinct(out_meanScore_sub$tumourType)),family='Helvetica',cex.main = 0.6)
    mtext(side=1,text = 'log expression FC per 5Mb',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
    mtext(side=2,text = 'Density',family='Helvetica',font = 1,cex = 0.7,line = 0.4)
  }
  saveFig(file.path(outdir,paste0('allSamples_',method,'_exprDistr')),plotFun,width = 2.4,height = 2.4,rawData = out_meanScore_sub) 
}









    
######################################################
## Setting arbitary threshold of 0.2
mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'
thresholds = seq(0.0001,0.5,0.0005)
roc_data_by_PDID = tibble()

all_out_meanScore = all_out_meanScore[!is.na(all_out_meanScore$meanScore),]
if(sum(all_out_meanScore$meanScore == 0) > 0){
  stop('Please check!... ')
}

for(threshold in thresholds){
  all_out_meanScore$accuracy = NA
  all_out_meanScore$accuracy = ifelse((all_out_meanScore$nTot == 2) & (abs(all_out_meanScore$meanScore) < threshold), 'TN',
                                  ifelse((all_out_meanScore$nTot == 2) & (abs(all_out_meanScore$meanScore) > threshold), 'FP',
                                         
                                         ifelse((all_out_meanScore$nTot != 2) & (abs(all_out_meanScore$meanScore) < threshold),'FN',
                                                #ifelse((all_out_meanScore$nTot != 2) & (abs(all_out_meanScore$meanScore) > threshold),'TP','Weird')))))
                                                
                                                ifelse((all_out_meanScore$nTot > 2) & (all_out_meanScore$meanScore > threshold),'TP',
                                                       ifelse((all_out_meanScore$nTot > 2) & (all_out_meanScore$meanScore < -threshold),'FP', # wrong class classification
                                                              
                                                              ifelse((all_out_meanScore$nTot < 2) & (all_out_meanScore$meanScore < -threshold),'TP',
                                                                     ifelse((all_out_meanScore$nTot < 2) & (all_out_meanScore$meanScore > threshold),'FP','Weird'))))))) # wrong class classification
  if(sum(all_out_meanScore$accuracy == 'Weird')>0){
    stop('Please check weird accuracy assignment...')
  }
  tmp = all_out_meanScore %>% group_by(method,accuracy,PDID) %>% summarise(n=n())
  tmp$threshold = threshold
    
  roc_data_by_PDID = rbind(roc_data_by_PDID,tmp)
}

write_delim(roc_data_by_PDID,file.path(mainDir,'sensitivity_analysis/ROC/expressionMethod_ROC_byPDID_rawdata_2.txt'),delim = '\t',col_names = T)

roc_data_allSamples = roc_data_by_PDID %>% group_by(method,accuracy,threshold) %>% summarise(total=sum(n))
write_delim(roc_data_allSamples,file.path(mainDir,'sensitivity_analysis/ROC/expressionMethod_ROC_rawdata_2.txt'),delim = '\t',col_names = T)

roc_data_by_PDID = pivot_wider(roc_data_by_PDID,names_from = 'accuracy',values_from = 'n')
roc_data_by_PDID$FP = ifelse(is.na(roc_data_by_PDID$FP),0,roc_data_by_PDID$FP)
roc_data_by_PDID$FN = ifelse(is.na(roc_data_by_PDID$FN),0,roc_data_by_PDID$FN)
roc_data_by_PDID$TP = ifelse(is.na(roc_data_by_PDID$TP),0,roc_data_by_PDID$TP)
roc_data_by_PDID$TN = ifelse(is.na(roc_data_by_PDID$TN),0,roc_data_by_PDID$TN)
# calculate y-axis True Positive Rate (TPR) = TP/(TP+FN)
roc_data_by_PDID$TPR = roc_data_by_PDID$TP/(roc_data_by_PDID$TP + roc_data_by_PDID$FN)
# calculate x-axis False Positive Rate (FPR) = FP/(FP+TN)
roc_data_by_PDID$FPR = roc_data_by_PDID$FP/(roc_data_by_PDID$FP + roc_data_by_PDID$TN)
write_delim(roc_data_by_PDID,file.path(mainDir,'sensitivity_analysis/ROC/expressionMethod_ROC_byPDID_dataForPlot_2.txt'),delim = '\t',col_names = T)

for(method in c('CK','inferCNV')){
  dat = roc_data_by_PDID[roc_data_by_PDID$method == method,]
  p = ggplot(dat,aes(x=FPR,y=TPR))+
    geom_point()+
    geom_line() +
    theme_bw()+
    ggtitle(method)
  
  pdf(file.path(mainDir,paste0('sensitivity_analysis/',method,'_ROC.pdf')))
  print(p)
  dev.off()
}
    

####### A different way to calculate ROC and AUC
library(glmnet)
library(pROC)





plot_modelAUC = function(model.results,model.names, selected.features, style='sep',
                         my.col=NULL,bg.col=grey(0.8), cex=1.5,line=NULL, auc=T){
  
  if(length(model.results) != length(model.names)){
    stop('Model.names and Model.results are not of the same length')
  }
  
  if(is.null(my.col)){
    if (length(model.results) < 13){
      require(RColorBrewer)
      my.col =  brewer.pal(length(model.results),'Dark2')
    }
  }
  #my.col = c(my.col[2],my.col[c(1,3:5)])
  if(is.null(line)){
    lines = seq(-20,(-20-3*length(model.results)+3*1),by=-2)
  }
  
  if (!style %in% c('sep','avg')) stop("style must be 'sep' or 'avg'")
  par(mar=c(5,10,5,10))
  plot(x=c(0,1),y=c(0,1),type='n',xaxs='i',yaxs='i',
       xlab="False Positive Rate",ylab="True Positive Rate", cex.axis = 1.7, cex.lab=1.9)
  grid.position <- seq(0.1, 0.9, 0.1)
  segments(x0=grid.position,y0=0,x1=grid.position,y1=1,lty=2,col=bg.col)
  segments(x0=0,y0=grid.position,x1=1,y1=grid.position,lty=2,col=bg.col)
  
  for(i in 1:length(model.results)){
    
    result = model.results[[i]]
    fg.col = my.col[i]
    
    if(model.names[i] == 'RF'){
      pred=list()
      for (j in 1:length(result$pred)){
        pred[[j]] = result$pred[[j]]$`1`
      }
      perf <- evaluate.classification(pred, result$truth)
    }else if(model.names[i] == 'SVM'){
      pred=list()
      for (j in 1:length(result$pred)){
        pred[[j]] = as.data.frame(attr(result$pred[[j]],"probabilities"))$`1`
      }
      perf <- evaluate.classification(pred, result$truth)
    }else{
      perf <- evaluate.classification(result$pred, result$truth)  
    }
    
    
    switch(EXPR=style,
           sep=plot(perf$roc,avg="none",spread.estimate="stderror",xaxs='i',yaxs='i',
                    xlab="False Positive Rate",ylab="True Positive Rate",col=fg.col,add=T),
           avg=plot(perf$roc,avg="vertical",spread.estimate="stderror",xaxs='i',yaxs='i',
                    xlab="False Positive Rate",ylab="True Positive Rate",col=fg.col,add=T)
    )
    
    if(auc){
      mtext(paste('AUC =',format(mean(perf$auc), digits=3), '±', format(sd(perf$auc), digits=3),
                  '  Model: ',model.names[i]),
            side=3, adj=1, line=lines[i], cex=cex, col=fg.col)}
    
  }
  abline(a=0,b=1,lty=1,col=bg.col)
  
  title(main=paste(c('10-fold cross validation', paste(selected.features, collapse = ' + ')),
                   collapse = '\n'), cex.main = 1.5)
}











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
binSize=25e6
threshold = 0.2
main.srat = nb.srat
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
      if(!chr %in% unique(segs$chr)){
        next
      }
      chrsegs = segs[segs$chr == chr,]
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
          #if(i == length(stops) - 1){
          #  segLen = start - stop +1
          #  segLen_next = start[i+1] - stop[i+1] +1
          #  if(segLen_next)
          #}
          
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

out$segLen = out$stop - out$start 
out = out[out$segLen > 20e6,]
d = out %>% group_by(PDID,celltype,cna,predCall) %>% tally()

d$truth = ifelse(d$celltype == 'Tumour',as.character(d$cna),FALSE)

d$result = ifelse(d$truth == TRUE & d$predCall == 'yes','TP',
                  ifelse(d$truth == TRUE & d$predCall == 'no','FN',
                         ifelse(d$truth == FALSE & d$predCall == 'yes','FP',
                                ifelse(d$truth == FALSE & d$predCall == 'no','TN','?'))))



saveRDS(list(out,d),'/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_nb.CK.Quant_v2.RDS')


## Calculate Sensitivity and Specificity #####
nb = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_nb.CK.Quant_v2.RDS')
rcc = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v2_rcc.CK.Quant.RDS')

nb[[2]]$dataset = 'NB'
rcc[[2]]$dataset = 'RCC'

data = rbind(nb[[2]],rcc[[2]])

data2 = data %>% group_by(dataset,result) %>% summarise(total = sum(n,na.rm = T))
data3 = data2 %>% filter(!is.na(result)) %>% group_by(dataset) %>% summarise(sensitivity = 100*total[result=='TP']/(total[result=='TP']+total[result=='FN']),
                                                           specificity = 100*total[result=='TN']/(total[result=='FP']+total[result=='TN']))
data3 = pivot_longer(data3,cols = c(2:3),names_to = 'type',values_to = 'value')

ggplot(data3,aes(dataset,y=value,fill=type))+
  geom_col(position = 'dodge')+
  scale_fill_manual(values = c('black','grey'))+
  theme_bw()



