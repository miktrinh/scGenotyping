# Plot correlation between inferCNV and copyKat

setwd('~/lustre_mt22/CN_methods/revision_2204_v2/')

#############
# Libraries #
#############
library(tidyverse)
library(readxl)
source('~/lustre_mt22/CN_methods/scripts/finalScripts/misc.R')

#########################
# Set Global parameters #
#########################
normREF = T
tgtChr = c(1:22)
# Import Manifest
projMani = read_excel("~/lustre_mt22/projectManifest.xlsx",sheet = "alleleIntegrator")
# Remove RCC normal samples from projMani
projMani = projMani[!grepl('^n_',projMani$SampleID),]  
mainDir = '~/lustre_mt22/CN_methods/revision_2204/'
# Set outdir
outdir = file.path('~/lustre_mt22/CN_methods/revision_2204_v2/sensitivity_analysis/')

if(!dir.exists(outdir)){
  print('Making new Dir')
  dir.create(outdir,recursive = T)
}



# For each tumour type
# loop through each PDID
corr_data = tibble()
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('NB','RCC','ATRT','Wilms','Ewings')){
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    
    for(current_PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Collecting inferCNV and CopyKat values for Sample %s - tumourType: %s',current_PDID,tumourType))
      donorMani = projMani[projMani$PDID == current_PDID,]
      srat.sub = subset(srat, subset = PDID == current_PDID)
      
      ###------------ Import inferCNV results ------------###
      #### 1. Import inferCNV output, ie. final heatmap matrix ####
      if(normREF){
        expr.mtx.fp = file.path(mainDir,'inferCNV_output/normREF',tumourType,current_PDID,'expr.infercnv.dat')
        gene_order_fp = file.path(mainDir,'inferCNV_output/normREF',tumourType,current_PDID,'run.final.infercnv_obj')
      }else{
        expr.mtx.fp = file.path(mainDir,'inferCNV_output/noNormREF',tumourType,current_PDID,'expr.infercnv.dat')
        gene_order_fp = file.path(mainDir,'inferCNV_output/noNormREF',tumourType,current_PDID,'run.final.infercnv_obj')
      }
      expr.mtx = read.delim(expr.mtx.fp,sep = '\t')
      # Import gene_order
      gene_order = readRDS(gene_order_fp)
      gene_order = gene_order@gene_order
      gene_order$chrom = gsub('chr','',gene_order$chr)
      gene_order = gene_order[gene_order$chrom %in% tgtChr,]
      # Subset expr.mtx to keep only genes on Chromosome of interest
      expr.mtx = expr.mtx[rownames(expr.mtx) %in% rownames(gene_order),]
      
      
      #### Calculate the mean expression shifts for each gene across all cells belonging to a given celltype
      inferCNVsummary_byCellType = data.frame()
      for(celltype in unique(srat.sub$annot)){
        CNA_mat_sub = expr.mtx[,c(which(gsub('^X','',colnames(expr.mtx)) %in% gsub('-','.',rownames(srat.sub@meta.data[srat.sub@meta.data$annot == celltype,]))))]
        chrom_tmp = data.frame(celltype = celltype,gene = rownames(CNA_mat_sub),mean_logCN = apply(CNA_mat_sub,MARGIN = 1,FUN = function(x){mean(log(x))}))
        inferCNVsummary_byCellType = rbind(inferCNVsummary_byCellType,chrom_tmp)
      }
      
      #inferCNVsummary_byCellType$mean_logCN = log(inferCNVsummary_byCellType$mean_logCN)
      
      
      
      
      #### 2. Import copyKat CNA matrix (per gene values) and prediction ####
      if(normREF){
        cna.mtx.fp = file.path(mainDir,'CopyKAT_output/normREF',tumourType,paste0(current_PDID,'_copykat_CNA_raw_results_gene_by_cell.txt'))
      }else{
        cna.mtx.fp = file.path(mainDir,'CopyKAT_output/noNormREF',tumourType,paste0(current_PDID,'_copykat_CNA_raw_results_gene_by_cell.txt'))
      }
      
      CNA_mat = read.delim(cna.mtx.fp)
      colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
      # Subset CNA_mat to keep only genes on Chromosome of interest
      CNA_mat = CNA_mat[CNA_mat$chromosome_name %in% tgtChr,]
      
      # Genes that are common in both inferCNV and CopyKat
      genesToKeep = intersect(CNA_mat$hgnc_symbol,rownames(expr.mtx))
      #sum(!CNA_mat$hgnc_symbol %in% rownames(expr.mtx))
      #sum(!rownames(expr.mtx) %in% CNA_mat$hgnc_symbol)
      
      CNA_summary_byCellType = data.frame()
      # subset by annotated cell type
      for(celltype in unique(srat.sub$annot)){
        CNA_mat_sub = CNA_mat[,c(1:7,which(colnames(CNA_mat) %in% gsub('-','.',rownames(srat.sub@meta.data[srat.sub@meta.data$annot == celltype,]))))]
        
        if(ncol(CNA_mat_sub) == 8){ # if there is only 1 cell
          chrom_tmp=data.frame(celltype = celltype,CNA_mat_sub)
          colnames(chrom_tmp)[9] = 'mean_logCN'
        }else if (ncol(CNA_mat_sub) > 8){
          chrom_tmp = data.frame(celltype = celltype,CNA_mat_sub[,c(1:7)],mean_logCN = apply(CNA_mat_sub[,-c(1:7)],MARGIN = 1,FUN = mean))
        }else if (ncol(CNA_mat_sub) < 8){
          chrom_tmp = data.frame(celltype = celltype,CNA_mat_sub[,c(1:7)],mean_logCN = NA)
        }
        CNA_summary_byCellType = rbind(CNA_summary_byCellType,chrom_tmp)
      }
      
      
      
      ###----------------------------------------------------------------#####
      ### Merge inferCNV and copyKAT average matrices ####
      inferCNVsummary_byCellType$celltype_gene = paste0(inferCNVsummary_byCellType$celltype,'_',inferCNVsummary_byCellType$gene)
      colnames(inferCNVsummary_byCellType)[3] = 'inferCNV_mean_logCN'
      CNA_summary_byCellType$celltype_gene = paste0(CNA_summary_byCellType$celltype,'_',CNA_summary_byCellType$hgnc_symbol)
      colnames(CNA_summary_byCellType)[9] = 'ck_mean_logCN'
      
      data = merge(inferCNVsummary_byCellType[,c(3,4)],CNA_summary_byCellType[,c(9,10)],by = 'celltype_gene')
      data$PDID = current_PDID
      corr_data = rbind(corr_data,data)
    }
  }
}

write_delim(corr_data,file.path(mainDir,'../revision_2204_v2/sensitivity_analysis/inferCNV_CK_corr_data.txt'),delim = '\t',col_names = T)

corr_data = read.delim(file.path(mainDir,'../revision_2204_v2/sensitivity_analysis/inferCNV_CK_corr_data.txt'))
library(stringr)
corr_data$celltype = sapply(strsplit(corr_data$celltype_gene,split = '_'),'[',1)
corr_data$celltype2 = ifelse(corr_data$celltype == 'Tumour','Tumour','Others')

data = corr_data %>% filter(!is.na(inferCNV_mean_logCN) & (!is.na(ck_mean_logCN)))
data = data[data$celltype != 'Leukocytes',]
corr_coeff = cor(data$inferCNV_mean_logCN,data$ck_mean_logCN)
#idx = sample(1:nrow(data),size = nrow(data)*0.2)
idx=c(1:nrow(data))
plotFun = function(noFrame=TRUE,noPlot=FALSE){
  par(mar=c(3,3,2.5,2.5),xpd=F)  
  noFrame=T
  y = data$inferCNV_mean_logCN[idx]
  x = data$ck_mean_logCN[idx]
  n_min = round(min(c(y,x)),2)
  n_max = round(max(c(y,x)),2)
  
  plot(x=x,y=y,
       #xlim = c(n_min,n_max),ylim=c(n_min,n_max),
       xlim = c(-0.05,0.05),ylim=c(-0.05,0.05),
       las=1,
       type='n',axes = F,
       xlab=ifelse(noFrame,'','CopyKAT'),
       ylab=ifelse(noFrame,'','inferCNV'),
       frame.plot=F,main = paste0('R = ',round(corr_coeff,3)),cex.main=0.8)
  
  points(x=x,y=y,pch=19,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=0.01)
  
  axis(side = 2,at =seq(-0.4,0.3,0.1),labels = round(seq(-0.4,0.3,0.1),1),tck=-0.015,lwd = 0.8,cex.axis=0.6,las=1,lwd.ticks = 0.5,hadj = 0.4)
  axis(side = 1,at =seq(-0.4,0.3,0.1),labels = round(seq(-0.4,0.3,0.1),1),tck=-0.015,lwd = 0.8,cex.axis=0.6,las=0,lwd.ticks = 0.5,padj=-2.5)
  
  mtext(side=1,text = 'CopyKAT',family='sans',font = 1,cex = 0.6,line = 1.7)
  mtext(side=2,text = 'inferCNV',family='sans',font = 1,cex = 0.6,line = 1.7)
}


saveFig(file.path(outdir,paste0('FigS6_inferCNV_CK_correlation_sample02_v3tmp')),plotFun,width = 5.0,height = 5,res=500,rawData = data)  














p=ggplot(corr_data,aes(x=ck_mean_logCN,y=inferCNV_mean_logCN))+
  #geom_point(size=0.5,alpha=0.1)+
  #geom_smooth(method = "lm", formula = y ~ x) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
  #facet_wrap(vars(celltype2))+
  xlim(-0.3,0.3) + ylim(-0.3,0.3)
  #geom_abline(slope = 1,intercept = 0)
p
a = lm(inferCNV_mean_logCN~ck_mean_logCN,data = corr_data)

eq <- as.expression(substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(a)[1], digits = 4),
                      b = format(coef(a)[2], digits = 4),
                      r2 = format(summary(a)$r.squared, digits = 3))))

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

dftext <- data.frame(x = 0.1, y = 0.25, eq = as.character(as.expression(eq)))

p + geom_text(x=0.2,y=0.25,label = eq, parse = TRUE)





