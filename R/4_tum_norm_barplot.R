# This scripts plot the tum_norm barplot and UMAPs based on the 
# tum-norm calls by alelleIntegrator and copyKat

setwd('~/lustre_mt22/CN_methods/')


#############
# Libraries #
#############
library(tidyverse)
library(readxl)
library(plotrix)
#########################
# Set Global parameters #
#########################
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204_v2/'


for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('RCC','NB','Wilms','ATRT','Ewings')){
    srat = readRDS(file.path(mainDir,'sc_seuratObjects',tumourType,paste0(tumourType,'_ann.RDS')))
    m = match(rownames(srat@meta.data),rownames(srat@reductions$umap@cell.embeddings))
    if(sum(is.na(m))>0){
      stop("Unidentified cells....")
    }
    # Extract CK output
    dd = cbind(srat@meta.data,srat@reductions$umap@cell.embeddings[m,])
    if('AIcall_normREF_v2' %in% colnames(dd)){
      dd$AIcall_normREF = as.character(dd$AIcall_normREF_v2)
    }
    # Define relevant cell types to keep
    if(tumourType %in% c('ATRT','Ewings','Wilms')){
      celltype_toKeep = c('Endothelium','Leukocytes','Tumour')
    }else if(tumourType == 'RCC'){
      celltype_toKeep = c('PTC','Leukocytes','Tumour')
    }else if(tumourType == 'NB'){
      celltype_toKeep = c('Endothelium','Mesenchyme','Leukocytes','Tumour')
    }
    
    dd$AIcall_normREF = ifelse(is.na(dd$AIcall_normREF),'Uncalled',
                               ifelse(dd$AIcall_normREF == 'abbFrac','Tumour',
                                      ifelse(dd$AIcall_normREF == 'normFrac','Normal','Uncalled')))
    dd$CKpred.normREF.default.100perc.2397 = ifelse(is.na(dd$CKpred.normREF.default.100perc.2397),'Uncalled',
                                                    ifelse(dd$CKpred.normREF.default.100perc.2397 == 'aneuploid','Aneuploid',
                                                           ifelse(dd$CKpred.normREF.default.100perc.2397 == 'diploid','Diploid','Uncalled')))
    dd$PDID = as.factor(dd$PDID)
    dd.sub = dd[dd$annot %in% celltype_toKeep,]
    dd.sub$annot = factor(dd.sub$annot, levels = celltype_toKeep)
    
    
    for(method in c('AI','CK')){
      if(method == 'AI'){
        col = 'AIcall_normREF'
        cols = c(Normal = '#c8c8c8',
                 Tumour = '#b02c46',
                 Uncalled = '#474646')  
        cols_umap = c(Normal = 'black',
                 Tumour = '#b02c46',
                 Uncalled = '#dedede')
        
        dd.sub = dd.sub[dd.sub$AIcall_normREF != 'Uncalled',]
        dd.sub$AIcall_normREF = factor(dd.sub$AIcall_normREF,levels = c('Tumour','Normal'))
        
      }else if(method == 'CK'){
        col = 'CKpred.normREF.default.100perc.2397'
        cols = c(Diploid = '#c8c8c8',
                 Aneuploid = '#b02c46',
                 Uncalled = '#474646')
        
        cols_umap = c(Diploid = 'black',
                 Aneuploid = '#b02c46',
                 Uncalled = '#dedede')
        
        dd.sub = dd.sub[dd.sub$CKpred.normREF.default.100perc.2397 != 'Uncalled',]
        dd.sub$CKpred.normREF.default.100perc.2397 = factor(dd.sub$CKpred.normREF.default.100perc.2397,levels = c('Aneuploid','Diploid'))
      }
      #Define the layout
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
        celltype = levels(dd.sub$annot)[1]
        for(celltype in levels(dd.sub$annot)){
          par(mar=c(0.05,1,0.2,0.1),xpd=TRUE)
          #par(mar=c(0.2,1.5,0.5,0.1),xpd=TRUE)
          
          tmp=as.matrix(table(dd.sub[dd.sub$annot == celltype,col],dd.sub$PDID[dd.sub$annot == celltype]))
          tmp = sweep(tmp,2,colSums(tmp),'/')
          
          if(celltype == levels(dd.sub$annot)[length(levels(dd.sub$annot))]){
            par(mar=c(2.5,1,0.5,0.1),xpd=TRUE)
            barplot(tmp,
                    col=cols[rownames(tmp)],
                    space=0.12,axes = FALSE,
                    las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
            text(x = seq(0.62,5.7,by = 1.13),y = -0.06,colnames(tmp),cex=0.6,family = 'Helvetica',font=1,srt=90,adj = 1)
            #text(x =-1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
            text(x = -0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
            text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$annot == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
          }else{
            barplot(tmp,
                    col=cols[rownames(tmp)],
                    space=0.12,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                    las = 1,main = '',border = F)
            text(x =-0.55,y = 0.5,celltype,cex=0.6,family = 'Helvetica',font=1,srt=90)
            text(x = -0.1,y = 0.5,paste0('n=',sum(dd.sub$annot == celltype)),cex=0.45,family = 'Helvetica',font=1,srt=90)
          }
          
        }
        
      }
      saveFig(file.path(mainDir,paste0('bar_plots/',tumourType,'_barplot_',method,'_noUncalled')),plotFun,width = (0.17+0.15*n_distinct(dd.sub$PDID)),height = 2.75,rawData = dd.sub,res=500)
      
      
      
      #####=========== plot UMAPs ===================####
      
      plotFun = function(noFrame=FALSE,noPlot=FALSE){
        layout(matrix(c(1,1),ncol=1))
        par(mar=c(0.6,0.6,0.8,0.1))
        
        plot(dd$UMAP_1,dd$UMAP_2,
             las=1,
             #type='n',
             #xlim=c(min(dd$UMAP_1),max(min(dd$UMAP_1))),
             #ylim=c(min(dd$UMAP_2),max(min(dd$UMAP_2))), 
             main=ifelse(noFrame,'',paste0(tumourType,'-',method)),
             cex.main = 0.45,xaxt='n',yaxt='n',
             col = cols_umap[dd[[col]]],
             #col = cols_umap[dd$AIcall_normREF_v2],
             pch = 19,cex=0.009,
             frame.plot=T)
        
        mtext('UMAP 1',side=1,line = -0.3,family = 'Helvetica',cex = 0.35)
        mtext('UMAP 2',side=2,line = 0.1,family = 'Helvetica',cex = 0.35)
      }
      saveFig(file.path(mainDir,'umap',paste0(tumourType,'_',method,'_umap')),plotFun,rawData=dd,width = 2,height = 2.1,res = 500)
    }
    
    
    
    ### Plot UMAP and colored by cell types
    plotFun = function(noFrame=FALSE,noPlot=FALSE){
      layout(matrix(c(1,1),ncol=1))
      par(mar=c(0.6,0.6,0.8,0.1))
      
      # Define color by PDID
      if(length(unique(dd$PDID)) > 4){
        colfunc <- colorRampPalette(c("#383838", "#E2E2E2"))  
      }else{
        colfunc <- colorRampPalette(c("#383838", "#C5C5C5"))  
      }
      
      ccs = colfunc(length(unique(dd$PDID)))
      names(ccs) = unique(dd$PDID)
      
      # Plot UMAP
      plot(dd$UMAP_1,dd$UMAP_2,
           las=1,
           #type='n',
           #xlim=c(min(dd$UMAP_1),max(min(dd$UMAP_1))),
           #ylim=c(min(dd$UMAP_2),max(min(dd$UMAP_2))), 
           main=ifelse(noFrame,'',tumourType),
           cex.main = 0.45,xaxt='n',yaxt='n',
           col = ccs[dd$PDID],
           pch = 19,cex=0.009,
           frame.plot=T)
      
      mtext('UMAP 1',side=1,line = -0.25,family = 'Helvetica',cex = 0.45)
      mtext('UMAP 2',side=2,line = 0.1,family = 'Helvetica',cex = 0.45)
      
      
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ annot,data=dd,FUN=mean)
      
      #Position tweaks
      
      #mids[mids$annot=='Leukocyte','UMAP_2'] = mids[mids$finalAnn=='Leukocyte','UMAP_2'] + 4.9
      #mids[mids$finalAnn=='Leukocyte','UMAP_1'] = mids[mids$finalAnn=='Leukocyte','UMAP_1'] - 0.3
      #mids[mids$finalAnn=='PTC','UMAP_2'] = mids[mids$finalAnn=='PTC','UMAP_2'] - 4.9
      #mids[mids$finalAnn=='PTC','UMAP_1'] = mids[mids$finalAnn=='PTC','UMAP_1'] -0.8
      #mids[mids$finalAnn=='Tumour','UMAP_2'] = mids[mids$finalAnn=='Tumour','UMAP_2'] + 3
      #mids[mids$finalAnn=='Tumour','UMAP_1'] = mids[mids$finalAnn=='Tumour','UMAP_1'] + 0.9
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$annot,cex = 0.3,xpad = 1.2,ypad = 1.9,border = T,
                   col=c(rep('black',sum(mids$annot != 'Tumour')),'#b02c46'),font=c(rep(1,sum(mids$annot != 'Tumour')),2))
      
      
      legend(x=3.5, y=-1,legend=names(ccs),fill = ccs,lwd = 0,cex = 0.3,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    }
    saveFig(file.path(mainDir,'umap',paste0(tumourType,'_celltype_umap')),plotFun,rawData=dd,width = 2,height = 2.1,res = 500)
  }
}
    
    



  