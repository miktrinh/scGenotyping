#### This scripts stores the old function and code in generatePlots.R 
# most of these I suspect are no longer relevant. However, just in case....

# Import RCC and NB datasets
rcc.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')
nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')
chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')



# Figure 2C - Example of CopyKat output for PTC and tumour RCC cells ####
sampleList = 'PD37228'
rccCK = function(sampleList='all'){
  rcc.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')
  rcc.srat$finalAnn = rcc.srat$finalAnn2
  
  # Import chromInfo
  chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')
  # Import Manifest
  projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = "RCC_mani")
  # CopyKat results
  copykat.results = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_rcc.CK.normREF.default.2397/v5_CKresults_default_normREF_80perc_2397.rds')
  
  if(sampleList != 'all'){
    rcc.srat = subset(rcc.srat,subset = PDID %in% sampleList)
  }
  # Extract data for PD36793 only (as examples)
  #----- Processing CopyKat results -------#
  for(i in 1:length(copykat.results)){
    # Get copyKat CNA matrix and prediction
    CNA_summary_byCellType = data.frame()
    CNA_mat = copykat.results[[i]]$CNAmat
    colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
    colnames(CNA_mat) = gsub('\\.','-',colnames(CNA_mat))
    pred = as.data.frame(copykat.results[[i]]$prediction)
    
    sample = unique(rcc.srat@meta.data[rownames(rcc.srat@meta.data) %in% rownames(pred),]$PDID)
    PDID=sample
    
    if(length(sample) > 1){
      message(paste0('More than 1 sample detected: i=',i,', samples are ',sample))
    }else if(length(sample) <1){
      next
    }else if(length(sample)==1){
      message(paste0('Checking sample ',sample))
    }
    
    # subset annrcc.srat object to keep only cells of that sample
    srat = subset(rcc.srat, subset = PDID == sample)
    
    
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
    
    
    # Remove X chromosome 
    CNA_summary_byCellType = CNA_summary_byCellType[CNA_summary_byCellType$chrom != 23,]
    CNA_summary_byCellType$celltype = factor(CNA_summary_byCellType$celltype,levels=c('Tumour','PTC','Leukocytes'))
    if(sampleList != 'all'){
      CNA_summary_byCellType = CNA_summary_byCellType[CNA_summary_byCellType$celltype != 'Leukocytes',]
      CNA_summary_byCellType$celltype = factor(CNA_summary_byCellType$celltype,levels=c('Tumour','PTC'))
    }
    
    
    ####------------------ Generate Battenberg CN summary file ----------------####
    # Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
    donorMani = projMani[projMani$PDID == sample,]
    btb.fp = unique(donorMani$battenbergFp[!grepl('^n_',donorMani$SampleID)])
    #----- Processing Battenberg data -------#
    dna.data = annotateBTB(btb.fp,subCl.minSegLen = 1e7,PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'totalCN')  
    # Remove chromosome X
    dna.data = dna.data[dna.data$Chr != 23,]
    
    # Add equivalent lines for other cell types
    dna.data$celltype = 'Tumour'
    tmp = rbind(data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype='PTC'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype='PTC'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype='Leukocytes'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype='Leukocytes'))
    
    dna.data = rbind(dna.data,tmp)
    dna.data$log_CNratio = log(as.numeric(dna.data$tumTot)/2)
    
    # Separate subclone and major clone CN profile
    subCl.dna = dna.data[dna.data$type == 'sub',]
    majCl.dna = dna.data[dna.data$type == 'maj',]
    
    
    
    #----- Plotting copykat results! -------#
    chromInfo2 = chromInfo[chromInfo$chrom != 23,]
    
    plotFun = function(noFrame=T,noPlot=FALSE){
      noFrame=T
      layout(mat=matrix(c(1:nlevels(CNA_summary_byCellType$celltype)),ncol=1),
             heights = rep(2,nlevels(CNA_summary_byCellType$celltype)))
      
      for(celltype in levels(CNA_summary_byCellType$celltype)){
        
        if(sampleList == 'all'){
          par(mar=c(0.2,0.9,0.2,0.4),xpd=TRUE)
          ylim = c(round(min(CNA_summary_byCellType$mean_logCN)-0.1,digits = 1),round(max(CNA_summary_byCellType$mean_logCN)+0.2,digits = 1))
          text.cex = 0.7
          ybottom=min(round(min(CNA_summary_byCellType$mean_logCN)-0.1,digits = 1),round(log(0.5)/2,2))
          
          if((round(max(dna.data$log_CNratio),2) - round(max(CNA_summary_byCellType$mean_logCN)+0.1,digits = 1))>0.2){
            ytop=max(round(max(CNA_summary_byCellType$mean_logCN)+0.1,digits = 1),round(max(dna.data$log_CNratio)/2,2))
            type=2
          }else{
            ytop=max(round(max(CNA_summary_byCellType$mean_logCN)+0.1,digits = 1),round(max(dna.data$log_CNratio),2))
            type=1
          }
          ytext = ytop + 0.1
          chrom.y = 0.85
        }else{
          par(mar=c(0.2,0.9,0.2,0.4),xpd=TRUE)
          ylim = c(round(min(CNA_summary_byCellType$mean_logCN)-0.1,digits = 1),round(max(CNA_summary_byCellType$mean_logCN)+0.2,digits = 1))
          text.cex = 0.7
          ybottom=min(round(min(CNA_summary_byCellType$mean_logCN)-0.1,digits = 1),round(log(0.5)/2,2))
          ytop=max(round(max(CNA_summary_byCellType$mean_logCN)+0.1,digits = 1),round(max(dna.data$log_CNratio),2))
          ytext = ytop + 0.1
          chrom.y = 0.7
        }
        
        tmp = CNA_summary_byCellType[CNA_summary_byCellType$celltype == celltype,]
        tmp = tmp[order(tmp$abspos,decreasing = F),]
        dna = majCl.dna[majCl.dna$celltype == celltype,]
        ncells = nrow(rcc.srat@meta.data[rcc.srat@meta.data$PDID == PDID & rcc.srat@meta.data$finalAnn == celltype,])
        
        # Plot main frame
        plot(CNA_summary_byCellType$abspos, CNA_summary_byCellType$mean_logCN,
             las=1,
             type='n',
             #xlim=c(-15.5,15),
             ylim=ylim,
             xlab=ifelse(noFrame,'','Genomic Position'),
             ylab=ifelse(noFrame,'',''),
             #main=ifelse(noFrame,'',''),
             xaxt=ifelse(noFrame,'n','s'),
             yaxt=ifelse(noFrame,'n','s'),
             frame.plot=F)
        
        
        #text(x=1e3,y=ytext,celltype,cex=text.cex[1],family = 'Helvetica',font=2,adj = 0)
        #text(x=2.86e9,y=ytext,paste0('n=',ncells),cex=text.cex[2],family = 'Helvetica',font=1,adj = 1)
        #axis(2,at=c(0),labels = c(0),las=1,pos = 0,tck = -.02,lwd = 0.3,cex.axis=0.65,hadj = -0.8,padj = 0.5)
        if(sampleList=='all'){
          #axis(2,at=c(0),labels = c(0),las=1,pos = 0,tck = -.02,lwd = 0.3,cex.axis=0.75,hadj = -0.8,padj = 0.5)
          axis(2,at=c(ybottom+0.05,0,ytop-0.05),labels = c('Low',0,'High'),las=1,pos = 0,tck = -.00,lwd = 0.3,cex.axis=0.6,hadj = 0.3,padj = 0.5)
          #if((type==1) & (round(log(0.5),2) < ybottom)){
          #  axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.4,cex.axis=0.6,hadj = 2.1,
          #       at=c(0,round(log(3/2),2),round(log(4/2),2)),col.axis = '#b02c46',
          #       labels = c(2,3,4))    
          #}else if((type==1) & (round(log(0.5),2) > ybottom)){
          #  axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.4,cex.axis=0.6,hadj = 2.1,
          #       at=c(round(log(0.5),2),0,round(log(3/2),2),round(log(4/2),2)),col.axis = '#b02c46',
          #       labels = c(1,2,3,4))    
          #}else if((type==2) & (round(log(0.5)/2,2) < ybottom)){
          #  axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.4,cex.axis=0.6,hadj = 2.1,
          #       at=c(0,round(log(3/2)/2,2),round(log(4/2)/2,2)),col.axis = '#b02c46',
          #       labels = c(2,3,4))
          #}else if((type==2) & (round(log(0.5)/2,2) > ybottom)){
          #  axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.4,cex.axis=0.6,hadj = 2.1,
          #       at=c(round(log(0.5)/2,2),0,round(log(3/2)/2,2),round(log(4/2)/2,2)),col.axis = '#b02c46',
          #       labels = c(1,2,3,4))
          axis4 = c(round(log(0.5)/2,2),0,round(log(3/2)/2,2),round(log(4/2)/2,2))
          names(axis4) = c(1,2,3,4)
          axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.4,cex.axis=0.6,hadj = 2.1,
               at=axis4[axis4>ybottom & ybottom < ytop],col.axis = '#b02c46',
               labels = names(axis4[axis4>ybottom & ybottom < ytop]))    
          
        }else{
          axis(2,at=c(ybottom+0.05,0,ytop-0.05),labels = c('Low',0,'High'),las=1,pos = 0,tck = -.00,lwd = 0.3,cex.axis=0.6,hadj = 0.1,padj = 0.5)
          axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.3,cex.axis=0.6,hadj = 1.9,col.axis = '#b02c46',
               at=c(round(log(0.5)/2,2),0,round(log(3/2)/2,2)),
               labels = c(1,2,3))
        }
        
        
        
        #Plot background chromosome
        xleft = c(0,chromInfo2[chromInfo2$arm == 'q' & chromInfo2$chrom!=22,]$abspos*1000)
        xright = c(chromInfo2[chromInfo2$arm == 'q',]$abspos*1000)
        
        text(x=1.2e9,y=ytext,paste0(celltype,'_',PDID,' (n=',ncells,')'),cex=text.cex[1],family = 'Helvetica',font=1,adj = 0)
        
        col = replicate(c('white','#dbdbdb'),n = 22/2)
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
             border = 'black',lwd = 0.5)
        #textSize = (xright-xleft)/xright[1]
        #text(x=(xleft+xright)/2,y = 0.9,labels = c(1:22),cex = c(1.5*textSize),font = 2)
        
        #text(x=(xleft+xright)/2,y = chrom.y,labels = c(1:22),cex = c(rep(0.7,10),rep(0.62,4),rep(0.5,4),rep(0.31,4)),font = 1)
        
        # Plot ground truth
        #segments(x0=min(xleft),x1 = max(xright), 
        #         y0=0, y1=0,
        #         col = 'black')
        rect(xleft = dna[dna$posType=='Start',]$abspos_kb*1000,
             xright = dna[dna$posType=='Stop',]$abspos_kb*1000,
             ybottom = 0,
             ytop=dna[dna$posType=='Start',]$log_CNratio/2,
             col=colAlpha('#b02c46',1),
             border=colAlpha('#b02c46',1),lwd = 1.0)
        
        # Plot subclone total CN
        if(celltype == 'Tumour' & nrow(subCl.dna) > 0){
          for(chr in unique(subCl.dna$Chr)){
            if(unique(subCl.dna[subCl.dna$Chr == chr,]$log_CNratio == 0)){
              lines(subCl.dna[subCl.dna$Chr == chr,]$abspos_kb*1000,
                    subCl.dna[subCl.dna$Chr == chr,]$log_CNratio/2,col='#4169E1',lwd=1.5)
            }else{
              rect(xleft = subCl.dna[subCl.dna$Chr == chr &subCl.dna$posType=='Start',]$abspos_kb*1000,
                   xright = subCl.dna[subCl.dna$Chr == chr &subCl.dna$posType=='Stop',]$abspos_kb*1000,
                   ybottom = 0,
                   ytop=subCl.dna[subCl.dna$Chr == chr & subCl.dna$posType=='Start',]$log_CNratio/2,
                   col=colAlpha('#4169E1',0.7),
                   border=colAlpha('#4169E1',0.7),lwd = 0.1)
            }
            
          }
          
        }
        
        
        # Plot CopyKat output
        
        lines(x=tmp$abspos,tmp$mean_logCN,col='black',lwd=ifelse(sampleList == 'all',0.8,0.5))
        segments(x0=min(xleft),x1 = max(xright), 
                 y0=0.2, y1=0.2,
                 col = 'darkgrey',lty = 'dashed',lwd = 0.4)
        segments(x0=min(xleft),x1 = max(xright), 
                 y0=-0.2, y1=-0.2,
                 col = 'darkgrey',lty = 'dashed',lwd = 0.4)
        
      }
      
    }
    
    
    if(sampleList != 'all'){
      saveFig(file.path(resd,paste0('RCC_CK_',PDID)),plotFun,width = 2.6,height = 2.1,res=500)
    }else{
      saveFig(file.path(res,paste0('FigS1_CK_',PDID,'_v1e7')),plotFun,width = 2.3,height = 3.5,res=500)  
    }
    
  }
}




# Figure 2Ba - Barplots of RCC CK classification ####
rccBarplot.CK = function(){
  rcc.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')
  # Extract CK output
  dd = rcc.srat@meta.data
  dd$finalAnn = dd$finalAnn2
  dd$finalAnn = factor(dd$finalAnn,levels = c('Leukocytes','PTC','Tumour'))
  dd$CKpred.normREF.default.80perc.2397 = ifelse(is.na(dd$CKpred.normREF.default.80perc.2397),'Uncalled',
                                                 ifelse(dd$CKpred.normREF.default.80perc.2397 == 'aneuploid','Aneuploid','Diploid'))
  dd = dd[dd$CKpred.normREF.default.80perc.2397 != 'Uncalled',]
  #dd$CKpred.normREF.default.80perc.2397 = factor(dd$CKpred.normREF.default.80perc.2397,levels = c('Aneuploid','Diploid','Uncalled'))
  dd$CKpred.normREF.default.80perc.2397 = factor(dd$CKpred.normREF.default.80perc.2397,levels = c('Aneuploid','Diploid'))
  
  
  #Define the layout
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    layout(matrix(c(1,2,3,4),ncol=1),heights = c(0.2,1,1,1.8))
    par(mar=c(0,0.6,1,0.6))
    plot(0, 0,
         las=1,
         type='n',frame.plot = F,axes = F)
    title('RCC',cex.main=1,family = 'Helvetica',font=2)
    for(celltype in levels(dd$finalAnn)){
      print(celltype)
      par(mar=c(0.2,1.5,0.5,0.1),xpd=TRUE)
      
      tmp=as.matrix(table(dd$CKpred.normREF.default.80perc.2397[dd$finalAnn == celltype],dd$PDID[dd$finalAnn == celltype]))
      tmp = sweep(tmp,2,colSums(tmp),'/')
      cols = c(Diploid = '#c8c8c8',
               Aneuploid = '#b02c46',
               Uncalled = '#474646')
      if(celltype == levels(dd$finalAnn)[length(levels(dd$finalAnn))]){
        par(mar=c(4.4,1.5,0.5,0.1),xpd=TRUE)
        barplot(tmp,
                col=cols[rownames(tmp)],
                space=0.1,axes = FALSE,
                las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
        #axis(2,at=c(1,0.5,0),las=1,pos = 0.1,tck=-0.02,cex.axis=0.8,lwd.ticks = 0,hadj = 0.5)
        text(x = seq(0.6,4.5,by = 1.105),y = -0.1,colnames(tmp),cex=0.7,family = 'Helvetica',font=1,srt=90,adj = 1)
        #text(x = 4.8,y = 0.5,celltype,cex=0.85,family = 'Helvetica',font=1,srt=270)
        text(x =-1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
      }else{
        barplot(tmp,
                col=cols[rownames(tmp)],
                space=0.1,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                las = 1,main = '',border = F)
        #axis(2,at=c(1,0.5,0),las=1,pos = 0.1,tck=-0.02,cex.axis =0.8,lwd.ticks = 0,hadj = 0.5,)
        #text(x = 4.8,y = 0.5,celltype,cex=0.85,family = 'Helvetica',font=1,srt=270)
        text(x =-1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
      }
      
    }
    
  }
  saveFig(file.path(res,paste0('Fig2Ba_RCC_CK_barplot_noUncalled')),plotFun,width = 0.6,height = 2.75,rawData = dd,res=500)
}





# Figure 2Ea - Barplots of NB CK classification ####
nbBarplot.CK = function(){
  nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')
  # Extract CK output
  dd = nb.srat@meta.data
  dd$finalAnn = factor(dd$finalAnn,levels = c('Leukocytes','Endothelium','Mesenchyme','Tumour'))
  dd$PDID = factor(dd$PD_ID,levels = c("PD42184","PD42752-1","PD42752-2","PD46693","PD43255"))
  dd$CKpred.normREF.default.80perc.2397 = ifelse(is.na(dd$CKpred.normREF.default.80perc.2397),'Uncalled',
                                                 ifelse(dd$CKpred.normREF.default.80perc.2397 == 'aneuploid','Aneuploid','Diploid'))
  dd = dd[dd$CKpred.normREF.default.80perc.2397 != 'Uncalled',]
  #dd$CKpred.normREF.default.80perc.2397 = factor(dd$CKpred.normREF.default.80perc.2397,levels = c('Aneuploid','Diploid','Uncalled'))
  dd$CKpred.normREF.default.80perc.2397 = factor(dd$CKpred.normREF.default.80perc.2397,levels = c('Aneuploid','Diploid'))
  
  #Define the layout
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    layout(matrix(c(1,2,3,4,5),ncol=1),heights = c(0.2,1,1,1,2))
    par(mar=c(0,0.6,0.8,0.6))
    plot(0,0,
         las=1,
         type='n',frame.plot = F,axes = F)
    title('Neuroblastoma',cex.main=1,family = 'Helvetica',font=2)
    for(celltype in levels(dd$finalAnn)){
      print(celltype)
      par(mar=c(0.05,1.5,0.2,0.1),xpd=TRUE)
      #Define the empty plot area
      #plot(0,0,las=1,
      #     type='n',
      #     xlab='',xaxt='n',
      #     ylab='',yaxt='n',frame.plot=T)
      
      tmp=as.matrix(table(dd$CKpred.normREF.default.80perc.2397[dd$finalAnn == celltype],dd$PDID[dd$finalAnn == celltype]))
      tmp = sweep(tmp,2,colSums(tmp),'/')
      cols = c(Diploid = '#c8c8c8',
               Aneuploid = '#b02c46',
               Uncalled = '#474646')
      if(celltype == levels(dd$finalAnn)[length(levels(dd$finalAnn))]){
        par(mar=c(4.9,1.5,0.2,0.1),xpd=TRUE)
        barplot(tmp,
                col=cols[rownames(tmp)],
                space=0.1,axes = FALSE,
                las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
        #axis(2,at=c(1,0.5,0),las=1,pos = 0.1,tck=-0.02,cex.axis=0.75,lwd.ticks = 0,hadj = 0.5)
        #text(x = seq(0.5,6,by = 1.2),y = -0.1,colnames(tmp),cex=0.8,family = 'Helvetica',font=1,srt=45,adj = 1)
        #text(x = 5.8,y = 0.5,celltype,cex=0.65,family = 'Helvetica',font=1,srt=270)
        
        #axis(2,at=c(1,0.5,0),las=1,pos = 0.1,tck=-0.02,cex.axis=0.8,lwd.ticks = 0,hadj = 0.5)
        #text(x = seq(0.5,6,by = 1.2),y = -0.1,colnames(tmp),cex=0.9,family = 'Helvetica',font=1,srt=90,adj = 0.95)
        text(x = seq(0.65,6,by = 1.105),y = -0.1,colnames(tmp),cex=0.7,family = 'Helvetica',font=1,srt=90,adj = 0.95)
        text(x = -1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
      }else{
        barplot(tmp,
                col=cols[rownames(tmp)],
                space=0.1,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                las = 1,main = '',border = F)
        #axis(2,at=c(1,0.5,0),las=1,pos = 0.1,tck=-0.02,cex.axis =0.75,lwd.ticks = 0,hadj = 0.5,)
        text(x =-1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
      }
      
    }
    
  }
  #saveFig(file.path(resDir,paste0('TEST')),plotFun,width = 1.7,height = 3.2,rawData = dd)
  saveFig(file.path(res,paste0('NB_CK_barplot_noUncalled')),plotFun,width = 0.7,height = 3.0,rawData = dd,res=500)
  
}

# Figure 2Eb - Barplots of NB AI classification ####
nbBarplot.AI = function(){
  nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')
  # Extract CK output
  dd = nb.srat@meta.data
  dd$finalAnn = factor(dd$finalAnn,levels = c('Leukocytes','Endothelium','Mesenchyme','Tumour'))
  dd$PDID = factor(dd$PD_ID,levels = c("PD42752-1","PD42752-2","PD46693","PD43255","PD42184"))
  dd$AI_sc_call = ifelse((is.na(dd$AI_sc_call)|dd$AI_sc_call == 'Uncalled'),'Uncalled',
                         ifelse(dd$AI_sc_call == 'abbFrac','Tumour','Normal'))
  dd = dd[dd$AI_sc_call != 'Uncalled',]
  #dd$AI_sc_call = factor(dd$AI_sc_call,levels = c('Tumour','Normal','Uncalled'))
  dd$AI_sc_call = factor(dd$AI_sc_call,levels = c('Tumour','Normal'))
  
  #Define the layout
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    layout(matrix(c(1,2,3,4,5),ncol=1),heights = c(0.2,1,1,1,2.2))
    par(mar=c(0,0.6,0.8,0.6))
    plot(0,0,
         las=1,
         type='n',frame.plot = F,axes = F)
    title('Neuroblastoma',cex.main=1,family = 'Helvetica',font=2)
    for(celltype in levels(dd$finalAnn)){
      par(mar=c(0.05,1.5,0.2,0.1),xpd=TRUE)
      
      tmp=as.matrix(table(dd$AI_sc_call[dd$finalAnn == celltype],dd$PDID[dd$finalAnn == celltype]))
      tmp = sweep(tmp,2,colSums(tmp),'/')
      cols = c(Normal = '#c8c8c8',
               Tumour = '#b02c46',
               Uncalled = '#474646')
      if(celltype == levels(dd$finalAnn)[length(levels(dd$finalAnn))]){
        par(mar=c(4.9,1.5,0.2,0.1),xpd=TRUE)
        barplot(tmp,
                col=cols[rownames(tmp)],
                space=0.1,axes = FALSE,
                las = 1,names.arg = rep(NA,ncol(tmp)),border = F)  
        
        text(x = seq(0.65,6,by = 1.105),y = -0.1,colnames(tmp),cex=0.7,family = 'Helvetica',font=1,srt=90,adj = 0.95)
        text(x = -1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
      }else{
        barplot(tmp,
                col=cols[rownames(tmp)],
                space=0.1,axes = FALSE,names.arg = rep(' ',ncol(tmp)),
                las = 1,main = '',border = F)
        text(x =-1.3,y = 0.5,celltype,cex=0.7,family = 'Helvetica',font=1,srt=90)
      }
      
    }
    
  }
  saveFig(file.path(res,paste0('Fig2Eb_NB_AI_barplot_noUncalled')),plotFun,width = 0.7,height = 3.0,rawData = dd,res=500)
  
}


# Figure 2F - Example of CopyKat output for normal and tumour NB cells ####
sampleList = c('PD42184','PD42752-2')
# Figure S1 NB - CopyKat output for normal and tumour NB cells ####
subCl.minSegLen = 2e7
nbCK = function(){
  nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')
  nb.srat$PDID = as.character(nb.srat$PD_ID)
  # Import chromInfo
  chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')
  # Import Manifest
  projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = "NB_mani")
  # CopyKat results
  copykat.results = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_nb.CK.normREF.default.2397/v5_CKresults_default_normREF_80perc_2397.rds')
  
  #----- Processing CopyKat results -------#
  for(i in 1:length(copykat.results)){
    # Get copyKat CNA matrix and prediction
    CNA_summary_byCellType = data.frame()
    CNA_mat = copykat.results[[i]]$CNAmat
    colnames(CNA_mat) = gsub('^X','',colnames(CNA_mat))
    colnames(CNA_mat) = gsub('\\.','-',colnames(CNA_mat))
    pred = as.data.frame(copykat.results[[i]]$prediction)
    
    sample = unique(nb.srat@meta.data[rownames(nb.srat@meta.data) %in% rownames(pred),]$PDID)
    PDID=sample
    
    if(length(sample) > 1){
      message(paste0('More than 1 sample detected: i=',i,', samples are ',sample))
    }else if(length(sample) <1){
      next
    }else if(length(sample)==1){
      message(paste0('Checking sample ',sample))
    }
    
    # subset annnb.srat object to keep only cells of that sample
    srat = subset(nb.srat, subset = PDID == sample)
    
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
    
    # Remove X chromosome 
    CNA_summary_byCellType = CNA_summary_byCellType[CNA_summary_byCellType$chrom != 23,]
    CNA_summary_byCellType$celltype = factor(CNA_summary_byCellType$celltype,levels=c("Tumour","Mesenchyme", "Endothelium",'Leukocytes'))
    
    
    ####------------------ Generate Battenberg CN summary file ----------------####
    donorMani = projMani[projMani$PDID == sample,]
    btb.fp = unique(donorMani$battenbergFp)
    #----- Processing Battenberg data -------#
    dna.data = annotateBTB(btb.fp,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=c(1:22),removeBalancedSegs=F,longFormat = T,method = 'totalCN')  
    # Remove X chromosome 
    dna.data = dna.data[dna.data$Chr != 23,]
    
    dna.data$celltype = 'Tumour'
    tmp = rbind(data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype='Endothelium'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype='Endothelium'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype='Leukocytes'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype='Leukocytes'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Start',pos=1,idx=0,abspos_kb=0,celltype='Mesenchyme'),
                data.frame(Idx=0,Chr=1,matNum=1,patNum=1,frac=1,segLen=0,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=0,posID=0,posType='Stop',pos=2874771,idx=0,abspos_kb=2874771,celltype='Mesenchyme'))
    
    
    dna.data = rbind(dna.data,tmp)
    dna.data$log_CNratio = log(as.numeric(dna.data$tumTot)/2)
    subCl.dna = dna.data[dna.data$type == 'sub',]
    majCl.dna = dna.data[dna.data$type == 'maj',]
    
    
    
    #----- Plotting copykat results! -------#
    chromInfo2 = chromInfo[chromInfo$chrom != 23,]
    
    plotFun = function(noFrame=TRUE,noPlot=FALSE){
      noFrame=T
      # Set layout
      layout(mat=matrix(c(1:nlevels(CNA_summary_byCellType$celltype)),ncol=1),
             heights = rep(2,nlevels(CNA_summary_byCellType$celltype)))
      
      for(celltype in levels(CNA_summary_byCellType$celltype)){
        par(mar=c(0.2,0.6,0.8,0.6),xpd=TRUE)
        tmp = CNA_summary_byCellType[CNA_summary_byCellType$celltype == celltype,]
        dna = majCl.dna[majCl.dna$celltype == celltype,]
        ncells = nrow(nb.srat@meta.data[nb.srat@meta.data$PDID == PDID & nb.srat@meta.data$finalAnn == celltype,])
        
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
        
        text(x=9e8,y=ytext,paste0(celltype,'_',PDID,' (n=',ncells,')'),cex=0.6,family = 'Helvetica',font=2,adj = 0)
        
        axis(2,at=c(ybottom+0.05,0,ytop-0.05),labels = c('Low',0,'High'),las=1,pos = 0,tck = -.00,lwd = 0.3,cex.axis=0.6,hadj = 0.3,padj = 0.5)
        axis(4,las=1,pos = max(chromInfo2$abspos_kb*1000),tck = -.02,lwd = 0.7,cex.axis=0.6,hadj = 1.5,col='black',
             at=c(round(log(0.5)/2,2),0,round(log(3/2)/2,2),round(log(4/2)/2,2)),col.axis = '#b02c46',
             labels = c(1,2,3,4))
        
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
             ytop=dna[dna$posType=='Start',]$log_CNratio/2,
             col=colAlpha('#b02c46',1),
             border=colAlpha('#b02c46',1),lwd = 0.8)
        
        # Plot subclone total CN
        if(celltype == 'Tumour' & nrow(subCl.dna) > 0){
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
    
    saveFig(file.path(res,paste0('FigS2_NB_CK_',PDID,'_',subCl.minSegLen)),plotFun,width = 2.0,height = 3.2,res=500)  
    
  }
}


# Figure 2C and S2 - Example of AlleleIntegrator output for normal and tumour RCC cells ####
subCl.minSegLen = 2e7
cov=500
tumour = 'NB'
sampleList = 'all'

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





# GOSH25 - bulkDNA BAF plot ####



# Plot Tumour bulkDNA BAF #### 
setwd('~/lustre_mt22/CN_methods/')
#############
# Libraries #
#############
library(alleleIntegrator)

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
minSegLen=1e6
subCl.minSegLen=2e7
skipIfExists = T
# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204'


refGenome = '/lustre/scratch119/realdata/mdt1/team78pipelines/reference/Human/GRCH37d5/genome.fa'
nParallel=25

#----------- Run AlleleIntegrator on NB dataset (5 samples)
#----------------------------------------------------------- 
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('Ewings','Wilms','ATRT')){
    for(PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Generating BAF plot from bulkDNA for Sample %s - tumourType: %s',PDID,tumourType))
      # Set output directory
      outDir = file.path(mainDir,'alleleIntegrator_output',tumourType,PDID)
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
      
      baf.out = generateCoverageAndBAF(BAM = tumourDNA,refGenome = refGenome,hSNPs=hSNPs,
                                       outPath = file.path(outDir,paste0(PDID,'_cov_BAF.RDS')),nParallel=nParallel)
      
      plotFun = function(noFrame=T,noPlot=FALSE,minCoverage=10){
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
        
        #axis(side=1, at=cumsum(c(0,chrLens[-length(chrLens)]))+chrLens/2, labels = chrs)
        axis(side=2, at=c(0,0.5,1),labels=c(0,0.5,1),las=1)
        abline(v=cumsum(chrLens),col='lightgrey')
        
      }
      
      saveFig(file.path(res,paste0('BAF_',PDID,'_',tumourType)),plotFun,width = 5.8,height = 2.2,res=1000)
    } 
  }
}










baf.out = generateCoverageAndBAF(BAM = tumourDNA,refGenome = refGenome,hSNPs=hSNPs,outPath = paste0('/lustre/scratch117/casm/team274/mt22/CN_methods/',PD46693,'_cov_BAF.RDS',nParallel=24)
                                 
                                 plotFun = function(noFrame=T,noPlot=FALSE,minCoverage=10){
                                   #Work out the chromosome boundaries
                                   chrsToPlot=c(1:22)
                                   chrs = chrsToPlot
                                   chrLens = seqlengths(filt)
                                   tmp = sapply(split(start(filt),as.character(seqnames(filt))),max)
                                   chrLens[is.na(chrLens)] = tmp[names(chrLens)[is.na(chrLens)]]
                                   chrLens = as.numeric(chrLens[chrs])
                                   
                                   #Filter to just the ones that we trust
                                   filt = hSNPs[hSNPs$coverage>=minCoverage,]
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
                                   
                                   #axis(side=1, at=cumsum(c(0,chrLens[-length(chrLens)]))+chrLens/2, labels = chrs)
                                   axis(side=2, at=c(0,0.5,1),labels=c(0,0.5,1),las=1)
                                   abline(v=cumsum(chrLens),col='lightgrey')
                                   
                                 }
                                 
                                 saveFig(file.path(resd,paste0('BAF_',PDID)),plotFun,width = 5.8,height = 2.2,res=1000)
                                 
                                 