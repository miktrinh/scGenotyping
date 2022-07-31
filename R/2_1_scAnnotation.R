# Process scRNAseq data for tumour samples
setwd('~/lustre_mt22/CN_methods/')
#############
# Libraries #
#############
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(readxl)
source('scripts/finalScripts/R/qc10X.R')


# Import Manifest
projMani = read_excel("../projectManifest.xlsx",sheet = "alleleIntegrator")
mainDir = '~/lustre_mt22/CN_methods/revision_2204'


#------------------------------------------------------# 
# ** Perform QC on single cell 10X raw output channels
#------------------------------------------------------#

#############
#   Params  #
#############
maxMT = 20
minGenes = 200
minUMIs = 600
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = F
scrubScoreMax = 0.5
scrubPath='../cleanCounts/scrubletScores.tsv'
scPath="../cleanCounts/strainedCounts"
#scrubPath=NULL
#scPath=NULL
doPlot=T
verbose = T
skipIfExists=F
keepMTCells=T


qc.summary = tibble()
for(tumourType in unique(projMani$TumourType)){
  if(tumourType %in% c('Wilms','ATRT','Ewings')){
    for(PDID in unique(projMani$PDID[projMani$TumourType == tumourType])){
      message(sprintf('Sample %s - %s: Processing raw 10X channels',PDID,tumourType))
      outDir = file.path(mainDir,'sc_seuratObjects',tumourType,PDID)
      
      if(!file.exists(outDir)){
        print('Making new Dir')
        dir.create(outDir,recursive = T)
      }
      donorMani = projMani[projMani$PDID == PDID,]
      
      # Get 10X cellranger dataDirs
      dataDirs = donorMani$cellrangerDir
      names(dataDirs) = donorMani$sangerSampleID
      names(dataDirs) = gsub('_','.',names(dataDirs))
      dataDirs=dataDirs[file.exists(dataDirs)]

      # Collect relevant metadata
      metadata = donorMani[,which(colnames(donorMani) %in% c("PDID",'SampleID','Sex',"TumourType",'sangerSampleID'))]
      metadata$orig.ident = gsub('_','.',metadata$sangerSampleID)
      # Run basicQC
      plotDir = file.path(outDir,paste0(PDID,'_'))
      outPath = file.path(outDir,paste0(PDID,'_'))
      
      QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,clusteringRes=clusteringRes,
                          skipScrub=skipScrub,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                          metadata=metadata,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                          doPlot=doPlot,plotDir=plotDir,verbose=verbose)
      
      df.out = QC.output[[2]]
      qc.summary=rbind(qc.summary,df.out)

    }
  }
}


write.csv(qc.summary,file.path(mainDir,'sc_seuratObjects/qc_summary.csv'))


#------------------------------------------------------# 
#             ** Clustering and Annotation
#------------------------------------------------------#
#############
# Libraries #
#############
library(SoupX)
source("/lustre/scratch117/casm/team274/mt22/generalScripts/utils/logisticRegression.R")
source("/lustre/scratch117/casm/team274/mt22/generalScripts/utils/runLR.R")

#############
#   Params  #
#############
keepMTCells=F
numPCs = 75
clusteringRes=1
seed = 2397

processSrat = function(srat,numPCs=75,clusteringRes=1,min.dist=0.3){
  srat = NormalizeData(srat,verbose=FALSE)
  srat = FindVariableFeatures(srat,verbose=FALSE)
  srat = ScaleData(srat,features = rownames(srat),verbose=FALSE)
  srat = RunPCA(srat,npcs=numPCs,approx=FALSE,verbose=FALSE)
  srat = FindNeighbors(srat,dims=seq(numPCs),verbose=FALSE)
  srat = FindClusters(srat,res=clusteringRes,verbose=FALSE)
  srat = RunUMAP(srat,dims = seq(numPCs),min.dist = min.dist, verbose=FALSE)
  
  return(srat)
}


#-----------#
#  Ewings   #####
#-----------#
if(keepMTCells){
  e1 = readRDS(file.path(mainDir,'sc_seuratObjects/Ewings/PD42181/PD42181_clean_withMTCells.RDS'))
  e2 = readRDS(file.path(mainDir,'sc_seuratObjects/Ewings/PD47706/PD47706_clean_withMTCells.RDS'))
}else{
  e1 = readRDS(file.path(mainDir,'sc_seuratObjects/Ewings/PD42181/PD42181_clean_noMTCells.RDS'))
  e2 = readRDS(file.path(mainDir,'sc_seuratObjects/Ewings/PD47706/PD47706_clean_noMTCells.RDS'))
}
e = merge(e1,e2)

# Standard Seurat scRNA-seq processing and clustering
e = processSrat(e,numPCs = 30,clusteringRes=1)

## Get marker for each cluster - soupX
ewings.markers = quickMarkers(toc = e@assays$RNA@counts,clusters = e$seurat_clusters,N = 30)
DoHeatmap(e,features = ewings.markers$gene,slot = 'data')

pdf(file.path(mainDir,'sc_seuratObjects/Ewings/Ewings_UMAP.pdf'))
DimPlot(e, group.by = 'orig.ident',label = T)
DimPlot(e, group.by = 'Phase',label = T)
DimPlot(e, group.by = 'seurat_clusters',label = T)


e$annot = as.character(e$seurat_clusters)
e$annot[!e$annot %in% c('7','6','9','10')] = 'Tumour'
e$annot[e$annot %in% c('7','6','9')] = 'Leukocytes'
e$annot[e$annot %in% c('10')] = 'Endothelium'
DimPlot(e, group.by = 'annot',label = T)

# Cluster 11 shows SCX (Tendon markers) + NNAT, NEUROD6(neural markers)

FeaturePlot(e,features = c('nFeature_RNA','nCount_RNA','percent.mt'))
FeaturePlot(e,features = c('PTPRC','TRAC','MS4A7','MS4A1','HLA-DRA','LYZ'))
FeaturePlot(e,features = c("PLVAP", "PECAM1", "KDR","PTPRB",'IGFBP7','COL4A1'))
FeaturePlot(e,features = c('CD99','EPCAM','ATP1A1', 'BCL11B', 'GLG1', 'FLI1','EWSR1','TFRC','NEUROD6'))

DotPlot(e,features = c('PTPRC','TRAC','MS4A7','MS4A1','CD99','PLVAP','ATP1A1', 'BCL11B', 'GLG1', 'FLI1','EWSR1','TFRC','NEUROD6'))+
  RotatedAxis()

dev.off()


saveRDS(e,file.path(mainDir,'sc_seuratObjects/Ewings/Ewings_ann.RDS'))





#-----------#
#   Wilms   ####
#-----------#
if(keepMTCells){
  wilms = readRDS(file.path(mainDir,'sc_seuratObjects/Wilms/PD48777/PD48777_clean_withMTCells.RDS'))
}else{
  wilms = readRDS(file.path(mainDir,'sc_seuratObjects/Wilms/PD48777/PD48777_clean_noMTCells.RDS'))
}

# Standard Seurat scRNA-seq processing and clustering
wilms = processSrat(wilms,numPCs = 40,clusteringRes=2,min.dist = 0.4)

## Get marker for each cluster - soupX
wilms_markers = quickMarkers(toc = wilms@assays$RNA@counts,clusters = wilms$seurat_clusters,N = 30)
DoHeatmap(wilms,features = wilms_markers$gene[wilms_markers$cluster %in% c(0:11)],cells = WhichCells(wilms,idents = c(0:11)),slot = 'data')
DoHeatmap(wilms,features = wilms_markers$gene[wilms_markers$cluster %in% c(12:15)],cells = WhichCells(wilms,idents = c(12:15)),slot = 'data')

pdf(file.path(mainDir,'sc_seuratObjects/Wilms/Wilms_UMAP.pdf'))
DimPlot(wilms, group.by = 'Phase',label = T)
DimPlot(wilms, group.by = 'orig.ident',label = T)
DimPlot(wilms, group.by = 'seurat_clusters',label = T)

wilms$annot2 = as.character(wilms$seurat_clusters)
wilms$annot2[!wilms$annot2 %in% c('19','21','4','7','9','2','23','22','14')] = 'Tumour'
wilms$annot2[wilms$annot2 %in% c('19','21')] = 'Leukocytes'
wilms$annot2[wilms$annot2 %in% c('4','7','9','2')] = 'Epithelium'
wilms$annot2[wilms$annot2 %in% c('23')] = 'Endothelium'
wilms$annot2[wilms$annot2 %in% c('22')] = 'Fetal Muscle'
wilms$annot2[wilms$annot2 %in% c('14')] = 'Smooth Muscle'

wilms$annot = as.character(wilms$seurat_clusters)
wilms$annot[!wilms$annot %in% c('19','21','23','22','14')] = 'Tumour'
wilms$annot[wilms$annot %in% c('19','21')] = 'Leukocytes'
#wilms$annot[wilms$annot %in% c('4','7','9','2')] = 'Epithelium'
wilms$annot[wilms$annot %in% c('23')] = 'Endothelium'
wilms$annot[wilms$annot %in% c('22','14')] = 'Muscle'
#wilms$annot[wilms$annot %in% c('14')] = 'Smooth Muscle'

DimPlot(wilms, group.by = 'annot',label = T)
DimPlot(wilms, group.by = 'annot2',label = T)
FeaturePlot(wilms,features = c('EPCAM','PAX2','ITGA8','HMGA1')) #+ ggtitle('Epithelium Markers')
FeaturePlot(wilms,features = c('PTPRC','MS4A7','LYZ','HLA-DRA','CD3D','GZMA')) #+ ggtitle('Leukocytes Markers')
FeaturePlot(wilms,features = c("PLVAP", "PECAM1", "KDR","PTPRB")) #+ ggtitle('Endothelium Markers')
FeaturePlot(wilms,features = c("SIX1", "SIX2", "WT1")) # Wilms
FeaturePlot(wilms,features = c("ACTC1", "MYL4", "MYOG",'TTN')) # Fetal muscle
FeaturePlot(wilms,features = c("COL8A1",'ACTA2','CRYAB','COX7A1')) # Smooth Muscle
FeaturePlot(wilms,features = c('nFeature_RNA','nCount_RNA','percent.mt'))

dev.off()

saveRDS(wilms,file.path(mainDir,'sc_seuratObjects/Wilms/Wilms_ann.RDS'))






#-----------#
#   ATRT   ####
#-----------#
if(keepMTCells){
  atrt = readRDS(file.path(mainDir,'sc_seuratObjects/ATRT/PD47705/PD47705_clean_withMTCells.RDS'))
}else{
  atrt = readRDS(file.path(mainDir,'sc_seuratObjects/ATRT/PD47705/PD47705_clean_noMTCells.RDS'))
}

# Standard Seurat scRNA-seq processing and clustering
atrt = processSrat(atrt,numPCs = 50,clusteringRes=1,min.dist = 0.4)

## Get marker for each cluster - soupX
atrt_markers = quickMarkers(toc = atrt@assays$RNA@counts,clusters = atrt$seurat_clusters,N = 30)
DoHeatmap(atrt,features = atrt_markers$gene,slot = 'data')
DoHeatmap(wilms,features = wilms_markers$gene[wilms_markers$cluster %in% c(12:15)],cells = WhichCells(wilms,idents = c(12:15)),slot = 'data')

pdf(file.path(mainDir,'sc_seuratObjects/ATRT/markers.pdf'),height = 40,width = 10)
DoHeatmap(atrt,features = m$gene[m$cluster %in% c('12','15','14','16','17')],slot = 'count',cells = WhichCells(atrt,idents = c('12','15','14','16','17')))
dev.off()

pdf(file.path(mainDir,'sc_seuratObjects/ATRT/ATRT_UMAP.pdf'))
DimPlot(atrt, group.by = 'Phase',label = T)
DimPlot(atrt, group.by = 'orig.ident',label = T)
DimPlot(atrt, group.by = 'seurat_clusters',label = T)

atrt$annot = as.character(atrt$seurat_clusters)
atrt$annot[!atrt$annot %in% c('12','15','14','16','17')] = 'Tumour'
atrt$annot[atrt$annot %in% c('17')] = 'Endothelium'
atrt$annot[atrt$annot %in% c('16')] = 'Fibroblast'
atrt$annot[atrt$annot %in% c('14')] = 'Leukocytes'
atrt$annot[atrt$annot %in% c('12','15')] = 'Epithelium'
DimPlot(atrt, group.by = 'annot',label = T)

FeaturePlot(atrt,features = c("SMARCB1"))
FeaturePlot(atrt,features = c("PLVAP", "PECAM1", "KDR","PTPRB")) # Endothelium
FeaturePlot(atrt,features = c('CD14','CD68','PTPRC','HLA-B')) # Leukocytes
FeaturePlot(atrt,features=c('HIGD1B','COX4I2','NOTCH3','ENPEP')) # Fibroblasts
FeaturePlot(atrt,features = c("EPCAM",'KRT18','HOXB9','PAX2')) # Epithelium
FeaturePlot(atrt,features = c('nFeature_RNA','nCount_RNA','percent.mt'))

dev.off()

saveRDS(atrt,file.path(mainDir,'sc_seuratObjects/ATRT/ATRT_ann.RDS'))
































#### 1. Import fetal REF data ####
kidREF = readRDS('/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/preProcess/scData/finalisedSeurat_fKid.RDS')
adrREF = readRDS('/lustre/scratch117/casm/team274/mt22/abnormal_karyotypes/adrenal/adrREF.rds')
adrREF$annot = adrREF$cell_type
REF.srat = merge(kidREF,adrREF)

#### 4. Logistic Regression ####
ref_annot='annot'
LR_level='both'
srat_annot=''
tissue = ''
minGeneMatch = 0.99
maxCells=4000

scLR_TGTtype = ''
scLR_REFtype = 'cell.labels'

plotPath = '~/lustre_mt22/CN_methods/sc_seuratObjects/ewings_adrKid_'
outPath = '~/lustre_mt22/CN_methods/sc_seuratObjects/ewings_adrKid_'

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=e,LR_level='both',srat_annot='',minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,scLR_TGTtype='',scLR_REFtype='')
e2 = outputs[[1]]
DimPlot(e2, group.by = 'LR_celltype.scLevel',label = T,repel = T)
DimPlot(e, group.by = 'seurat_clusters',label = T)
DimPlot(e, group.by = 'orig.ident',label = T)
LR_outputs = outputs[[2]]
