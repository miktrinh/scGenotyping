# Run Copy Kat
setwd('/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/')

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
source('../scripts/finalScripts/misc.R')


runCopyKat = function(seuratObject,seed=2397,ident,n.cores=20, normREF=F,normREF.ident,norm.cells.ref = "",normREF_p = 1,title=''){
  if(normREF){
    title = paste0(title,'.normREF.default.',seed)
    col_to_add = paste0('CKpred.normREF.default.',normREF_p*100,'perc.',seed) 
    seuratObject[[paste0('CK.normREF.',seed,'.',normREF_p*100,'perc')]] = NA
  }else{
    title = paste0(title,'.default.',seed)
    col_to_add = paste0('CKpred.default.',seed)  
  }
  
  if(!dir.exists(file.path('./',title))){
    dir.create(file.path('./',title))  
  }else{
    message('Existing output directory found... Please check!')
    return()
  }
  
  setwd(file.path('./',title))
  
  seuratObject[[col_to_add]] = NA
  
  copykat.results = list()
  
  for (i in 1:nrow(unique(seuratObject[[ident]]))){
    sample = unique(seuratObject[[ident]])[i,1]
    print(paste0('Analyzing sample ',sample))
    srat = subset(seuratObject,subset = PD_ID == sample)
    mtx = srat@assays$RNA@counts  
    
    if(normREF){
      norm.cells = rownames(srat@meta.data[srat@meta.data$cell_type == normREF.ident,])  
      set.seed(seed)
      norm.cells.ref = sample(norm.cells,size = length(norm.cells)*normREF_p)
      m = match(norm.cells.ref,rownames(seuratObject@meta.data))
      sum(is.na(m))
      seuratObject[[paste0('CK.normREF.',seed,'.',normREF_p*100,'perc')]][m,1] = TRUE  
    }
    
    set.seed(seed)
    copykat.results[[i]] = copykat_inhouse(rawmat=mtx, id.type="S", sam.name=sample, n.cores=n.cores, norm.cell.names = norm.cells.ref, seed = seed)
    names(copykat.results)[i] = sample
    
    pred = as.data.frame(copykat.results[[i]]$prediction)
    m = match(rownames(pred),rownames(seuratObject@meta.data))
    sum(is.na(m))
    
    seuratObject[[col_to_add]][m,1] = pred[,'copykat.pred']  
  }
  
  return(list(seuratObject,copykat.results))
}


###------- Run CopyKat on Neuroblastoma data ----------###
#--------------------------------------------------------- 
nb.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann.rds')
#nb.srat@meta.data = nb.srat@meta.data[,-c(33:34)]

# Import Manifest
projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = "NB_mani")
View(projMani)

# 2. Run CopyKat on individual samples - WITHOUT NORMAL REF
# output is a list of:
# 1. seuratObject 
# 2. copykat.results
REF_output_list = runCopyKat(seuratObject = nb.srat,ident = 'PD_ID',title = 'nb.CKnoREF')
srat=REF_output_list[[1]]
copykat.results = REF_output_list[[2]]



# 3. Run CopyKat on individual samples - WITH NORMAL REF
REF_output_list = runCopyKat(seuratObject = nb.srat,ident = 'PD_ID',normREF = T,
                             normREF.ident='Leukocytes',title = 'v6_nb.CK')
srat=REF_output_list[[1]]
norm.copykat.results = REF_output_list[[2]]


# 4. Save results ####
m = match(rownames(nb.srat@meta.data),rownames(srat@meta.data))
sum(is.na(m))

nb.srat@meta.data$CK.normREF.2397.100perc = srat@meta.data$CK.normREF.2397.100perc[m]
nb.srat@meta.data$CKpred.normREF.default.100perc.2397 = srat@meta.data$CKpred.normREF.default.100perc.2397[m]


saveRDS(nb.srat,'/lustre/scratch117/casm/team274/mt22/CN_methods/NB_ann_CK.100pctNorm.rds')
saveRDS(copykat.results, '/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/nb.CKnoREF.default.2397/CKresults_default_2397.rds')
saveRDS(norm.copykat.results, '/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v6_nb.CK.normREF.default.2397/v6_CKresults_default_normREF_100perc_2397.rds')


###------------ Run CopyKat on RCC data ---------------###
#--------------------------------------------------------- 
rcc.srat = readRDS('/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')
rcc.srat = subset(rcc.srat,subset = finalAnn %in% c('Renal_cell_carcinoma','Proximal_tubuluar_cells','Leukocytes'))
message(sprintf('Total number of input cells: %d',dim(rcc.srat)[2]))

# Import Manifest
projMani = read_excel("/lustre/scratch117/casm/team274/mt22/projectManifest.xlsx",sheet = "RCC_mani")
View(projMani)


# 2. Run CopyKat on individual samples - WITHOUT NORMAL REF
# output is a list of:
# 1. seuratObject 
# 2. copykat.results
kid.noREF_output_list = runCopyKat(seuratObject = rcc.srat,ident = 'PDID',title = 'rcc.CKnoREF')
rcc.srat=kid.noREF_output_list[[1]]
kid.copykat.results = kid.noREF_output_list[[2]]




# 3. Run CopyKat on individual samples - WITH NORMAL REF
kid.REF_output_list = runCopyKat(seuratObject = rcc.srat,ident = 'PDID',normREF = T,
                                 normREF.ident='Leukocytes',title = 'v5_rcc.CK')
rcc.srat=kid.REF_output_list[[1]]
kid.norm.copykat.results = kid.REF_output_list[[2]]


# 4. Save results ####
m = match(rownames(rcc.srat@meta.data),rownames(rcc.srat@meta.data))
sum(is.na(m))
rcc.srat@meta.data$CK.normREF.2397.80perc = rcc.srat@meta.data$CK.normREF.2397.80perc[m]
rcc.srat@meta.data$CKpred.normREF.default.80perc.2397 = rcc.srat@meta.data$CKpred.normREF.default.80perc.2397[m]

saveRDS(rcc.srat,'/lustre/scratch117/casm/team274/mt22/CN_methods/RCC_PCT_ann3sub.rds')
saveRDS(kid.copykat.results, '/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/rcc.CKnoREF.default.2397/CKresults_default_2397.rds')
saveRDS(kid.norm.copykat.results, '/lustre/scratch117/casm/team274/mt22/CN_methods/CopyKAT_output/v5_rcc.CK.normREF.default.2397/v5_CKresults_default_normREF_80perc_2397.rds')


