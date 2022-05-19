#############
# Libraries #
#############
# Load libraries
#library(tidyverse)
library(Seurat)
library(SoupX)
library(RColorBrewer)


#' Run scrublet
#' 
#' Uses a crude call to python to run scrublet and pull in the results.
#'
#' @param dat10X Directory containing 10X matrix and genes files.
#' @param nPCs Number of PCs to use.
#' @return A data.frame with the scrublet results
runScrublet = function(dat10X,nPCs){
  fNom = tempfile()
  cmd = sprintf(
    "import scrublet as scr
import scipy.io
import numpy as np
import os
import gzip
os.chdir('%s')
input_dir = os.path.expanduser('%s')
output_file = '%s'
try:
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
  genes = np.array(scr.load_genes(input_dir + '/genes.tsv', delimiter='\\t', column=1))
except FileNotFoundError:
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
  #Copy relevant part of scr.load_genes
  gene_list = []
  gene_dict = {}
  with gzip.open(input_dir + '/features.tsv.gz','rt') as f:
    for l in f:
      gene = l.strip('\\n').split('\\t')[1]
      if gene in gene_dict:
        gene_dict[gene] += 1
        gene_list.append(gene + '__' + str(gene_dict[gene]))
        if gene_dict[gene] == 2:
          i = gene_list.index(gene)
          gene_list[i] = gene + '__1'
      else: 
        gene_dict[gene] = 1
        gene_list.append(gene)
  genes = np.array(gene_list)
#Convert to integers.  Will do nothing most of the time.
counts_matrix.data = np.round(counts_matrix.data).astype('int32')
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                        min_cells=3, 
                                                        min_gene_variability_pctl=85, 
                                                        n_prin_comps=%d)
with open(output_file,'w') as f:
  #Header row
  f.write('score\\tcall\\tsimScore1\\tsimScore2\\n')
  #As nSim = 2*nObs, we can split into two columns 
  offset = len(doublet_scores)
  f.write('\\n'.join([str(x)+'\\t'+str(predicted_doublets[i]).upper()+'\\t'+str(scrub.doublet_scores_sim_[i])+'\\t'+str(scrub.doublet_scores_sim_[offset+i]) for i,x in enumerate(doublet_scores)]))
  f.close()",getwd(),dat10X,fNom,nPCs)
  system(paste0('python -c "',cmd,'"'))
  out = read.table(fNom,header=TRUE,sep='\t')
  if(file.exists(file.path(dat10X,'barcodes.tsv'))){
    tmp = read.table(file.path(dat10X,'barcodes.tsv'),header=FALSE)
  }else{
    tmp = read.table(file.path(dat10X,'barcodes.tsv.gz'),header=FALSE)
  }
  out$barcode = tmp[,1]
  return(out)
}



# To overwrite existing output, just specify scPath = NULL
cleanWithSoupX = function(cleanSrat,dataDirs,scPath,numPCs=75,clusteringRes=10,numPCs_soupX=75,clusteringRes_soupX=1,plotDir=NULL){
  #### =================================================== #
  #### Run SoupX to remove ambiant reads/mRNA
  # Call doublets with Scrublet
  if(is.null(dataDirs)){
    warning("Data directories not named, generating automatic names.")
    names(dataDirs) = paste0('DataSource',seq_along(dataDirs))
  }
  
  
  if(verbose)
    message("Loading or running soupX")
  if(!is.null(scPath)){
    soupRun = !file.exists(file.path(dataDirs,scPath))
  }else{
    soupRun = rep(TRUE,length(dataDirs))
  }
  
  
  bigSrat = list()
  toRemove = c()
  for(i in seq_along(dataDirs)){
    print(i)
    
    if(soupRun[i]){
      if(verbose)
        message(sprintf('Running soupX for channel %s',names(dataDirs)[i]))
      
      
      setwd(file.path(dataDirs[i],'../'))
      if(!file.exists(file.path(dataDirs[i],'../cleanCounts'))){
        if(verbose)
          message(sprintf('Creating Clean Count dir for channel %s',names(dataDirs)[i]))
        dir.create(file.path(dataDirs[i],'../cleanCounts'))
      }
      # Clustering
      if(sum(cleanSrat$orig.ident == names(dataDirs)[i]) <= 10){
        message(sprintf('PDID %s - channel %s : %d cells identified for SoupX - Skipping',PDID,names(dataDirs)[i],sum(cleanSrat$orig.ident == names(dataDirs)[i])))
        toRemove = c(toRemove,i)
        next()
      }
      srat = subset(cleanSrat,subset = orig.ident == names(dataDirs)[i])
      srat = NormalizeData(srat,verbose=FALSE)
      srat = FindVariableFeatures(srat,verbose=FALSE)
      srat = ScaleData(srat,verbose=FALSE)
      srat = RunPCA(srat,npcs=numPCs_soupX,approx=FALSE,verbose=FALSE)
      srat = FindNeighbors(srat,dims=seq(numPCs_soupX),verbose=FALSE)
      srat = FindClusters(srat,res=clusteringRes_soupX,verbose=FALSE)
      srat = RunUMAP(srat,dims=seq(numPCs_soupX),verbose=FALSE)
      
      srat@meta.data$RD1 = srat@reductions$umap@cell.embeddings[,1]
      srat@meta.data$RD2 = srat@reductions$umap@cell.embeddings[,2]
      
      #DimPlot(srat)
      #if(names(dataDirs)[i] %in% c('adr66')){
      #  srat.souped = srat
      #  bigSrat[[i]] = srat.souped
      #  names(bigSrat)[i] = as.character(unique(srat.souped@meta.data$orig.ident))  
      #}
      datadir = file.path(dataDirs[i],'../')
      names(datadir) = names(dataDirs)[i]
      srat.souped = runSoupX(datadir = datadir,toc = srat@assays$RNA@counts,metadata = srat@meta.data,plotDir=plotDir,scPath=scPath)
      #srat.souped = srat.out
      
      
    }else{
      # If results do exist and do not overwrite existing output
      message(sprintf('Existing data found for %s - Loading...',names(dataDirs[i])))
      srat.souped = Read10X(file.path(dataDirs[i],scPath))
      srat.souped = CreateSeuratObject(srat.souped)
      
    }
    
    # Add clean count to bigSrat
    
    bigSrat[[i]] = srat.souped
    names(bigSrat)[i] = as.character(unique(srat.souped@meta.data$orig.ident))  
    
  }
  
  if(!is.null(toRemove)){
    bigSrat = bigSrat[-toRemove]  
  }
  
  # Merge all clean channel srat objects into a big object
  if(length(bigSrat) >1){
    bigSrat2 = merge(bigSrat[[1]],bigSrat[-1])  
  }else{
    bigSrat2 = bigSrat[[1]]
  }
  
  # Check if all cells in cleanSeurat are included
  if(sum(!rownames(cleanSrat@meta.data) %in% rownames(bigSrat2@meta.data)) >0){
    message(sprintf('%d cells in the provided cleanSrat are not included in soupedSrat output',sum(!rownames(cleanSrat@meta.data) %in% rownames(bigSrat2@meta.data))))  
  }
  
  # Clustering
  #bigSrat2 = NormalizeData(bigSrat2,verbose=FALSE)
  #bigSrat2 = FindVariableFeatures(bigSrat2,verbose=FALSE)
  #bigSrat2 = ScaleData(bigSrat2,verbose=FALSE)
  #bigSrat2 = RunPCA(bigSrat2,npcs=numPCs,approx=FALSE,verbose=FALSE)
  #bigSrat2 = FindNeighbors(bigSrat2,dims=seq(numPCs),verbose=FALSE)
  #bigSrat2 = FindClusters(bigSrat2,res=clusteringRes,verbose=FALSE)
  #bigSrat2 = RunUMAP(bigSrat2,dims=seq(numPCs),verbose=FALSE)
  
  return(bigSrat2)
  
}


#datadir = file.path(dataDirs[i], "../")
#toc = srat@assays$RNA@counts
#metadata = srat@meta.data

#' Run soupX
#'
#' @param datadir Directory containing 10X cellranger output.
#' @param toc Table of UMI Count 
#' @param metadata Metadata table for the list of cells under investigation 
#' @return A data.frame with the scrublet results
runSoupX = function(datadir,toc = NULL,metadata=NULL,plotDir=NULL,scPath=scPath){
  
  if(is.null(metadata)){
    message('Please provide metadata!')
  }
  # Read in table of counts
  if(is.null(toc)){
    toc = Seurat::Read10X(file.path(datadir, "filtered_feature_bc_matrices"))
  }
  # Read in table of droplets
  tod = Seurat::Read10X(file.path(datadir, "raw_feature_bc_matrix"))
  rownames(tod) = gsub('_','-',rownames(tod))
  
  sc = SoupChannel(tod, toc)
  
  # Add clustering metadata
  sc = setClusters(sc, setNames(metadata$seurat_clusters, rownames(metadata)))
  sc = setDR(sc, metadata[colnames(sc$toc), c("RD1", "RD2")])
  
  if(verbose)
    message("Making diagnostic plots.")
  
  #Make diagnostic plots!
  if(!is.null(plotDir)){
    pdf(file.path(paste0(plotDir,names(datadir),'_soupXplots.pdf')),width = 9.5, height = 8, paper = 'a4r') 
  }
  # Visual sanity checks
  require(ggplot2)
  dd = metadata[colnames(sc$toc), ]
  mids = aggregate(cbind(RD1, RD2) ~ seurat_clusters, data = dd, FUN = mean)
  gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = seurat_clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = seurat_clusters)) + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
  plot(gg)
  
  gg = plotMarkerMap(sc, c("HBA1",'HBB','HBA2'))
  plot(gg)
  
  
  
  if(names(datadir) %in% c('adr66')){
    sc = setContaminationFraction(sc, 0.1)
    out = adjustCounts(sc)
    #  out = toc
    #}else if (names(datadir) == 'adr36'){
    #  sc = autoEstCont(sc,tfidfMin = 0.7)  
    # Clean the data
    #  out = adjustCounts(sc)
  }else if(names(datadir) == 'CG.SB.NB8715409'){
    # Estimate contamination/soup fraction
    sc = autoEstCont(sc,tfidfMin = 0.7)  
    # Clean the data
    out = adjustCounts(sc)
  }else if(names(datadir) == 'CG.SB.NB8715411'){
    # Estimate contamination/soup fraction
    sc = autoEstCont(sc,tfidfMin = 0.8)  
    # Clean the data
    out = adjustCounts(sc)
  }else{
    # Estimate contamination/soup fraction
    sc = autoEstCont(sc)  
    # Clean the data
    out = adjustCounts(sc)
  }
  
  
  
  # Investigating changes in expression
  cntSoggy = rowSums(sc$toc > 0)
  cntStrained = rowSums(out > 0)
  mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
  mostZeroed
  tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
  
  # Visualising expression distribution
  plotChangeMap(sc, out, "MT-ND3")
  if(is.null(scPath)){
    DropletUtils:::write10xCounts(file.path(datadir,"cleanCounts/strainedCounts"), out,overwrite = T)  
  }else{
    DropletUtils:::write10xCounts(file.path(datadir,gsub('^../','',scPath)), out,overwrite = T)  
  }
  
  
  srat.out = CreateSeuratObject(out)
  # Add metadata to srat.out 
  metadata$cellID = rownames(metadata)
  srat.out@meta.data$cellID = rownames(srat.out@meta.data)
  m = match(srat.out@meta.data$cellID,metadata$cellID)
  sum(is.na(m))
  srat.out@meta.data =  merge(srat.out@meta.data,metadata[,!colnames(metadata) %in% colnames(srat.out@meta.data)],by=0,all.x = T)
  rownames(srat.out@meta.data) =  srat.out@meta.data$cellID
  
  if(!is.null(plotDir)){
    dev.off()
  }
  
  return(srat.out)
}



#' Initial filtering and metadata
#'
#' Calculates metadata, drop cells based on QC metrics.  Doesn't require the pre-processing structure of loading things into unswappedCounts and cleanCounts, along with the corresponding automatically generated files.  But if these are given, the extra files and meta-data will be used and loaded into the resulting Seurat object.
#'
#' @param dataDirs Data directories to load, passed to Read10X.
#' @param excludeGenes Genes that will be dropped from the count matrix.
#' @param geneSets Sets of genes to aggregate counts for and store as metadata.  Name by list element names.  Must include an entry named "mtGenes".
#' @param sPhaseGenes Genes used to identify S phase genes.  If NULL, will be loaded from Seurat.
#' @param g2mPhaseGenes Genes used to identify G2M phase genes.  If NULL, will be loaded from Seurat.
#' @param maxMT Maximum fraction of expression from Mitochondrial genes allow per cell.
#' @param minGenes Minimum number of genes for a cell to express and pass QC.
#' @param minUMIs Minimum number of UMIs for a cell to have and pass QC.
#' @param maxBadFrac Maximum fraction of cells not passing QC in a cluster before we dump the whole cluster.
#' @param skipScrub Should we skip running scrublet, even if it's not available.
#' @param skipSoup Should we skip running scrublet, even if it's not available.
#' @param scrubScoreMax Mark a cell as a doublet if its scrublet score exceeds this value, even if scrublet does not call it as such.
#' @param numPCs Number of PCs to use for the rough clustering and doublet identification.  This is just used for generating very fine grained clustering and it shouldn't matter much what you set it to.
#' @param clusteringRes Resolution for fine clustering.  Should be very large.
#' @param scrubPath Relative path to file containing scrublet run.  If NULL, scrublet will be run on data. 
#' @param scPath Path to where SoupChannel objects created by SoupX are saved in RDS format.  If not found, ignored.
#' @param cellCallPath Path to file detailing the reason each barcode is called as a cell. 
#' @param doPlot Make some QC summary plots.
#' @param verbose Shut up or not?
#' @param ... Parameters passed to \code{Read10X}
#' @return A Seurat object with cells dropped that fail the QC filters.
basicQC = function(dataDirs,maxMT=20,minGenes=300,minUMIs=1000,maxBadFrac=0.5,numPCs=75,clusteringRes=10.0,
                   skipScrub=FALSE,skipSoup=FALSE,scrubScoreMax = 0.5,scrubPath='../cleanCounts/scrubletScores.tsv',
                   metadata=NULL,matchBy=NULL,scPath="../cleanCounts/strainedCounts",outPath,skipIfExists=F,keepMTCells=T,
                   excludeGenes = c(),geneSets=list(),mtPattern = '^MT-',sPhaseGenes=NULL,g2mPhaseGenes=NULL,doPlot=TRUE,plotDir=NULL,verbose=TRUE,...){
  require(tidyverse)
  require(Seurat)
  require(ggVennDiagram)
  require(RColorBrewer)
  
  # Check if plotDir is provided
  if(doPlot & is.null(plotDir)){
    warning('Please provide plotDir or set doPlot=F')
    return()
  }
  if(is.null(dataDirs)){
    warning('Please provide dataDirs')
    return()
  }
  
  if(!is.null(outPath) & file.exists(paste0(outPath,'clean_noMTCells.RDS')) & skipIfExists==T){
    # Read in existing output
    message(sprintf('Existing clean output detected:\noutPath = %s \nTo overwite existing file, please specify skipIfExists = F',paste0(outPath,'clean_noMTCells.RDS')))
    
    if(keepMTCells & file.exists(paste0(outPath,'clean_withMTCells.RDS'))){
      message(sprintf('Existing clean output WITH MT Cells detected:\noutPath = %s \nReading this instead of clean_noMTCells output',paste0(outPath,'clean_withMTCells.RDS')))
      srat = readRDS(paste0(outPath,'clean_withMTCells.RDS'))
    }else if(keepMTCells & !file.exists(paste0(outPath,'clean_withMTCells.RDS'))){
      message(sprintf('NO clean output WITH MT Cells detected. Run analysis...\nTo import existing clean_noMTCells.RDS, please specify keepMTCells = F'))
      srat = readRDS(paste0(outPath,'clean_noMTCells.RDS'))
    }else if(!keepMTCells){
      srat = readRDS(paste0(outPath,'clean_noMTCells.RDS'))
    }
    #srat = read.RDS(outPath)
    
    df = srat@misc$preQC@meta.data
    df.out = df %>% group_by(orig.ident) %>% summarise(nCells_prefilter=n(),medMT=median(percent.mt),muMT=mean(percent.mt),nPASS = sum(PASS),nDoublet = sum(PASS_doublet==FALSE),medCnt = median(nCount_RNA),medGenes=median(nFeature_RNA),nBadClust = sum(PASS_cluster==FALSE),
                                                       lowCnt = sum(PASS_nCounts==F),highMT = sum(PASS_MT==F),lowGenes = sum(PASS_nGenes==F))
    return(list(srat,df.out))
  }
  
  if(is.null(names(dataDirs))){
    warning("Data directories not named, generating automatic names.")
    names(dataDirs) = paste0('DataSource',seq_along(dataDirs))
  }
  
  if(verbose)
    message('Loading data.')
  srat = Read10X(dataDirs,...)
  srat = CreateSeuratObject(srat)
  
  #### =================================================== #
  # Load and store the gene objects
  if(verbose)
    message("Loading gene meta-data.")
  srat@misc$geneMap = lapply(dataDirs,function(e) {
    if(file.exists(file.path(e,'features.tsv.gz'))){
      x = read.table(file.path(e,'features.tsv.gz'),sep='\t',header=FALSE)
    }else{
      x = read.table(file.path(e,'genes.tsv'),sep='\t',header=FALSE)
    }
    #Do the Seurat gene name conversion.  Assume column 2.
    x[,2] = make.unique(gsub('_','-',x[,2]))
    return(x)
  })
  
  #Check if they're all the same.
  srat@misc$geneMapCommon=TRUE
  for(i in seq_along(dataDirs)){
    if(i==1){
      base = srat@misc$geneMap[[i]]
    }else{
      tst = srat@misc$geneMap[[i]]
      #Check they have the same dimensions
      if(!all(dim(base)==dim(tst))){
        srat@misc$geneMapCommon=FALSE
        next
      }
      #Check they're all the same
      for(j in seq(ncol(base)))
        if(!all(base[,j]==tst[,j]))
          srat@misc$geneMapCommon=FALSE
    }
  }
  #If they're all the same, store the same one
  if(srat@misc$geneMapCommon){
    srat@misc$geneMapCommon = base
    rownames(srat@misc$geneMapCommon) = base[,2]
  }else{
    srat@misc$geneMapCommon = NULL
  }
  
  
  
  #### =================================================== #
  srat = NormalizeData(srat,verbose=FALSE)
  # Get cell cycle scores
  data('cc.genes.updated.2019',package='Seurat')
  sPhaseGenes = cc.genes.updated.2019$s.genes
  g2mPhaseGenes = cc.genes.updated.2019$g2m.genes
  #The seed parameter magic is needed because the Seurat authors are dicks
  srat = CellCycleScoring(srat, s.features=sPhaseGenes,g2m.features=g2mPhaseGenes,seed=sample(1e9,1))
  
  # Get MT content
  if(!is.null(mtPattern)){
    srat[['percent.mt']] = PercentageFeatureSet(srat, pattern = mtPattern)
  }
  
  # Get content of any other gene sets
  if(is.null(names(geneSets))){
    if(length(geneSets)==0){
      message('No extra gene sets provided.')
    }else if (length(geneSets) > 0){
      warning("No names given for gene sets.  Setting to uninformative auto-generated values.")
      names(geneSets) = paste0('geneSet',seq_along(geneSets))  
    }
    
  }
  if(verbose)
    message("Generating extra meta-data.")
  for(nom in names(geneSets)){
    tmp = geneSets[[nom]]
    tmp = tmp[tmp %in% rownames(srat@assays$RNA@data)]
    if(length(tmp)==0)
      next
    srat@meta.data[,nom] = colSums(srat@assays$RNA@counts[tmp,,drop=FALSE])/srat@meta.data$nCount_RNA
  }
  
  
  #### =================================================== #
  # Call doublets with Scrublet
  if(verbose)
    message("Loading or running scrublet.")
  srat@meta.data$scrubScore=NA
  srat@meta.data$scrubCall=NA
  srat@meta.data$scrubSim1=NA
  srat@meta.data$scrubSim2=NA
  if(skipScrub){
    scrubRun = rep(FALSE,length(dataDirs))
  }else{
    if(!is.null(scrubPath)){
      scrubRun = !file.exists(file.path(dataDirs,scrubPath))
    }else{
      scrubRun = rep(TRUE,length(dataDirs))
    }
    
    for(i in seq_along(dataDirs)){
      if(scrubRun[i]){
        if(verbose)
          message(sprintf('Running for channel %s',names(dataDirs)[i]))
        setwd(dataDirs[i])
        if(!file.exists(file.path(dataDirs[i],'../cleanCounts'))){
          if(verbose)
            message(sprintf('Creating Clean Count dir for channel %s',names(dataDirs)[i]))
          dir.create(file.path(dataDirs[i],'../cleanCounts'))
        }
        
        scrubScores = runScrublet(dataDirs[i],nPCs=numPCs)
        write_delim(scrubScores,file.path(dataDirs[i],'../cleanCounts/scrubletScores.tsv'),delim = '\t')
      }else{
        scrubScores = read.table(file.path(dataDirs[i],scrubPath),sep='\t',header=TRUE)
      }
      
      #Save the relevant information in meta.data
      cellIDs = rownames(srat@meta.data)
      if(sum(grepl('-\\d$',cellIDs))>0 & sum(grepl('-\\d$',scrubScores$barcode))==0){
        cellIDs = gsub('-\\d$','',cellIDs)
      }
      
      m = match(paste0(names(dataDirs)[i],'_',scrubScores$barcode),cellIDs)
      o = !is.na(m)
      if(sum(o) == 0){
        m = match(scrubScores$barcode,cellIDs)
        o = !is.na(m)
      }
      if(sum(o) == 0){
        warning(sprintf('Skipping Scrublet for channel %s as scrublet cell barcodes do not match with seurat object. Please check!',names(dataDirs)[i]))
      }else{
        srat@meta.data$scrubScore[m[o]] = scrubScores$score[o]
        srat@meta.data$scrubCall[m[o]] = scrubScores$call[o]
        srat@meta.data$scrubSim1[m[o]] = scrubScores$simScore1[o]
        srat@meta.data$scrubSim2[m[o]] = scrubScores$simScore2[o]
      }
    }
  }
  
  
  #### =================================================== #
  #Do basic processing and clustering
  if(verbose)
    message("Basic processing to complete QC.")
  srat = FindVariableFeatures(srat,verbose=FALSE)
  srat = ScaleData(srat,verbose=FALSE)
  srat = RunPCA(srat,npcs=numPCs,approx=FALSE,verbose=FALSE)
  srat = FindNeighbors(srat,dims=seq(numPCs),verbose=FALSE)
  srat = FindClusters(srat,res=clusteringRes,verbose=FALSE)
  srat = RunUMAP(srat,dims=seq(numPCs),verbose=FALSE)
  
  #### =================================================== #
  #Set the basic filters. Default to TRUE if NA
  #srat@meta.data$PASS_MT = !(srat@meta.data$percent.mt >= maxMT)
  srat@meta.data$PASS_MT = TRUE
  srat@meta.data$PASS_nGenes = !(srat@meta.data$nFeature_RNA <= minGenes)
  srat@meta.data$PASS_nCounts = !(srat@meta.data$nCount_RNA <= minUMIs)
  srat@meta.data$PASS_doublet = !(!is.na(srat@meta.data$scrubCall) & (srat@meta.data$scrubCall | srat@meta.data$scrubScore > scrubScoreMax))
  table(srat@meta.data$PASS_doublet)
  #Work out which ones don't pass because of clustering
  # this gives Percentage of cells in each clusters that did not pass at least one of these criteria
  tt = aggregate(!PASS_MT | !PASS_nGenes | !PASS_nCounts | !PASS_doublet ~ seurat_clusters,FUN=mean,data=srat@meta.data) 
  srat@meta.data$fracBadClust = tt[match(srat@meta.data$seurat_clusters,tt[,1]),2]
  # Removing clusters containing too many bad cells
  srat@meta.data$PASS_cluster = !(srat@meta.data$fracBadClust > maxBadFrac)
  #Work out the final filter
  srat@meta.data$PASS = with(srat@meta.data,PASS_MT & PASS_nGenes & PASS_nCounts & PASS_doublet & PASS_cluster)
  #Make a reason for fail variable
  nBasicFail = with(srat@meta.data,4 - (PASS_MT + PASS_nGenes + PASS_nCounts + PASS_doublet))
  tmp = rep('Pass',length(nBasicFail))
  tmp[!srat@meta.data$PASS_MT] = 'highMT'
  tmp[!srat@meta.data$PASS_nGenes] = '#genesLow'
  tmp[!srat@meta.data$PASS_nCounts] = '#countsLow'
  tmp[!srat@meta.data$PASS_doublet] = 'doublet'
  tmp[nBasicFail>1] = 'multiple'
  #Don't want "multiple" to encompass cluster otherwise there'll be heaps where cluster + bad obscures the true reason
  tmp[tmp=='Pass' & !srat@meta.data$PASS_cluster] = 'badCluster'
  srat@meta.data$reasonForFail = tmp
  
  #### =================================================== #
  #### Run SoupX to remove ambiant reads/mRNA
  if(!skipSoup){
    preSoupSrat = subset(srat,subset = PASS)
    soupedSrat = cleanWithSoupX(cleanSrat = preSoupSrat,dataDirs,scPath = scPath,plotDir = plotDir) 
    
    # Add metadata to srat
    srat@meta.data$PASS_soupX = F
    m=match(rownames(soupedSrat@meta.data),rownames(srat@meta.data))
    srat@meta.data$PASS_soupX[m] = !(soupedSrat@meta.data$nCount_RNA <= minUMIs)
    srat@meta.data$postSoupX_nCount=srat@meta.data$nCount_RNA
    srat@meta.data$postSoupX_nCount[m] = soupedSrat@meta.data$nCount_RNA
    
    m = which((!srat@meta.data$PASS_soupX) & (srat@meta.data$PASS))
    srat@meta.data$PASS[m] = F
    srat@meta.data$reasonForFail[m] = 'lowcount_postSoupX'
    
  }else{
    soupedSrat = srat
  }
  
  # Add Real MT PASS filter
  srat@meta.data$PASS_withMT = srat@meta.data$PASS
  
  srat@meta.data$PASS_MT = !(srat@meta.data$percent.mt >= maxMT)
  m = which((!srat@meta.data$PASS_MT) & (srat@meta.data$PASS))
  srat@meta.data$PASS[m] = F
  srat@meta.data$reasonForFail[m] = 'highMT'
  
  #Move PASS column to the end
  w = which(colnames(srat@meta.data)=='PASS')
  srat@meta.data = srat@meta.data[,c(1:(w-1),(w+1):ncol(srat@meta.data),w)]
  
  
  # Add metadata if provided
  srat[['cellID']] = rownames(srat@meta.data)
  
  if('HCA' %in% unique(srat@meta.data$orig.ident)){
    srat$orig.ident = paste0(srat$orig.ident,'_',sapply(strsplit(srat$cellID,split = '_'),'[',2))
  }
  
  if(!is.null(metadata)){
    if(!is.null(matchBy) & length(matchBy) == 2){
      m = match(srat[[matchBy[1]]],metadata[[matchBy[2]]])
      o = !is.na(m)
    }else{
      # Default is to match metadata by orig.ident
      m = match(srat@meta.data$orig.ident,metadata$orig.ident)
      o = !is.na(m)
    }
    
    if(sum(o)==0){
      warning('No matching found between seurat object and metadata provided. Please check!')
    }else{
      srat@meta.data = merge(srat@meta.data,metadata,by='orig.ident',all.x=T)
      rownames(srat@meta.data) = srat@meta.data$cellID
    }
  }
  
  
  #### =================================================== #
  #### Filter
  if(verbose)
    message("Filtering and finalizing.")
  goodCells = rownames(srat@meta.data[srat@meta.data$PASS_withMT,])
  w = which(srat@meta.data$PASS_withMT)
  #Make the new version, store the old one in it
  srat@meta.data$preSoupX_nCount = srat@meta.data$nCount_RNA
  sratOld = srat
  #Restart without them
  mDat = srat@meta.data[w,]
  #Drop old clustering
  mDat = mDat[,!colnames(mDat) %in% c('seurat_clusters',grep('^RNA_snn_res',colnames(mDat),value=TRUE))]
  #And drop genes that we don't care about
  genesToKeep = rownames(srat@assays$RNA@counts)
  genesToKeep = genesToKeep[!genesToKeep %in% excludeGenes]
  srat = CreateSeuratObject(soupedSrat@assays$RNA@counts[genesToKeep,colnames(soupedSrat@assays$RNA@counts) %in% goodCells])
  #Merge in metadata
  m=match(rownames(srat@meta.data),rownames(mDat))
  srat@meta.data = cbind(srat@meta.data,mDat[m,!(colnames(mDat) %in% colnames(srat@meta.data)),drop=FALSE])
  
  srat = NormalizeData(srat,verbose=FALSE)
  srat = FindVariableFeatures(srat,verbose=FALSE)
  srat@misc$preQC = sratOld
  srat@misc$geneMap = sratOld@misc$geneMap
  srat@misc$geneMapCommon = sratOld@misc$geneMapCommon
  
  
  #### =================================================== #
  #Make some QC plots before dropping things
  if(doPlot){
    if(verbose)
      message("Making diagnostic plots.")
    
    #Plot distributions of basic things
    if(!is.null(plotDir)){
      pdf(file.path(paste0(plotDir,'QCplots.pdf')),width = 9.5, height = 8, paper = 'a4r') 
    }
    
    df = sratOld@meta.data
    df.out = df %>% group_by(orig.ident) %>% summarise(nCells_prefilter=n(),medMT=median(percent.mt),muMT=mean(percent.mt),nPASS = sum(PASS),nDoublet = sum(PASS_doublet==FALSE),medCnt = median(nCount_RNA),medGenes=median(nFeature_RNA),nBadClust = sum(PASS_cluster==FALSE),
                                                       lowCnt = sum(PASS_nCounts==F),highMT = sum(PASS_MT==F),lowGenes = sum(PASS_nGenes==F))
    
    ## QC Part 1: How many cells prefiltered? ##
    counts_preprocess = df %>% group_by(orig.ident) %>% summarise(nCells_prefilter = n())
    p = ggplot(data=counts_preprocess, aes(x=orig.ident, y=nCells_prefilter)) + 
      geom_bar(stat="identity") + 
      #scale_fill_manual(values = c("#85a040","#8774ca","#ca7040","#4bae8d","#ca5688"))+
      theme_bw() + 
      theme(text=element_text(size=11), axis.text.x=element_text(angle=45,hjust = 1)) + 
      ggtitle('Prefiltered Counts') + xlab('') + ylab('# cells')
    plot(p)
    #ggsave(file="/home/jovyan/Aneuploidy/Dec21/adrenal/prefiltered_counts.pdf", dpi=300, height=5, width=7)
    
    #QC Part 2: Mitochondrial Content ##
    p=ggplot(data=df, aes(x=percent.mt)) + 
      geom_histogram(bins=500,color='black',fill='white') +
      facet_grid(orig.ident~., scales="free") + 
      #scale_fill_manual(values = c("#85a040","#8774ca","#ca7040","#4bae8d","#ca5688"))+
      theme_bw() + theme(text=element_text(size=10)) + xlab('')+
      geom_vline(xintercept=maxMT, linetype="dashed", color="red")+
      ggtitle('MT content')
    plot(p)
    #ggsave(file="/home/jovyan/Aneuploidy/Dec21/adrenal/MT_Distribution.pdf", dpi=300, height=10, width=15)
    
    #QC Part 3: nFeature and nUMI Counts ##
    p=ggplot(data=df, aes(x=log10(nCount_RNA), y=nFeature_RNA, color=percent.mt)) + 
      geom_point(stat="identity",size=0.2) + 
      geom_hline(yintercept = minGenes, linetype="dashed", color="black") + 
      geom_vline(xintercept = log10(minUMIs), linetype="dashed", color="black") + 
      facet_grid(.~orig.ident) +
      theme_bw() + 
      scale_y_log10()+
      theme(text=element_text(size=15)) + 
      #geom_hline(yintercept=9000, linetype="dashed", color="red") + 
      scale_color_gradient(high="#cc2b5e", low="grey") + 
      ggtitle("UMI Count vs Features")
    plot(p)
    #ggsave(file="/home/jovyan/Aneuploidy/Dec21/adrenal/nFeature_vs_nUMI.pdf", dpi=300, height=10, width=15)
    
    p=ggplot(df,aes(x = scrubScore))+
      geom_histogram(bins=300,color='black',fill='white') +
      facet_grid(orig.ident~., scales="free") + 
      #scale_fill_manual(values = c("#85a040","#8774ca","#ca7040","#4bae8d","#ca5688"))+
      theme_bw() + theme(text=element_text(size=10)) + xlab('')+
      geom_vline(xintercept=scrubScoreMax, linetype="dashed", color="red")+
      ggtitle('Scrub Score')
    plot(p)
    
    #QC Part 4: How many cells are being removed from each filter? ##
    nCount_removed <- rownames(df[df$nCount_RNA <= minUMIs,]) #1
    nFeature_removed <- rownames(df[df$nFeature_RNA <= minGenes,]) #2
    MT_Content_removed <- rownames(df[df$percent.mt > maxMT,]) #3
    #badCluster <- rownames(df[df$PASS_cluster==F,]) #4
    #doublets <- rownames(df[df$PASS_doublet==F,]) #4
    
    x <- list(nCount = nCount_removed, nFeature = nFeature_removed, MT_Content = MT_Content_removed)#,
    #badClusters = badCluster, doublets = doublets)
    #pdf(file="/home/jovyan/Aneuploidy/Dec21/adrenal/cellsRemoved_Venn.pdf", height=6, width=8)
    p=ggVennDiagram(x, label_size=4, color="black") + scale_fill_gradient(low = "white", high = "white")
    plot(p)
    #dev.off()
    
    ## QC Part 5: Plot UMI count before (post Scrublet) and after SoupX ##
    df2 = sratOld@meta.data
    p=ggplot(data=df2, aes(x=log10(nCount_RNA), y=log10(preSoupX_nCount))) + 
      geom_point(stat="identity",size=0.01) + 
      geom_hline(yintercept = log10(minUMIs), linetype="dashed", color="red",size=0.1) + 
      geom_vline(xintercept = log10(minUMIs), linetype="dashed", color="red",size=0.1) + 
      facet_grid(.~orig.ident) +
      theme_bw() + 
      theme(text=element_text(size=15)) + 
      #geom_hline(yintercept=9000, linetype="dashed", color="red") + 
      scale_color_gradient(high="#cc2b5e", low="grey") + 
      ggtitle("UMI count pre/post SoupX")
    plot(p)
    
    ## QC Part 6: How many cells pre vs post-filtered, Scrublet and SoupX? ##
    counts_postprocess = srat@meta.data %>% group_by(orig.ident) %>% summarise(nCells_postfilter = n())
    counts_postprocess = merge(counts_postprocess,counts_preprocess,by='orig.ident')
    counts_postprocess = pivot_longer(counts_postprocess,cols = c(2:3),names_to = 'type',values_to = 'nCells')
    counts_postprocess$type = factor(counts_postprocess$type,levels = c('nCells_prefilter','nCells_postfilter'))
    p = ggplot(data=counts_postprocess, aes(x=orig.ident, y=nCells,fill=type)) + 
      geom_bar(stat="identity",position = 'dodge') + 
      scale_fill_manual(values = c("darkgrey","black"))+
      theme_bw() + 
      theme(text=element_text(size=11), axis.text.x=element_text(angle=45,hjust = 1)) + 
      ggtitle('Postfiltered Cell Counts') + xlab('') + ylab('# cells')
    plot(p)
    
    #Standard QC things
    plot(DimPlot(sratOld,group.by='seurat_clusters',label=TRUE,repel=TRUE)+guides(colour="none"))
    
    #A hack for another of Seurat's stupid decisions
    if('labels' %in% names(as.list(args(DimPlot)))){
      plot(DimPlot(sratOld,group.by=c('Phase','orig.ident'),labels=c('Phase','Source')))
    }else{
      gg = DimPlot(sratOld,group.by=c('Phase','orig.ident'),combine=FALSE,cols = c(brewer.pal(12,'Paired'),brewer.pal(12,'Dark2'),brewer.pal(12,'Set1')))
      labs = c('Phase','Source')
      for(i in seq_along(labs)){
        gg[[i]] = gg[[i]]  + ggtitle(labs[i])
        plot(gg[[i]])
      }
    }
    plot(FeaturePlot(sratOld,c('percent.mt','nFeature_RNA','nCount_RNA','fracBadClust')))
    #Which cells pass based on what
    if('labels' %in% names(as.list(args(DimPlot)))){
      plot(DimPlot(sratOld,group.by=c('PASS_MT','PASS_nGenes','PASS_nCounts','PASS_doublet','PASS_cluster','PASS_soupX','PASS'),labels=c('Pass MT Frac?','Pass #Genes?','Pass #Counts?','Pass Doublet?','Pass cluster','Pass soupX?','PASS'),legend='none'))
    }else{
      gg = DimPlot(sratOld,group.by=c('PASS_MT','PASS_nGenes','PASS_nCounts','PASS_doublet','PASS_cluster','PASS'),combine=FALSE)
      labs =c('Pass MT Frac?','Pass #Genes?','Pass #Counts?','Pass Doublet?','Pass cluster?','PASS?')
      for(i in seq_along(labs)){
        gg[[i]] = gg[[i]] +guides(colour="none") + ggtitle(labs[i])
      }
      plot(patchwork::wrap_plots(gg))
    }
    plot(DimPlot(sratOld,group.by=c('reasonForFail'),cols = c(brewer.pal(9,'Paired')[-c(1,6)],'lightgrey')))
    
    if(!is.null(plotDir)){
      dev.off()
    }
  }
  
  
  if(!is.null(outPath)){
    if(keepMTCells){
      # Save the results with all highMT cells
      message(sprintf('Saving cleanSrat with hightMT cells to %s',outPath))
      saveRDS(srat,paste0(outPath,'clean_withMTCells.RDS'))   
    }
    
    # Save the results without highMT cells
    message(sprintf('Saving cleanSrat to %s',outPath))
    print(dim(srat))
    print(table(srat$PASS))
    srat = subset(srat, subset = PASS)
    print(dim(srat))
    print(table(srat$PASS))
    srat = NormalizeData(srat,verbose=FALSE)
    srat = FindVariableFeatures(srat,verbose=FALSE)
    saveRDS(srat,paste0(outPath,'clean_noMTCells.RDS')) 
  }
  
  return(list(srat,df.out))
}







