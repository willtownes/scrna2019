#functions for filtering genes by various informativeness criteria
library(SingleCellExperiment)
library(Seurat)
source("./util/functions.R") #needed for compute_gene_info

filterDev<-function(sce,nkeep=nrow(sce),dev=c("binomial","multinomial","poisson","geometric")){
  dev<-match.arg(dev)
  gm<-compute_gene_info(counts(sce),gmeta=rowData(sce),mod=dev)
  o<-order(gm$deviance,decreasing=TRUE,na.last=FALSE)
  #NA deviance => badly fitting null model=> highly variable gene
  res<-sce[o[1:nkeep],]
  res[,colSums(counts(res))>0]
}

filterExpr<-function(sce, nkeep=nrow(sce)){
  #modified from https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/filtering/filterExpr.R
  exprsn<-rowMeans(logcounts(sce))
  o<-order(exprsn,decreasing=TRUE)
  res<-sce[o[1:nkeep],]
  res[,colSums(counts(res))>0]
}

filterHVG<-function(sce, nkeep=nrow(sce)){
  #modified from https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/filtering/filterHVG.R
  seurat<-CreateSeuratObject(counts(sce),meta.data=as.data.frame(colData(sce)),min.cells=0,min.genes=0,project = "scRNAseq")
  seurat<-NormalizeData(seurat,display.progress=FALSE)
  seurat<-ScaleData(seurat,vars.to.regress="nUMI",display.progress=FALSE)
  seurat<-FindVariableGenes(seurat,do.plot=FALSE,sort.results=TRUE,selection.method="dispersion",top.genes=nkeep,display.progress=FALSE)
  res<-sce[seurat@var.genes, ]
  res[,colSums(counts(res))>0]
}