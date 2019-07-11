#functions for filtering genes by various informativeness criteria
suppressPackageStartupMessages(library(SingleCellExperiment))
library(Seurat)
source("./util/functions.R") #needed for compute_gene_info

filterDev<-function(sce,nkeep=nrow(sce),dev=c("binomial","poisson","geometric"),ret=c("sce","ranks")){
  dev<-match.arg(dev)
  ret<-match.arg(ret)
  gm<-compute_gene_info(counts(sce),gmeta=rowData(sce),mod=dev)
  o<-order(gm$deviance,decreasing=TRUE,na.last=FALSE)[1:nkeep]
  #NA deviance => badly fitting null model=> highly variable gene
  if(ret=="sce"){
    res<-sce[o,]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

filterExpr<-function(sce, nkeep=nrow(sce),ret=c("sce","ranks")){
  #modified from https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/filtering/filterExpr.R
  ret<-match.arg(ret)
  exprsn<-rowMeans(logcounts(sce))
  o<-order(exprsn,decreasing=TRUE)[1:nkeep]
  if(ret=="sce"){
    res<-sce[o,]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

filterHVG<-function(sce, nkeep=nrow(sce), total_umi="nUMI", ret=c("sce","ranks")){
  #modified from https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/filtering/filterHVG.R
  ret<-match.arg(ret)
  if(!(total_umi %in% colnames(colData(sce)))){ stop("total_umi must be in coldata") }
  seu<-CreateSeuratObject(counts(sce),meta.data=as.data.frame(colData(sce)),min.cells=0,min.features=0,project="scRNAseq")
  seu<-NormalizeData(seu,verbose=FALSE)
  seu<-ScaleData(seu,vars.to.regress=total_umi,verbose=FALSE)
  seu<-FindVariableFeatures(seu,selection.method="dispersion",nfeatures=nkeep,verbose=FALSE)
  vf<-head(VariableFeatures(seu), nkeep)
  #sometimes Seurat coerces rownames, so use numeric IDs instead
  o<-match(vf,rownames(seu))
  if(ret=="sce"){
    res<-sce[o, ]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

rank_all_genes<-function(sce, total_umi="nUMI"){
  #sce=SingleCellExperiment object with assays "counts" and "logcounts"
  #returns a dataframe with same rownames as sce
  #columns: dev,hvg,expr
  #each column is the rank order of genes by each criteria
  #low rank=highly informative gene
  gg<-rownames(sce)
  gfs<-list()
  gfs$dev<-filterDev(sce,ret="ranks")
  gfs$expr<-filterExpr(sce,ret="ranks")
  gfs$hvg<-filterHVG(sce,total_umi=total_umi,ret="ranks")
  res<-list()
  for(i in names(gfs)){
    rk<-seq_along(gfs[[i]])
    names(rk)<-gfs[[i]]
    res[[i]]<-rk[gg]
  }
  res<-as.data.frame(res)
  rownames(res)<-gg
  res
}

meanvar_plotdata<-function(sce,G=1000){
  #sce=SingleCellExperiment with assays "counts" and "logcounts"
  #logcounts = log2(1+counts/scran size factor)
  #G=number of genes to be "highly informative"
  #assumes rowData(sce) contains columns "dev","hvg" with ranks
  #rank=1 means gene is most informative according to the criterion
  #returns a data frame to use for making a mean/variance plot
  gm<-as.data.frame(rowData(sce))
  rk<-gm[,c("dev","expr","hvg")]
  rk<-rk[rownames(sce),]
  Y<-2^(logcounts(sce))-1
  pd<-data.frame(m=rowMeans(Y),v=apply(Y,1,var))
  pd$vmr=pd$v/pd$m
  pd<-cbind(pd,rk)
  pd<-subset(pd,vmr>0)
  xt<-table(pd$dev<=G,pd$hvg<=G)
  #formatting plot labels function
  f<-function(s,n){paste0(s," (",prettyNum(n,big.mark=","),")")}
  pd$criteria<-f("neither",xt[1,1])
  pd$criteria[pd$dev<=G & pd$hvg<=G]<-f("both",xt[2,2])
  pd$criteria[pd$dev<=G & pd$hvg>G]<-f("high deviance",xt[2,1])
  pd$criteria[pd$dev>G & pd$hvg<=G]<-f("highly variable",xt[1,2])
  #sort so that the "neither" category doesn't obscure the "highly deviant"
  pd[order(pd$criteria,decreasing=TRUE),]
}
