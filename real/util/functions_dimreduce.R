#install.packages("glmpca")
library(glmpca)
source("./util/functions.R")
source("./algs/existing.R")

logcpm<-function(Y){
  Ycpm<-1e6*t(t(Y)/colSums(Y))
  log2(1+Ycpm)
}

resd<-function(Y){
  null_residuals(Y,type="deviance",mod="binomial")
}

resp<-function(Y){
  null_residuals(Y,type="pearson",mod="binomial")
}

dimreduce_all<-function(Y,L=2,cm=NULL,...){
  #Y a counts matrix with rows=genes, cols=cells
  #L=number of latent dimensions
  #cm a data.frame with column metadata, eg prob of zero, cluster IDs for each cell
  #... additional args passed to glmpca, for example, increase L2 penalty
  res<-list()
  tt<-list()
  tt$pca_log<-system.time(res$pca_log<-pca(logcpm(Y),L)) 
  tt$pca_rd<-system.time(res$pca_rd<-pca(resd(Y),L))
  tt$pca_rp<-system.time(res$pca_rp<-pca(resp(Y),L))
  # tt$glmpca_poi<-system.time(
  #   res$glmpca_poi<-tryCatch(
  #     glmpca(Y,L,fam="poi")$factors,
  #     error=function(e){
  #       d<-rep(NA,L)
  #       names(d)<-paste0("dim",seq_along(d))
  #       as.data.frame(as.list(d))
  #     }
  #   )
  # )
  tt$glmpca_poi<-system.time(res$glmpca_poi<-glmpca(Y,L,fam="poi",...)$factors)
  tt$glmpca_nb<-system.time(res$glmpca_nb<-glmpca(Y,L,fam="nb",...)$factors)
  tt$zinbwave<-system.time(res$zinbwave<-zinbwave(Y,L))
  for(i in names(res)){
    res[[i]]$cell_id<-rownames(res[[i]])
    res[[i]]$dimreduce<-i
    if(!is.null(cm)){ res[[i]]<-cbind(res[[i]],cm) }
  }
  pd<-do.call(rbind,res)
  rownames(pd)<-NULL
  tt<-vapply(tt,function(t){t[["elapsed"]]},FUN.VALUE=1.0)
  list(factors=pd,elapsed=tt)
}

tfunc<-function(tt){
  #helper function for elapsed time extraction (see format_elapsed)
  data.frame(dimreduce=names(tt),elapsed_sec=as.numeric(tt))
}

format_elapsed<-function(elapsed_dev,elapsed_hvg){
  #convenience function for formatting the elapsed times of
  #deviance and hvg filtered dimreduce_all runs
  #returns a dataframe
  tt1<-tfunc(elapsed_dev)
  tt1$genefilter<-"dev"
  tt2<-tfunc(elapsed_hvg)
  tt2$genefilter<-"hvg"
  tt<-rbind(tt1,tt2)
  tt[,c(3,1,2)]
}
