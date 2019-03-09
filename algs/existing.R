# Functions for running existing methods such as tSNE, PCA etc
library(Matrix)

rm_zero_rowcol<-function(Y){
  #remove all rows and columns containing all zeros
  Y<-Y[rowSums(Y>0)>0,] #remove rows with zeros all the way across
  Y<-Y[,colSums(Y>0)>0]
  Y
}

norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}

pca<-function(Y,L=2,center=TRUE,scale=TRUE,rmzero=TRUE,ret_obj=FALSE){
  Y<-as.matrix(Y)
  if(rmzero==TRUE) Y<-rm_zero_rowcol(Y)
  res<-prcomp(as.matrix(t(Y)),center=center,scale.=scale,rank.=L)
  factors<-as.data.frame(res$x)
  colnames(factors)<-paste0("dim",1:L)
  if(ret_obj){
    return(list(factors=factors,obj=res))
  } else{
    return(factors)
  }
}

tsne<-function(Y,L=2,center=TRUE,scale=TRUE,rmzero=TRUE,method="Rtsne",...){
  if(rmzero==TRUE){
    Yt<-t(rm_zero_rowcol(Y))
  } else {
    Yt<-t(Y)
  }
  if(center || scale){
    Yt<-scale(Yt,center=center,scale=scale)
  }
  Yt<-as.matrix(Yt)
  if(method=="Rtsne"){
    fit<-Rtsne::Rtsne(Yt,dims=L,...)$Y
  } else if(method=="tsne"){
    fit<-tsne::tsne(Yt,k=L,...)
  }
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}

zinbwave<-function(Y,L=2,parallel=FALSE){
  #Y is unnormalized counts not log transformed
  if(parallel){
    bp<-BiocParallel::bpparam()
  } else {
    bp<-BiocParallel::SerialParam()
  }
  suppressWarnings(fit<-zinbwave::zinbFit(as.matrix(Y), K=L, BPPARAM=bp))
  factors<-as.data.frame(zinbwave::getW(fit))
  colnames(factors)<-paste0("dim",1:L)
  rownames(factors)<-colnames(Y)
  factors
}