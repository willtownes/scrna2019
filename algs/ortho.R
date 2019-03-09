#functions for orthogonalizing factors and loadings

norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}

ortho<-function(U,V,A,B,W=rep(1,nrow(V)),X=rep(1,nrow(U)),ret=c("m","df")){
  #U is NxL matrix of cell factors, V is GxL matrix of loadings onto genes
  #W is GxJ matrix of gene specific covariates
  #A is NxJ matrix of coefficients of W
  #X is NxK matrix of cell specific covariates
  #B is GxK matrix of coefficients of X
  #ret= return either data frame or matrix format
  #imputed expression: E[Y] ~= R = VU'+WA'+BX'
  ret<-match.arg(ret)
  L<-ncol(U)
  W<-as.matrix(W,nrow=nrow(V))
  X<-as.matrix(X,nrow=nrow(U))
  if(!is.null(B) && B!=0){
    reg<-lm(U~X-1)
    factors<-residuals(reg)
    B<-B+tcrossprod(V,coef(reg))
  } else {
    factors<-U
  }
  if(!is.null(A) && A!=0){
    reg<-lm(V~W-1)
    loadings<-residuals(reg)
    A<-A+tcrossprod(factors,coef(reg))
  } else {
    loadings<-V
  }
  svdres<-svd(loadings)
  loadings<-svdres$u
  factors<-t(t(factors%*%svdres$v)*svdres$d)
  o<-order(colNorms(factors),decreasing=TRUE)
  factors<-factors[,o]
  loadings<-loadings[,o]
  colnames(loadings)<-colnames(factors)<-paste0("dim",1:L)
  if(ret=="df"){
    loadings<-as.data.frame(loadings)
    factors<-as.data.frame(factors)
    if(!is.null(A)) A<-as.data.frame(A)
    if(!is.null(B)) B<-as.data.frame(B)
  }
  return(list(factors=factors,loadings=loadings,A=A,B=B))
}
