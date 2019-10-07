#Miscellaneous functions used by other scripts

##### Data Loading functions #####

get_10x_readcounts<-function(counts_folder,mol_info_h5){
  #provide two file paths:
  #counts_folder, containing "barcodes.tsv","genes.tsv","matrix.mtx"
  #mol_info_h5, an hdf5 file with the 10x molecule information file
  #output: a SingleCellExperiment with two sparse Matrix assays:
  #"counts" (UMI counts), and "read_counts"
  sce0<-DropletUtils::read10xCounts(counts_folder)
  colnames(sce0)<-colData(sce0)$Barcode
  m<-assay(sce0,"counts")
  gg<-Matrix::rowSums(m)>0 #remove genes that are all zero
  sce<-sce0[gg,]
  m<-m[gg,]
  cm<-colData(sce0)
  bc<-substr(cm$Barcode,1,nchar(cm$Barcode)-2) #assumes barcode ends in "-1"
  cm$bc_enc<-DropletUtils::encodeSequences(bc)
  mi0<-as.data.frame(rhdf5::h5dump(mol_info_h5))
  if(!all(cm$bc_enc %in% mi0$barcode)){
    stop("Some count matrix barcodes not found in molecule information file!")
  }
  gidx<-read.table(file.path(counts_folder,"genes.tsv"),stringsAsFactors=FALSE)[,1]
  mi<-subset(mi0, barcode %in% cm$bc_enc & gene<length(gidx))
  mi$gene_symbol<-gidx[mi$gene+1]
  cm2<-cm[,c("bc_enc","Barcode")]
  colnames(cm2)<-c("barcode","barcode_str")
  mi<-merge(mi,cm2,by="barcode")
  rc<-DropletUtils::makeCountMatrix(mi$gene_symbol,mi$barcode_str,value=mi$reads)
  rc<-rc[rownames(m),colnames(m)]
  #umi<-DropletUtils::makeCountMatrix(mi$gene_symbol,mi$barcode_str)
  #umi<-umi[rownames(m),colnames(m)]
  #umi counts from molecule info file consistent with sce
  #all(m==umi)
  if(!all((rc>0) == (m>0) )){
    stop("zero pattern inconsistent between counts matrix and molecule info")
  }
  if(!all(rc>=m)){
    stop("Read counts not all >= UMI counts")
  }
  assay(sce,"read_counts")<-rc
  sce
}

##### ID conversion functions #####

ensembl2symbol<-function(genes,sp=c("hsapiens","mmusculus")){
  #genes is a character vector of ensembl IDs
  #sp=the species
  #returns a data frame with the ensembl ID mapped to symbols
  #if the gene didn't have a corresponding symbol, it is removed from the data frame.
  sp<-match.arg(sp)
  symb<-list(hsapiens="hgnc",mmusculus="mgi")
  bm<-biomaRt::useMart("ensembl")
  mart<-biomaRt::useDataset(paste0(sp,"_gene_ensembl"),bm)
  bml<-biomaRt::getBM(filters="ensembl_gene_id",attributes=c("ensembl_gene_id",paste0(symb[[sp]],"_symbol")),values=genes,mart=mart)
  #make sure order matches the input genes order
  #if no matching symbol was found, use NA
  genes<-as.data.frame(genes,stringsAsFactors=FALSE)
  colnames(genes)<-"ensembl_gene_id"
  plyr::join(genes,bml,by="ensembl_gene_id",type="left",match="first")
}

##### Downsampling functions #####

Down_Sample_Matrix<-function(expr_mat,min_lib_size=NULL){
  #adapted from https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalisations
  min_sz<-min(colSums(expr_mat))
  if(is.null(min_lib_size)){
    min_lib_size<-min_sz
  } else {
    stopifnot(min_lib_size<=min_sz)
  }
  down_sample<-function(x){
    prob <- min_lib_size/sum(x)
    unlist(lapply(x,function(y){rbinom(1, y, prob)}))
  }
  apply(expr_mat, 2, down_sample)
}

##### Deviance functions #####

poisson_deviance<-function(x,mu,sz){
  #assumes log link and size factor sz on the same scale as x (not logged)
  #stopifnot(all(x>=0 & sz>0))
  2*sum(x*log(x/(sz*mu)),na.rm=TRUE)-2*sum(x-sz*mu)
}

multinomial_deviance<-function(x,p){
  -2*sum(x*log(p))
}

binomial_deviance<-function(x,p,n){
  term1<-sum(x*log(x/(n*p)), na.rm=TRUE)
  nx<-n-x
  term2<-sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
  2*(term1+term2)
}

gof_func<-function(x,sz,mod=c("binomial","multinomial","poisson","geometric")){
  #Let n=colSums(original matrix where x is a row)
  #if binomial, assumes sz=n, required! So sz>0 for whole vector
  #if poisson, assumes sz=n/geometric_mean(n), so again all of sz>0
  #if geometric, assumes sz=log(n/geometric_mean(n)) which helps numerical stability. Here sz can be <>0
  #note sum(x)/sum(sz) is the (scalar) MLE for "mu" in Poisson and "p" in Binomial
  mod<-match.arg(mod)
  fit<-list(deviance=0,df.residual=length(x)-1,converged=TRUE)
  if(mod=="multinomial"){
    fit$deviance<-multinomial_deviance(x,sum(x)/sum(sz))
  } else if(mod=="binomial"){
    fit$deviance<-binomial_deviance(x,sum(x)/sum(sz),sz)
  } else if(mod=="poisson"){
    fit$deviance<-poisson_deviance(x,sum(x)/sum(sz),sz)
  } else if(mod=="geometric"){
    if(any(x>0)) {
      fit<-glm(x~offset(sz),family=MASS::negative.binomial(theta=1))
    }
  } else { stop("invalid model") }
  if(fit$converged){
    dev<-fit$deviance
    df<-fit$df.residual #length(x)-1
    pval<-pchisq(dev,df,lower.tail=FALSE)
    res<-c(dev,pval)
  } else {
    res<-rep(NA,2)
  }
  names(res)<-c("deviance","pval")
  res
}

compute_size_factors<-function(m,mod=c("binomial","multinomial","poisson","geometric")){
  #given matrix m with samples in the columns
  #compute size factors suitable for the discrete model in 'mod'
  mod<-match.arg(mod)
  sz<-Matrix::colSums(m) #base case, multinomial or binomial
  if(mod %in% c("multinomial","binomial")){ return(sz) }
  sz<-log(sz)
  sz<-sz - mean(sz) #make geometric mean of sz be 1 for poisson, geometric
  if(mod=="poisson"){ return(exp(sz)) }
  sz #geometric, use log scale size factors
}

compute_gene_info<-function(m,gmeta=NULL,mod=c("binomial","multinomial","poisson","geometric")){
  #m a data matrix with genes=rows
  #gmeta a pre-existing data frame with gene-level metadata
  mod<-match.arg(mod)
  if(!is.null(gmeta)){ stopifnot(nrow(m)==nrow(gmeta)) }
  gnz<-Matrix::rowSums(m>0)
  sz<-compute_size_factors(m,mod)
  gof<-function(g){ gof_func(m[g,],sz,mod) }
  gof<-as.data.frame(t(vapply(1:nrow(m),gof,FUN.VALUE=rep(0.0,2))))
  #colnames(gof)<-c("deviance","pval")
  gmu<-Matrix::rowMeans(m)
  gvar<-apply(m,1,var)
  gfano<-ifelse(gvar>0 & gmu>0, gvar/gmu, 0)
  res<-cbind(nzsum=gnz,fano=gfano,gof)
  res$pval_fdr<-p.adjust(res$pval,"BH")
  if(is.null(gmeta)){ return(res) } else { return(cbind(gmeta,res)) }
}

##### Null Residuals functions #####

poisson_deviance_residuals<-function(x,xhat){
  #x,xhat assumed to be same dimension
  #sz<-exp(offsets)
  #xhat<-mu*sz
  term1<-x*log(x/xhat)
  term1[is.nan(term1)]<-0 #0*log(0)=0
  s2<-2*(term1-(x-xhat))
  sign(x-xhat)*sqrt(abs(s2))
}

binomial_deviance_residuals<-function(X,p,n){
  #X a matrix, n is vector of length ncol(X)
  #if p is matrix, must have same dims as X
  #if p is vector, its length must match nrow(X)
  if(length(p)==nrow(X)){
    mu<-outer(p,n)
  } else if(!is.null(dim(p)) && dim(p)==dim(X)){
    mu<-t(t(p)*n)
  } else { stop("dimensions of p and X must match!") }
  term1<-X*log(X/mu)
  # term1[is.nan(term1)]<-0 #0*log(0)=0
  term1[is.na(term1)]<-0 #0*log(0)=0
  nx<- t(n-t(X))
  term2<-nx*log(nx/outer(1-p,n))
  term2<-nx*log(nx/outer(1-p,n))
  # term2[is.nan(term2)]<-0
  sign(X-mu)*sqrt(2*(term1+term2))
}

# multinomial_deviance_residuals<-function(X,p,n){
#   #not clear if this actually makes sense. Don't use this function!
#   #X a matrix, n is vector of length ncol(X)
#   #if p is matrix, must have same dims as X
#   #if p is vector, its length must match nrow(X)
#   if(length(p)==nrow(X)){
#     mu<-outer(p,n)
#   } else if(!is.null(dim(p)) && dim(p)==dim(X)){
#     mu<-t(t(p)*n)
#   } else { stop("dimensions of p and X must match!") }
#   sign(X-mu)*sqrt(-2*X*log(p))
# }

null_residuals<-function(m,mod=c("binomial","multinomial","poisson","geometric"),type=c("deviance","pearson")){
  mod<-match.arg(mod)
  type<-match.arg(type)
  sz<-compute_size_factors(m,mod)
  if(mod %in% c("multinomial","binomial")) {
    phat<-Matrix::rowSums(m)/sum(sz)
    if(type=="pearson"){ #pearson resids same for multinomial as binomial
      mhat<-outer(phat,sz)
      return((m-mhat)/sqrt(mhat*(1-phat)))
    } else { #deviance residuals
      if(mod=="multinomial"){
        stop("multinomial deviance residuals not implemented yet")
        #return(multinomial_deviance_residuals(m,phat,sz))
      } else { #binomial
        return(binomial_deviance_residuals(m,phat,sz))
      }
    }
  } else if(mod=="poisson"){
    mhat<-outer(Matrix::rowSums(m)/sum(sz), sz) #first argument is "lambda hat" (MLE)
    if(type=="deviance"){ 
      return(poisson_deviance_residuals(m,mhat))
    } else { #pearson residuals
      return((m-mhat)/sqrt(mhat))
    }
  } else { #geometric
    gfunc<-function(g){
      fit<-glm(m[g,]~offset(sz),family=MASS::negative.binomial(theta=1))
      residuals(fit,type=type)
    }
    return(t(vapply(1:nrow(m),gfunc,rep(0.0,ncol(m)))))
  }
}

##### BIC functions #####

#throughout, m is the expression matrix with features=rows, samples=cols

mult_bic<-function(m){
  #multinomial model
  n<-colSums(m)
  p<-rowSums(m)/sum(n)
  ll<-sum(apply(m,2,dmultinom,prob=p,log=TRUE))
  df<-length(p)-1
  -2*ll+df*log(prod(dim(m)))
}

dmn_bic<-function(m){
  fit<-DirichletMultinomial::dmn(t(m),1)
  alpha<-drop(fit@fit$Estimate)
  ll<-sum(extraDistr::ddirmnom(t(m), colSums(m), alpha, log = TRUE))
  -2*ll+nrow(m)*log(prod(dim(m)))
}

poi_fit<-function(m,X=NULL,sz=NULL,maxit=100){
  #poisson
  if(is.null(sz)){ sz<-log(colMeans(m)) }
  if(is.null(X)){ #no covariates => closed form solution
    sz<-exp(sz)
    lam<-rowSums(m)/sum(sz)
    mu<-outer(lam,sz)
    ll<-matrix(dpois(m,mu,log=TRUE),nrow=nrow(m))
    return(data.frame(ll=rowSums(ll),converged=TRUE))
  } else { #covariates included
    stopifnot(all(X[,1]==1))
    k<-ncol(X)
    fam<-poisson()
    #default maxit is 25 for glm.fit, but some fits take longer to converge
    ctl<-list(maxit=maxit) 
    f<-function(y){
      fit<-glm.fit(X,y,offset=sz,control=ctl,family=fam)
      c(ll=k-fit$aic/2, converged=fit$converged)
    }
    res<-as.data.frame(t(apply(m,1,f)))
    res$converged<-as.logical(res$converged)
    return(res)
  }
}

poi_bic<-function(m,X=NULL,prefit=NULL,maxit=100){
  #poisson. prefit should be a data frame with column "ll" for log-likelihood
  k<-if(is.null(X)){ 1 } else { ncol(X) }
  df<-k*nrow(m)
  if(is.null(prefit)){
    prefit<-poi_fit(m,X,maxit=maxit)
  }
  ll<-sum(prefit$ll)
  #compute BIC: -2*loglik+df*log(n_obs)
  -2*ll+df*log(prod(dim(m)))
}

nb_fit<-function(m,X=NULL,sz=NULL){
  #neg binom
  if(is.null(sz)){ sz<-log(colMeans(m)) }
  if(is.null(X)){
    f<-function(y){
      mgcv::gam(y~1,family=mgcv::nb,offset=sz,method="ML")
    }
  } else {
    f<-function(y){
      mgcv::gam(y~X-1,family=mgcv::nb,offset=sz,method="ML")
    }
  }
  g<-function(y){
    fit<-f(y)
    th<-fit$family$getTheta(TRUE) #FALSE (default) gives log(theta)
    ll<-as.numeric(logLik(fit))
    c(theta=th,ll=ll,converged=fit$converged)
  }
  res<-as.data.frame(t(apply(m,1,g)))
  res$converged<-as.logical(res$converged)
  res
}

nb_bic<-function(m,X=NULL,prefit=NULL){
  if(is.null(prefit)){ 
    prefit<-nb_fit(m,X)
  } else {
    stopifnot(nrow(m)==nrow(prefit))
  }
  ll<-sum(prefit$ll)
  k<-if(is.null(X)){ 2 } else { ncol(X)+1 }
  -2*ll+(k*nrow(m))*log(prod(dim(m)))
}

zip_fit<-function(m,X=NULL,sz=NULL){
  #zero inflated poisson
  if(is.null(sz)){ sz<-log(colMeans(m)) }
  if(is.null(X)){
    f<-function(y){
      pscl::zeroinfl(y~offset(sz) | 1,dist="poisson")
    }
  } else {
    f<-function(y){
      pscl::zeroinfl(y~offset(sz)+X-1 | 1,dist="poisson")
    }
  }
  g<-function(y){
    fit<-f(y)
    c(pz=plogis(coef(fit)[2]),ll=logLik(fit),converged=fit$converged)
  }
  res<-as.data.frame(t(apply(m,1,g)))
  res$converged<-as.logical(res$converged)
  res
}

zip_bic<-function(m,X=NULL,prefit=NULL){
  if(is.null(prefit)){ 
    prefit<-zip_fit(m,X) 
  } else {
    stopifnot(nrow(m)==nrow(prefit))
  }
  ll<-sum(prefit$ll)
  k<-if(is.null(X)){ 2 } else { ncol(X)+1 }
  -2*ll+(k*nrow(m))*log(prod(dim(m)))
}

normal_bic<-function(m){
  sz<-colMeans(m)
  lam<-rowSums(m)/sum(sz)
  mu<-outer(lam,sz)
  r<-m-mu
  s<-apply(r,1,sd)
  z<-r/s
  ll<-sum(dnorm(z,mean=0,sd=1,log=TRUE))
  df<-2*nrow(m)
  -2*ll+df*log(prod(dim(m)))
}

lnorm_nz_fit<-function(m,X=NULL){
  sz<-colMeans(m)
  lra<-t(t(log(m))-log(sz))
  lra[is.infinite(lra)]<-NA
  if(is.null(X)){
    mu<-rowMeans(lra,na.rm=TRUE)
    r<-lra-mu #matrix with NAs
    s<-apply(r,1,sd,na.rm=TRUE)
    z<-r/s #matrix with NAs
    return(sum(dnorm(z[!is.na(z)],mean=0,sd=1,log=TRUE)))
  } else {
    f<-function(y){
      r<-residuals(lm(y~X-1,na.action="na.omit"))
      z<-r/sd(r)
      sum(dnorm(z,mean=0,sd=1,log=TRUE),na.rm=TRUE)
    }
    return(sum(apply(lra,1,f)))
  }
}

ziln_bic<-function(m,X=NULL){
  #zero inflated log-normal (hurdle model)
  Z<-m>0
  p<-rowMeans(Z)
  ll1<-sum(dbinom(Z,1,p,log=TRUE))
  ll2<-lnorm_nz_fit(m,X)
  #per gene df is 1 for pzero, 1 for sd, 1 for mean+number of covars
  if(is.null(X) || is.null(ncol(X))){ #X is a vector or NULL
    k<-3
  } else { #X is a matrix
    k<-2+ncol(X)
  }
  N<-prod(dim(m))
  -2*(ll1+ll2)+(k*nrow(m))*log(N+sum(Z))
}

bic_all<-function(m,X=NULL,liks=c("mult","dmn","poi","nb","zip","normal","ziln")){
  covar_liks<-c("poi","nb","zip","ziln")
  if(!is.null(X)){ liks<-intersect(liks,covar_liks) }
  f<-function(x){eval(call(paste0(x,"_bic"),m,X))}
  sapply(liks,f)
}