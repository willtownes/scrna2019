#Miscellaneous functions used by other scripts

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

poisson_deviance_residuals<-function(x,xhat){
  #x,xhat assumed to be same dimension
  #sz<-exp(offsets)
  #xhat<-mu*sz
  term1<-x*log(x/xhat)
  term1[is.nan(term1)]<-0 #0*log(0)=0
  s2<-2*(term1-(x-xhat))
  sign(x-xhat)*sqrt(abs(s2))
}

multinomial_deviance_residuals<-function(X,p,n){
  #X a matrix, n is vector of length ncol(X)
  #if p is matrix, must have same dims as X
  #if p is vector, its length must match nrow(X)
  if(length(p)==nrow(X)){
    mu<-outer(p,n)
  } else if(!is.null(dim(p)) && dim(p)==dim(X)){
    mu<-t(t(p)*n)
  } else { stop("dimensions of p and X must match!") }
  sign(X-mu)*sqrt(-2*X*log(p))
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
  term1[is.nan(term1)]<-0 #0*log(0)=0
  nx<- t(n-t(X))
  term2<-nx*log(nx/outer(1-p,n))
  term2[is.nan(term2)]<-0
  sign(X-mu)*sqrt(2*(term1+term2))
}

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
        return(multinomial_deviance_residuals(m,phat,sz))
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
