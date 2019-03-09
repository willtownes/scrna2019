# GLM-PCA
source("./algs/ortho.R")

mat_binom_dev<-function(X,P,n){
  #binomial deviance for two matrices
  #X,P are GxN matrices
  #n is vector of length N (same as cols of X,P)
  #nz<-X>0
  X<-t(X); P<-t(P)
  term1<-sum(X*log(X/(n*P)), na.rm=TRUE)
  #nn<-x<n
  nx<-n-X
  term2<-sum(nx*log(nx/(n*(1-P))), na.rm=TRUE)
  2*(term1+term2)
}

glmpca<-function(Y,L,fam=c("poi","nb","mult","bern"),ctl=list(maxIter=100,eps=1e-4),penalty=1,verbose=FALSE,init=list(factors=NULL,loadings=NULL),nb_theta=1){
  #Y is data with genes=rows, samples=cols
  #L is number of desired latent dimensions
  #family same as in glm() or glm.nb()
  #for negative binomial family, dispersion param theta must be provided
  #note, colSum(Y) offsets included automatically for count data,
  #and as weights for family=binomial()
  #if the data matrix is binary, binomial means bernoulli (no weights)
  #if data matrix is counts, we assume each column is a multinomial sample,
  #and the binomial likelihood is used as an approximation to multinomial,
  #the total counts in a column is the "n" in the binomial.
  #Link functions: log for poi, nb; logit for mult,bern
  N<-ncol(Y); G<-nrow(Y)
  fam<-match.arg(fam)
  #sanity check inputs
  stopifnot(min(Y)>=0)
  if(fam=="bern"){ stopifnot(max(Y)<=1) }
  #create GLM family object
  if(fam=="poi"){
    family<-poisson()
  } else if(fam=="nb"){
    family<-MASS::negative.binomial(theta=nb_theta)
  } else {
    family<-binomial()
  }
  #variance function, determined by GLM family
  vfunc<-family$variance
  #mu as a function of eta
  ilfunc<-family$linkinv
  #derivative of inverse link function, dmu/deta
  dfunc<-family$mu.eta 
  if(fam %in% c("poi","nb")){
    sz<-colMeans(Y) #size factors
    offsets<-family$linkfun(sz)
    etafunc<-function(U,V){ t(offsets+tcrossprod(U,V)) } #linear predictor
    v1<-family$linkfun(rowSums(Y)/sum(sz)) 
    if(fam=="poi"){
      infograd<-function(eta){
        M<-ilfunc(eta)
        list(grad=(Y-M),info=M)
      }
    } else if(fam=="nb"){
      infograd<-function(eta){
        M<-ilfunc(eta)
        W<-1/vfunc(M)
        list(grad=(Y-M)*W*M, info=W*M^2)
      }
    }
  } else if(fam == "mult"){
    sz<-colSums(Y) #n in the multinomial
    etafunc<-function(U,V){ tcrossprod(V,U) }
    v1<-family$linkfun(rowSums(Y)/sum(sz))
    infograd<-function(eta){
      #eta is GxN matrix
      Pt<-t(ilfunc(eta)) #ilfunc=expit, Pt very small probabilities
      list(grad=Y-t(sz*Pt), info=t(sz*vfunc(Pt)))
    }
  } else { #no offsets in bernoulli data
    etafunc<-function(U,V){ tcrossprod(V,U) }
    v1<-family$linkfun(rowMeans(Y))
    if(fam=="bern"){
      infograd<-function(eta){
        P<-ilfunc(eta)
        list(grad=Y-P, info=vfunc(P))
      }
    } else { #this is not actually used but keeping for future reference
      #this is most generic formula for GLM but computationally slow
      stop("invalid fam")
      infograd<-function(eta){
        M<-ilfunc(eta)
        W<-1/vfunc(M)
        D<-dfunc(eta)
        list(grad=(Y-M)*W*D, info=W*D^2)
      }
    }
  }
  #create deviance function
  if(fam=="mult"){
    dev_func<-function(U,V){
      #note the global variables Y and sz
      mat_binom_dev(Y,ilfunc(etafunc(U,V)),sz)
    }
  } else {
    dev_func<-function(U,V){
      #global variable: Y
      #Note dev.resids function gives the square of actual residuals
      #so the sum of these is the deviance
      sum(family$dev.resids(Y,ilfunc(etafunc(U,V)),1)) 
    }
  }
  #initialize U,V, with row-specific intercept term
  U<-cbind(1, matrix(1e-5/L*rnorm(N*L),nrow=N))
  if(!is.null(init$factors)){
    print("initialize factors")
    L0<-min(L,ncol(init$factors))
    U[,2:(L0+1)]<-as.matrix(init$factors[,1:L0])
  }
  #v1 = naive MLE for gene intercept only
  V<-cbind(v1,matrix(1e-5/L*rnorm(G*L),nrow=G))
  if(!is.null(init$loadings)){
    print("initialize loadings")
    L0<-min(L,ncol(init$loadings))
    V[,2:(L0+1)]<-as.matrix(init$loadings[,1:L0])
  }
  #M<-ilfunc(etafunc(U,V))
  vid<-2:(L+1)
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){
    #rmse[t]<-sd(Y-ilfunc(etafunc(U,V)))
    dev[t]<-dev_func(U,V)
    if(t>5 && abs(dev[t]-dev[t-1])/(0.1+abs(dev[t-1]))<ctl$eps){
      break
    }
    if(verbose){ 
      dev_format<-format(dev[t],scientific=TRUE,digits=4)
      print(paste0("Iteration: ",t," | deviance=",dev_format)) 
    }
    for(l in 1:(L+1)){
      ig<- infograd(etafunc(U,V))
      #no penalty on intercept terms
      grads<- (ig$grad)%*%U[,l] - penalty*V[,l]*(l>1) 
      infos<- (ig$info) %*% U[,l]^2 + penalty*(l>1)
      V[,l]<-V[,l]+grads/infos
    }
    for(l in vid){
      ig<- infograd(etafunc(U,V))
      grads<- crossprod(ig$grad, V[,l]) - penalty*U[,l]
      infos<- crossprod(ig$info, V[,l]^2) + penalty
      U[,l]<-U[,l]+grads/infos
    }
  }
  res<-ortho(U[,vid],V[,vid],0,V[,1],ret="df")
  rownames(res$factors)<-colnames(Y)
  rownames(res$loadings)<-rownames(Y)
  res$dev=dev[1:t]; res$fam<-fam
  res
}
