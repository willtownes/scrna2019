# time Rscript ./real_benchmarking/02_dimreduce.R
# compute cell embeddings
library(SingleCellExperiment)
source("./util/functions.R") #null_residuals function
source("./algs/existing.R") #pca and zinbwave functions
source("./algs/glmpca.R")

fp<-file.path
mkdir<-function(d){
  if(!dir.exists(d)) dir.create(d,recursive=TRUE)
}

dimreduce<-function(Y,dims,mth=c("pca_log","pca_rp","pca_rd","glmpca","zinbwave")){
  #Y a count matrix with cells in the columns, genes in the rows
  mth<-match.arg(mth)
  Y<-rm_zero_rowcol(Y)
  sz<-colSums(Y)
  if(mth=="pca_log"){
    Y2<-log2(1+1e6*t(t(Y)/sz))
  } else if(mth=="pca_rp") {
    Y2<-null_residuals(Y,type="pearson",mod="binomial")
  } else if(mth=="pca_rd") {
    Y2<-null_residuals(Y,type="deviance",mod="binomial")
  } else if(mth=="glmpca") {
    return(glmpca(Y,dims)$factors)
  } else if(mth=="zinbwave") {
    return(zinbwave(Y,dims))
    #stop("not implemented yet")
  }
  pca(Y2,dims,center=TRUE,scale=TRUE)
}

#### Script Part ####

dats<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
fils<-c("dev","expr","hvg")
ep<-"./real_benchmarking/embeddings"
dp<-"./real_benchmarking/data"

master<-function(dat,fil,mths=c("pca_log","pca_rp","pca_rd","glmpca","zinbwave"),dims=c(10,30),genes=c(60,300,1500),verbose=TRUE){
  pth<-fp(ep,dat,fil)
  #mkdir(pth)
  sce<-readRDS(fp(dp,dat,paste0(fil,".rds")))
  Y0<-assay(sce,"counts")
  is_pca<-grepl("^pca_",mths)
  for(G in genes){
    Y<-Y0[1:G,]
    gpth<-fp(pth,paste0("G",G))
    mkdir(gpth)
    for(m in mths[is_pca]){
      fname<-fp(gpth,paste0(m,".txt"))
      if(!file.exists(fname)){
        if(verbose) print(fname)
        res<-dimreduce(Y,max(dims),m)
        write.table(res,file=fname,quote=FALSE,col.names=FALSE)
      }
    }
    for(m in mths[!is_pca]){
      for(L in dims){
        fname<-fp(gpth,paste0(m,"_L",L,".txt"))
        if(!file.exists(fname)){
          if(verbose) print(fname)
          res<-dimreduce(Y,L,m)
          write.table(res,file=fname,quote=FALSE,col.names=FALSE)
        }
      }
    }
  }
}

atlas<-expand.grid(dats=dats,fils=fils)
parallel::mcmapply(master,atlas$dats,atlas$fils,mc.cores=4)
