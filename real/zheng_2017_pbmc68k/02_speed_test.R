# Rscript ./real/zheng_2017_pbmc68k/02_speed_test.R
# Dimension Reduction Speed Tests
# Assumes already ran 01_data_loading.Rmd
# Uses the Zheng 2017 PBMC 68K dataset
library(SingleCellExperiment)
library(glmpca)
source("./util/functions.R") #null_residuals function
source("./algs/existing.R") #pca and zinbwave functions
fp<-file.path

funcs<-list()
funcs$pca_log<-function(Y,dims){
  sz<-colSums(Y)
  Y2<-log2(1+1e6*t(t(Y)/sz))
  pca(Y2,dims,center=TRUE,scale=TRUE)
}
funcs$pca_rp<-function(Y,dims){
  Y2<-null_residuals(Y,type="pearson",mod="binomial")
  pca(Y2,dims,center=TRUE,scale=TRUE)
}
funcs$pca_rd<-function(Y,dims){
  Y2<-null_residuals(Y,type="deviance",mod="binomial")
  pca(Y2,dims,center=TRUE,scale=TRUE)
}
funcs$glmpca<-function(Y,dims){
  glmpca(Y,dims)$factors
}
funcs$zinbwave<-function(Y,dims){
  zinbwave(Y,dims)
}

pth<-"./real/zheng_2017_pbmc68k"
G<-600; L<-2
sce<-readRDS(fp(pth,"data/01_sce_all_genes_all_cells.rds"))
odir<-fp(pth,"results")
if(!dir.exists(fp(odir,"factors"))){
  dir.create(fp(odir,"factors"),recursive=TRUE)
}

Y0<-rm_zero_rowcol(counts(sce[1:G,]))

Nvec<-c(680,6800,68000)
spf<-fp(odir,"speed.txt")
if(!file.exists(spf)){
  res<-expand.grid(mths=names(funcs),N=Nvec)
  res$elapsed<-res$system<-res$user<-NA
} else {
  res<-read.table(spf,header=TRUE)
}
for(N in Nvec){
  Y<-as.matrix(Y0[,1:N])
  for(mth in names(funcs)){
    id<- which(res$mths==mth & res$N==N)
    if(any(is.na(res[id,3:5]))){
      print(paste(mth,N))
      a<-list(Y,L)
      f<-funcs[[mth]]
      res[id,3:5]<-system.time(factors<-do.call(f,a))[1:3]
      fname<-paste0(mth,"_N",N,"_G",G,"_L",L,".txt")
      ofile<-fp(odir,"factors",fname)
      write.table(factors,file=ofile,quote=FALSE,col.names=FALSE)
      write.table(res,file=fp(odir,"speed.txt"),quote=FALSE,row.names=FALSE)
    }
  }
}
