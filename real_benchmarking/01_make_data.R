# Create all filtered datasets for comparisons
library(DuoClustering2018)
source("./util/functions_genefilter.R")

fp<-file.path
mkdir<-function(d){
  if(!dir.exists(d)) dir.create(d,recursive=TRUE)
}

bp<-"./real_benchmarking/data"
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
G<-1500

for(sn in sce_names){
  pth<-fp(bp,sn)
  mkdir(pth)
  sce<-do.call(paste0("sce_full_",sn),list())
  #sce<-sce_full_Zhengmix4eq()
  saveRDS(filterExpr(sce,G),file=fp(pth,"expr.rds"))
  saveRDS(filterHVG(sce,G),file=fp(pth,"hvg.rds"))
  saveRDS(filterDev(sce,G),file=fp(pth,"dev.rds"))
}
