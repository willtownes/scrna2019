library(dplyr)
library(SingleCellExperiment)
ari<-mclust::adjustedRandIndex
#jac<-function(x,y){clusteval::cluster_similarity(x,y,similarity="jaccard")}
fp<-file.path
mkdir<-function(d){
  if(!dir.exists(d)) dir.create(d,recursive=TRUE)
}

#### Script Part ####

dp<-"./real_benchmarking/data"
cp<-"./real_benchmarking/clusterings"
bp<-"./real_benchmarking/results"
mkdir(bp)
dats<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
fils<-c("dev","expr","hvg")
genes=c(1500,300,60)
mths<-c("pca_log","pca_rd","pca_rp","glmpca","zinbwave")
#dims<-c(10,20)
#resolutions<-c(0.1,0.3,1.0,1.1,1.3,1.4,1.6,2.0)
#ks<-seq(from=2,to=14,by=2)

ans<-list()
for(d in dats){
  for(f in fils){
    sce<-readRDS(fp(dp,d,paste0(f,".rds")))
    tru<-colData(sce)["phenoid"] #retain rowids
    #opth<-fp(bp,d,f)
    #mkdir(opth)
    for(G in genes){
      for(m in mths){
        key<-paste(d,f,paste0("G",G),m,sep="-")
        print(key)
        fname<-fp(cp,d,f,paste0("G",G),paste0(m,".txt"))
        cl<-read.table(fname,header=TRUE)
        cl$tru<-tru[cl$cell,"phenoid"]
        cl2<-cl %>% group_by(k,resolution,dims,method) %>% summarise(ari=ari(cluster,tru))#,jac=jac(cluster,tru))
        cl2$dat<-d; cl2$genefilter<-f; cl2$num_genes<-G; cl2$dimreduce<-m
        ans[[key]]<-cl2
      }
    }
  }
}
ans<-do.call(rbind,ans)
write.table(ans,fp(bp,"cluster_accuracy.txt"),quote=FALSE,row.names=FALSE)
