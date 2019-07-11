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
ep<-"./real_benchmarking/embeddings"
cp<-"./real_benchmarking/clusterings"
bp<-"./real_benchmarking/results"
mkdir(bp)
dats<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
fils<-c("dev","expr","hvg")
genes=c(1500,300,60)
mths<-c("pca_log","pca_rd","pca_rp","glmpca","zinbwave")
dims<-c(10,30)
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
      ipth<-fp(ep,d,f,paste0("G",G))
      for(m in mths){
        is_pca<-grepl("^pca_",m)
        fname<-fp(cp,d,f,paste0("G",G),paste0(m,".txt"))
        cl<-read.table(fname,header=TRUE)
        cl$tru<-tru[cl$cell,"phenoid"]
        if(is_pca){ embed0<-read.table(fp(ipth,paste0(m,".txt")),row.names=1) }
        for(L in dims){
          key<-paste(d,f,paste0("G",G),m,paste0("L",L),sep="-")
          print(key)
          if(is_pca){ 
            embed<-embed0[,1:L] 
          } else {
            embed<-read.table(fp(ipth,paste0(m,"_L",L,".txt")),row.names=1)
          }
          DM<-dist(embed)
          silfunc<-function(clust){
            #case where only one cluster, set silhouette to zero
            if(length(unique(clust))==1){ return(0) }
            sil<-cluster::silhouette(as.numeric(clust),DM)
            tryCatch(mean(sil[,3]),error=function(e){NA})
          }
          cl2<-cl %>% subset(dims==L) %>% group_by(k,resolution,method) %>% summarise(ari=ari(cluster,tru),silhouette=silfunc(cluster))#,jac=jac(cluster,tru))
          #cl2$dat<-d; cl2$genefilter<-f; cl2$num_genes<-G; cl2$dimreduce<-m; cl2$dims==L
          ans[[key]]<-cl2 %>% mutate(dat=d,genefilter=f,num_genes=G,dimreduce=m,dims=L)
        }
      }
    }
  }
}
ans<-do.call(rbind,ans)
ans2<-ans[,c(6:ncol(ans),c(3,1,2,4,5))]
ans2<-ans2 %>% arrange(dat,genefilter,num_genes,dimreduce,dims,method,k,resolution)
write.table(ans2,fp(bp,"cluster_accuracy.txt"),quote=FALSE,row.names=FALSE)
