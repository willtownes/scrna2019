# time Rscript ./real_benchmarking/03_clustering.R
source("./util/clustering.R")

fp<-file.path
mkdir<-function(d){
  if(!dir.exists(d)) dir.create(d,recursive=TRUE)
}

#### Script Part ####

ep<-"./real_benchmarking/embeddings"
cp<-"./real_benchmarking/clusterings"
dats<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
fils<-c("dev","expr","hvg")
#mths<-c("pca_log","pca_rp","pca_rd","glmpca")
#is_pca<-grepl("^pca_",mths)
#dims<-c(10,30)
resolutions<-c(0.05,0.1,0.2,0.5,0.8,1.0,1.2,1.5,2.0)
ks<-c(2,3,4,5,6,7,8,9,10)

master<-function(dat,fil,mths=c("pca_log","pca_rp","pca_rd","glmpca","zinbwave"),dims=c(10,30),genes=c(1500,300,60),verbose=TRUE){
  #ipth<-fp(ep,dat,fil)
  #opth<-fp(cp,dat,fil)
  #mkdir(opth)
  is_pca<-grepl("^pca_",mths)
  for(G in genes){
    ipth<-fp(ep,dat,fil,paste0("G",G))
    opth<-fp(cp,dat,fil,paste0("G",G))
    mkdir(opth)
    for(m in mths[is_pca]){
      fname_out<-fp(opth,paste0(m,".txt"))
      if(!file.exists(fname_out)){
        if(verbose) print(fname_out)
        fname_in<-fp(ipth,paste0(m,".txt"))
        embed<-read.table(fname_in,row.names=1)
        km<-kmeans_cluster(embed,dims=dims,k=ks)
        mc<-mclust_cluster(embed,dims=dims[dims<=20],k=ks)
        sc<-seurat_cluster(embed,dims=dims,res=resolutions)
        ans<-rbind(km,mc,sc)
        write.table(ans,file=fname_out,quote=FALSE,row.names=FALSE)
      }
    }
    for(m in mths[!is_pca]){
      fname_out<-fp(opth,paste0(m,".txt"))
      if(!file.exists(fname_out)){
        if(verbose) print(fname_out)
        func<-function(L){
          fname_in<-fp(ipth,paste0(m,"_L",L,".txt"))
          embed<-read.table(fname_in,row.names=1)
          km<-kmeans_cluster(embed,dims=L,k=ks)
          sc<-seurat_cluster(embed,dims=L,res=resolutions)
          ans<-rbind(km,sc)
          if(L<=20){
            ans<-rbind(ans,mclust_cluster(embed,dims=L,k=ks))
          }
          ans
        }
        ans<-do.call(rbind,lapply(dims,func))
        write.table(ans,file=fname_out,quote=FALSE,row.names=FALSE)
      }
    }
  }
}

atlas<-expand.grid(dats=dats,fils=fils)
parallel::mcmapply(master,atlas$dats,atlas$fils,mc.cores=4)
