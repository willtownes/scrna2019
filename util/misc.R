#miscellaneous functions

parcolapply<-function(m,f,transpose=TRUE,chunksize=50,cores=max(parallel::detectCores()-1,1),...){
  #same thing as apply(m,2,f) but using parallel processing for speed
  #also if m is sparse matrix this avoids casting the whole m to dense
  #m, a matrix for which we want to apply a function to each column
  #f, the function to be applied to each column of m
  #if transpose=TRUE, we assume f returns a data frame or matrix
  #whose rows match cols of m
  #if transpose=FALSE, we assume f returns a matrix whose columns match cols of m
  #chunks, number of sub-matrices to split m into
  #cores, number of parallel processors
  #... additional arg for f
  # https://stackoverflow.com/questions/45198194/partition-matrix-into-n-equally-sized-chunks-with-r
  chunks<-ceiling(ncol(m)/chunksize)
  cc<-seq_len(ncol(m))
  cf<-cut(cc,pretty(cc,chunks))
  ml<-lapply(split(cc,cf),function(ii){m[,ii]})
  res<-parallel::mclapply(ml,f,...,mc.cores=cores)
  names(res)<-NULL
  if(transpose){
    return(do.call(rbind,res))
  } else {
    return(do.call(cbind,res))
  }
}