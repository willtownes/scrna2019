#miscellaneous functions useful for text parsing, etc

split2mat<-function(x,splt,cnames=NULL){
  #takes a character vector x and string splits on 'splt" character
  #builds a matrix with nrow=length(x) and ncol=number of substrings
  #only works if each element of x has the splt character appearing same number of times
  x<-strsplit(as.character(x),splt)
  x<-t(matrix(unlist(x), nrow=length(x[[1]])))
  if(!is.null(cnames)) colnames(x)<-cnames
  x
}

rm_dupcols<-function(df){
  #takes a dataframe and removes all columns that have the same value in all rows
  fn<-function(x){all(is.na(df[,x])) || nlevels(as.factor(df[,x]))==1}
  dups<-vapply(1:ncol(df),fn,TRUE)
  df[,!dups]
}