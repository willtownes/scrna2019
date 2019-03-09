#!/usr/bin/env Rscript
#example usage
#Rscript ../../algs/zinbwave_script.R -i data/MGH26/batch_separability/subsets_raw/1.txt -o data/MGH26/batch_separability/zinbwave_out/1.tsv

library(modules)
import_package("optparse",attach=TRUE)

opt_list<-list(
  make_option(c("-i","--infile"), type="character", default=NULL, help="Required. Filename containing the (dense) matrix to analyze. First row assumed to be header. Assumes observations are columns and rows are features."),
  make_option(c("-o","--outfile"), type="character", default=NULL, help="Required. Filename for output of the inferred latent factors. A header is included by default"),
  make_option(c("-d","--dims"), type="integer", default=2, help="Desired number of latent dimensions. Default is 2.")
)
opt_parser<-OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)
if(is.null(opt$infile)){
  print_help(opt_parser)
  stop("Input file required.")
}
if(is.null(opt$outfile)){
  print_help(opt_parser)
  stop("Output file required.")
}

Y<-as.matrix(read.table(opt$infile))
suppressWarnings(fit<-zinbwave::zinbFit(Y, K=opt$dims, BPPARAM=BiocParallel::SerialParam()))
factors<-as.data.frame(zinbwave::getW(fit))
colnames(factors)<-paste0("dim",1:opt$dims)
#res<-as.data.frame(zinbwave::zinbFit(Y,K=opt$dims)@W)
#colnames(res)<-paste0("dim",1:opt$dims)
write.table(factors,file=opt$outfile,row.names=FALSE,quot=FALSE)