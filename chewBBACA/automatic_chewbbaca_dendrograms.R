#!/usr/bin/env Rscript
## --------------------------------------------------------------------------------
## libraries


suppressWarnings(suppressPackageStartupMessages(library(dendextend)))
suppressWarnings(suppressPackageStartupMessages(library(knitr)))
suppressWarnings(suppressPackageStartupMessages(library(zCompositions)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(foreach)))
suppressWarnings(suppressPackageStartupMessages(library(readr)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(matrixStats)))
suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(viridis)))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel)))
suppressWarnings(suppressPackageStartupMessages(library(MASS)))
suppressWarnings(suppressPackageStartupMessages(library(maptpx)))
suppressWarnings(suppressPackageStartupMessages(library(gplots)))
suppressWarnings(suppressPackageStartupMessages(library(grid)))
suppressWarnings(suppressPackageStartupMessages(library(gridExtra)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(ggforce)))
suppressWarnings(suppressPackageStartupMessages(library(Matrix)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(optparse)))

options(expressions = 5e5)

#possible options for input and output
option_list = list(
  make_option("--sampleinfo", type="character", default=NULL, 
              help="Path to file with sample ids and colors", metavar="character"), 
  make_option("--dist", type="character", default=NULL, 
              help="Path to distance matrix", metavar="character"),
  make_option("--out", type="character", default=NULL, 
              help="Output filename", metavar="character"),
  make_option("--titel", type="character", default="", 
              help="Titel for the plot", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##read dist mat from cgmlst-dists
chew<-data.matrix(read.table(opt$dist, row.names = 1, header = TRUE))

##rename columns and rows 
colnames(chew)<-colnames(chew)%>%gsub('^X', '', .)

##build dendrogram 
d_chew<-as.dendrogram(hclust(as.dist(chew), method="single"))

##cluster colors
if(!is.null(opt$sample)){
  sampleInfo<-read_table(opt$sampleinfo)
  labels_colors(d_chew) <- sampleInfo$color[match(labels(d_chew), sampleInfo$sample)]
} else{
  labels_colors(d_chew) <- "black"
}

pdf(opt$out)
plot(d_chew, main=opt$titel)
for (i in 1:25) {
  abline(h=i, col="black", lty='dotted' )
}
dev.off()







