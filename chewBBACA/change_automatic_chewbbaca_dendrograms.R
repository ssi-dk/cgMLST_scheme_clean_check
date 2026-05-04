#!/usr/bin/env Rscript
## --------------------------------------------------------------------------------
## libraries


suppressWarnings(suppressPackageStartupMessages(library(dendextend)))
suppressWarnings(suppressPackageStartupMessages(library(knitr)))
suppressWarnings(suppressPackageStartupMessages(library(zCompositions)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
#suppressWarnings(suppressPackageStartupMessages(library(grDevices)))
#suppressWarnings(suppressPackageStartupMessages(library(propr)))
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


#possible options for input and output
option_list = list(
  make_option("--sampleinfo", type="character", default=NULL, 
              help="Path to file with sample ids and colors", metavar="character"), 
  make_option("--dist", type="character", default=NULL, 
              help="Path to distance matrix", metavar="character"),
  make_option("--out", type="character", default=NULL, 
              help="Output filename", metavar="character"),
  make_option("--height", type="integer", default=NULL, 
              help="Height of plot, default 8", metavar="integer"),
  make_option("--width", type="integer", default=NULL, 
              help="Width of plot, default 12", metavar="integer"),
  make_option("--titel", type="character", default="", 
              help="Titel for the plot", metavar="character"),
  make_option("--lab", type="numeric", default=3, 
              help="Text size for sample labels, default 3", metavar="number"),
  make_option("--ylim", type="integer", default=100, 
              help="max y-axis, default 100", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



#setwd("//sshfs.r/niesof@login.ugerm.dksund.dk/users/data/Projects/FBI_species/proj/yersinia/cgMLST_test/chewBBACA")

##read dist mat from cgmlst-dists
#chew<-data.matrix(read.table("dist.mat.tsv", row.names = 1, header = TRUE, check.names = FALSE))
chew<-data.matrix(read.table(opt$dist, row.names = 1, header = TRUE, check.names = FALSE))


##rename columns and rows 
colnames(chew)<-colnames(chew)%>%gsub('X', '', .)

##build dendrogram 
d_chew<-as.dendrogram(hclust(as.dist(chew), method="single"))

#as.dendrogram(hclust(as.dist(chew), method="single"))%>%head(n=60) #as text with differences written 

##cluster colors
#sampleInfo<-read_table("//sshfs.r/niesof@login.ugerm.dksund.dk/users/data/Projects/FBI_species/proj/yersinia/reads/sampleInfo.tsv")
if(!is.null(opt$sample)){
  sampleInfo<-read_table(opt$sampleinfo)
  labels_colors(d_chew) <- sampleInfo$color[match(labels(d_chew), sampleInfo$sample)]
} else{
  labels_colors(d_chew) <- "black"
}

#change labels 
#labels(d_chew) <- sampleInfo$`ST-sample`[match(labels(d_chew), sampleInfo$sample)]


pdf(opt$out, width=opt$width, height=opt$height)
par(mar=c(7,4,4,2), cex.lab=opt$lab, cex=opt$lab)
plot(d_chew, main=opt$titel, ylim=c(0,opt$ylim), panel.first={
  for (i in 1:25) {
    abline(h=i, col="cornsilk3", lty='dotted')
  }
})
dev.off()


# 
# pdf("dendrogram_ST.pdf", height=60, width=160)
# par(mar=c(7,4,4,2), cex.lab=1.5, cex=1.5)
# plot(d_chew, main="Yersinia cgMLST", ylim=c(0,150), panel.first={
#   for (i in 1:25) {
#     abline(h=i, col="cornsilk3", lty='dotted')
#   }
# })
# dev.off()





