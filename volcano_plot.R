#######################################################
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot-r/tutorial.html
#########################################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]


#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#library(dplyr)
#library(ggplot2)
#library(ggrepel)
#install.packages('ggrepel')
##########################
filename="AdultvsBaby.txt"
#filename="Comp1_filtered_stat.xls"
A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
colnames(A)

#tiff(paste0(filename,".p05.volcano.tiff"), width = 3000, height = 3000, res=300)
tiff(paste0(filename,".volcano.tiff"), width = 3000, height = 3000, res=300)
EnhancedVolcano(A,
                lab = A[,2],
                x = 'logFC',
                y = 'P.Value',
                #pCutoff=0.05,
                title = filename,
                pointSize = 3.0,
                labSize = 5.0
                )

dev.off()

