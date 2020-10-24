rm(list=ls())
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################

#setwd("/Users/ping/Desktop/BabySkin/gene_fig")
#filename="Adult_Baby_marker.xls"

truefc<-function(VVV){
	XXX=VVV
	if(VVV==0){
	    XXX=NA
   	}else if(VVV<1){
	    XXX=-1/VVV
    	}
	return(XXX)	
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=log10(A[1:d[1], 5:d[2]])
rownames(B)=A[,2]
Cname=colnames(A)[5:d[2]]
#colnames(A)[1:4]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)
ID=rep(NA,Clen)
for(mm in  1:Clen ){
    Group1[mm]=splitname[[mm]][1]
    Group2[mm]=splitname[[mm]][2]
    Group3[mm]=splitname[[mm]][3]
    ID[mm]=splitname[[mm]][4]
}
#https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
#png(filename=paste0(filename,".png"), height=2800, width=3600, res=300)
#install.packages("pheatmap")
library(pheatmap)
#pheatmap(B)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(B, 1, cal_z_score))
#pheatmap(data_subset_norm)

my_hclust_gene <- hclust(dist(data_subset_norm), method = "complete")
#install.packages("dendextend")
library(dendextend)

#as.dendrogram(my_hclust_gene) %>%
#  plot(horiz = TRUE)
#dev.off()
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
#set.seed(1984)
#my_random <- as.factor(sample(x = 1:2, size = nrow(my_gene_col), replace = TRUE))
#my_gene_col$random <- my_random

#head(my_gene_col)
my_sample_col =data.frame(cbind(Group1,Group2, Group3))
# change order
#my_sample_col$sample <- factor(my_sample_col$Group1, levels = c("Baby", "Adult"))
#change color
#my_colour = list(
#  sample = c(normal = "#5977ff", tumour = "#f74747"),
#  random = c(random1 = "#82ed82", random2 = "#9e82ed"),
#  cluster = c(cluster1 = "#e89829", cluster2 = "#cc4ee0")
#)

row.names(my_sample_col) <- colnames(B)
png(filename=paste0(filename,".png"), height=4800, width=2800, res=120)

p1<-pheatmap(data_subset_norm,
         #annotation_colors = my_colour,
         annotation_row = my_gene_col,
         annotation_col = my_sample_col,
         cutree_rows = 2,
         cutree_cols = 2)
print(p1)
dev.off()

####retrieve result from the heatmap
# use silent = TRUE to suppress the plot
#my_heatmap <- pheatmap(data_subset_norm, silent = TRUE)
# results are stored as a list
#class(my_heatmap)
#[1] "list"
#names(my_heatmap)
#[1] "tree_row" "tree_col" "kmeans"   "gtable"
#my_heatmap$tree_row %>%
#  as.dendrogram() %>%
#  plot(horiz = TRUE)

#save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
#  png(filename, width = width, height = height, res = res)
#  grid::grid.newpage()
#  grid::grid.draw(x$gtable)
#  dev.off()
#}

#save_pheatmap_png(my_heatmap, "my_heatmap.png")