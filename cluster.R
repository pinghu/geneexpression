rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################
#install.packages("heatmaply")
#library("heatmaply")
#install.packages('autoplotly')
#library("autoplotly")
filename="GSS2656_data_wAnn"
#p <- autoplotly(prcomp(iris[c(1, 2, 3, 4)]), data = iris,
#                colour = 'Species', frame = TRUE)


A<-read.table(filename, sep="\t", header=TRUE)
d<-dim(A)
B=log10(A[1:d[1], 4:d[2]])
rownames(B)= paste0(A[,1], ".",A[,2])
Cname=colnames(A)[4:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)

GSSID=rep(NA, Clen)
for(mm in  1:Clen ){
    Group1[mm]=splitname[[mm]][1]
    Group2[mm]=splitname[[mm]][2]
    Group3[mm]=splitname[[mm]][3]
    GSSID[mm]=splitname[[mm]][4]
}

C=t(B)
df_pca <- prcomp(C)

df_out <- as.data.frame(df_pca$x)
df_out$group2=Group2
df_out$group1=Group1
df_out$group3=Group3
df_out$name=Cname
#library("autoplotly")
#iris[c(1, 2, 3, 4)]
#prcomp(iris[c(1, 2, 3, 4)])
#p <- autoplotly(df_out$PC1, df_out$PC2, data = df_out,
#                colour = group1, frame = TRUE)
#p

rm(A)

library(ggplot2)
library(grid)
library(gridExtra)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


####Generate a dot plot##########
library(RColorBrewer)

####with label
df_pca <- prcomp(C)
df_out <- as.data.frame(df_pca$x)
df_out$group2=Group2
df_out$group1=Group1
df_out$group3=Group3
df_out$name=Cname
install.packages("ggpmisc")
library("ggpmisc")
p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group1, label=name ))+geom_point(size=6)+ scale_color_manual(values=brewer.pal(nlevels(as.factor(Group1)), "Set1")) +stat_chull(aes(color=Group1, fill=Group1), alpha=0.1, geom="polygon") #+ stat_ellipse()
p1

df_pca2 <- prcomp(C[Group2 !="NotUse",])
df_out2 <-as.data.frame(df_pca2$x)
group2=Group2[Group2 !="NotUse"]
name2=Cname[Group2 !="NotUse"]

p2<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group2, label=name2 ))+geom_point(size=6) + scale_color_manual(values=brewer.pal(nlevels(as.factor(group2)), "Set2")) + stat_chull(aes(color=group2, fill=group2), alpha=0.01, geom="polygon") #+ stat_ellipse()

p2

df_pca3 <- prcomp(C[Group3 !="NotUse",])
df_out3 <-as.data.frame(df_pca3$x)
group3=Group3[Group3 !="NotUse"]
name3=Cname[Group3 !="NotUse"]
p3<-ggplot(df_out3,aes(x=PC1,y=PC2,color=group3, label=name3 )) +geom_point(size=6)+ scale_color_manual(values=brewer.pal(nlevels(as.factor(group3)), "Dark2")) + stat_chull(aes(color=group3, fill=group3), alpha=0.01, geom="polygon")
p3

library("ggpubr")
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
figure <- ggarrange(p1, p2, p3,
                    labels = c("A", "B", "C", "D"),
                    ncol = 3, nrow = 1)
png(paste0(filename,".PCA.png"), width=1800, height=800, res=120)
figure
dev.off()



library(RColorBrewer)
library(limma)
library(Glimma)
col.group <- as.factor(paste0(Group1, Group2, Group3))
levels(col.group) <-  brewer.pal(nlevels(col.group), "Dark2")
col.group <- as.character(col.group)

#plotMDS(B, labels=Group, col=col.group) Generate the web page for MDS data
#https://rdrr.io/bioc/Glimma/man/glMDSPlot.default.html abiut other options
glMDSPlot(B, labels=Cname, groups=paste0(Group1,Group2,Group3), launch=FALSE, main=paste0(filename, ".MDS_Plot"),html=paste0(filename, ".MDS_Plot"), folder = paste0(filename,".MDS_Plot") )

library(factoextra)
library(NbClust)
#kmeans, pam, clara and hcut (for hierarchical clustering).
c=t(B)
fviz_nbclust(c, kmeans, method = c("silhouette", "wss", "gap_stat"))
#https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
  # Elbow method
  jpeg(paste0(filename, ".kmeans_elbow.jpg"), width = 350, height = 350)
  fviz_nbclust(c, kmeans, method = "wss") +
      geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")
  dev.off()

  # Silhouette method
  jpeg(paste0(filename,".kmeans_Silhouette.jpg"), width = 350, height = 350)
  fviz_nbclust(c, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method")
  dev.off()

  # Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.

  jpeg(paste0(filename,".kmeans_Gap_Statistic.jpg"), width = 350, height = 350)
  set.seed(123)
  fviz_nbclust(c, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
    labs(subtitle = "Gap statistic method")
  dev.off()
  #data: matrix
#diss: dissimilarity matrix to be used. By default, diss=NULL, but if it is replaced by a dissimilarity matrix, distance should be “NULL”
#distance: the distance measure to be used to compute the dissimilarity matrix. Possible values include “euclidean”, “manhattan” or “NULL”.
#min.nc, max.nc: minimal and maximal number of clusters, respectively
  #method: The cluster analysis method to be used including “ward.D”, “ward.D2”, “single”, “complete”, “average”, “kmeans” and more.
  jpeg(paste0(filename, ".NbClust_euclidean_2_15.jpg"), width = 350, height = 350)
  NbClust(data = c, distance = "euclidean",
          min.nc = 2, max.nc = 15, method = kmeans)
  dev.off()

#once determin how many cluster to use, then
  set.seed(123)
  km.res <- kmeans(c, 3, nstart = 25)
  jpeg(paste0(filename,".3kmeans.jpg"), width = 350, height = 350)
  fviz_cluster(km.res, data = c,
               ellipse.type = "convex",
               palette = "jco",
               ggtheme = theme_minimal())
  dev.off()
