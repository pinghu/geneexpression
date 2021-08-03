rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
##########################
#Ping Hu
#1-25-2021
###################################
filename="GSS2656_Foreskin_Other_signal.txt"


A<-read.table(filename, sep="\t", header=TRUE)
d<-dim(A)
B=log10(A[1:d[1], 2:d[2]])
rownames(B)= paste0(A[,1], ".",A[,2])
Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)
Group4=rep(NA, Clen)
Gender=rep(NA, Clen)
Site=rep(NA, Clen)
GSSID=rep(NA, Clen)
SID=rep(NA, Clen)
for(mm in  1:Clen ){
    Group1[mm]=splitname[[mm]][1]
    Group2[mm]=splitname[[mm]][2]
    Group3[mm]=splitname[[mm]][3]
    Group4[mm]=splitname[[mm]][4]
    Site[mm]=splitname[[mm]][6]
    Gender[mm]=splitname[[mm]][5]
    SID[mm]=splitname[[mm]][7]
    GSSID[mm]=splitname[[mm]][6]
}


Site[Site=="Butt"]="buttocks"
C=t(B)
df_pca <- prcomp(C)

df_out <- as.data.frame(df_pca$x)
df_out$group2=Group2
df_out$group1=Group1
df_out$group3=Group3
df_out$name=Cname

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
df_out$gender=Gender
df_out$site=Site
#install.packages("ggpmisc")
library("ggpmisc")
library('ggpubr')

percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group1, label=name ))+  geom_point(aes(shape=Site,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(Group1)), "Set1"))+theme_bw()+ xlab(percentage[1]) + ylab(percentage[2])+stat_chull(aes(color=Group1, fill=Group1), alpha=0.1, geom="polygon")+geom_text(aes(label=ifelse(SID %in% c(23,30,40, 44,46,67),as.character(SID),'')),hjust=0,vjust=0, color="black")

  #geom_label_repel(aes(label = ifelse(SID %in% c(23,30,40, 44,46,67),as.character(SID),'') ),
  #                 box.padding   = 0.35, 
  #                 point.padding = 0.5,
  #                 segment.color = 'grey50')
  #
png(paste0(filename,".PCA_all.png"), width=800, height=800, res=120)
p1
dev.off()


p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group1, label=name ))+  geom_point(aes(shape=gender, size=5))+scale_color_manual(values=brewer.pal(nlevels(as.factor(Group1)), "Set1"))+theme_bw()+ xlab(percentage[1]) + ylab(percentage[2])+stat_chull(aes(color=Group1, fill=Group1), alpha=0.1, geom="polygon")+geom_text(aes(label=ifelse(SID %in% c(23,30,40, 44,46,67),as.character(SID),'')),hjust=0,vjust=0, color="black")

#geom_label_repel(aes(label = ifelse(SID %in% c(23,30,40, 44,46,67),as.character(SID),'') ),
#                 box.padding   = 0.35, 
#                 point.padding = 0.5,
#                 segment.color = 'grey50')
#
jiff(paste0(filename,".PCA_all_gender.jiff"), width=3000, height=3000, res=300)
p1
dev.off()

p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group1, label=name ))+  geom_point(aes(shape=Site, size=5))+scale_color_manual(values=brewer.pal(nlevels(as.factor(Group1)), "Set1"))+theme_bw()+ xlab(percentage[1]) + ylab(percentage[2])+stat_chull(aes(color=Group1, fill=Group1), alpha=0.1, geom="polygon")+geom_text(aes(label=ifelse(SID %in% c(23,30,40, 44,46,67),as.character(SID),'')),hjust=0,vjust=0, color="black")

#geom_label_repel(aes(label = ifelse(SID %in% c(23,30,40, 44,46,67),as.character(SID),'') ),
#                 box.padding   = 0.35, 
#                 point.padding = 0.5,
#                 segment.color = 'grey50')
#
png(paste0(filename,".PCA_all_site.png"), width=800, height=800, res=120)
p1
dev.off()

#############################all things below not needed#################################################



df_pca2 <- prcomp(C[Group2 !="NotUse",])
df_out2 <-as.data.frame(df_pca2$x)
group=Group2[Group2 !="NotUse"]
name2=Cname[Group1 !="InfantOutlier"]
p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(aes(shape=BodySite,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()
png(paste0(filename,".PCA_NoOutlier.png"), width=800, height=800, res=120)
p1
dev.off()

p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(aes(shape=BodySite,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_bw()
png(paste0(filename,".PCA_NoOutlier_BW.png"), width=800, height=800, res=120)
p1
dev.off()

df_pca2 <- prcomp(C[Group1 !="InfantOutlier" & Group1 !="Infant_foreskin",])
df_out2 <-as.data.frame(df_pca2$x)
group=Group1[Group1 !="InfantOutlier" & Group1 != "Infant_foreskin"]
name2=Cname[Group1 !="InfantOutlier"& Group1 != 'Infant_foreskin']
p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()
png(paste0(filename,".PCA_NoForeskin.png"), width=800, height=800, res=120)
p1
dev.off()



GroupAB=group
GroupAB[GroupAB=='Adult20s'| GroupAB=='Adult60s']='Adult'

p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=GroupAB, label=name2 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(GroupAB)), "Set1"))+stat_chull(aes(color=GroupAB, fill=GroupAB), alpha=0.1, geom="polygon") +theme_pubclean()
png(paste0(filename,".PCA_NoForeskin_AB.png"), width=800, height=800, res=120)
p1
dev.off()

BodySite=Site[Group1 !="InfantOutlier" & Group1 != "Infant_foreskin"]
gender=Gender[Group1 !="InfantOutlier" & Group1 != "Infant_foreskin"]
percentage <- round(df_pca2$sdev / sum(df_pca2$sdev) * 100, 2)
percentage <- paste( colnames(df_out2), "(", paste( as.character(percentage), "%", ")", sep="") )

#p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
#p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])
p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(aes(shape=BodySite,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()+ xlab(percentage[1]) + ylab(percentage[2])
tiff(paste0(filename,".PCA_NoForeskin_gender_bodysite_percent.tiff"), width=2800, height=1600, res=300)
p1
dev.off()

p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(aes(shape=BodySite,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()
tiff(paste0(filename,".PCA_NoForeskin_gender_bodysite.tiff"), width=2800, height=1600, res=300)
p1
dev.off()

p1<-ggplot(df_out2,aes(x=PC3,y=PC4,color=group, label=name2 ))+ geom_point(aes(shape=BodySite,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()+ xlab(percentage[3]) + ylab(percentage[4])
png(paste0(filename,".PCA_NoForeskin_gender_bodysite_percent_PC34.png"), width=2800, height=1600, res=300)
p1
dev.off()

p1<-ggplot(df_out2,aes(x=PC5,y=PC6,color=group, label=name2 ))+ geom_point(aes(shape=BodySite,  size=gender))+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()+ xlab(percentage[5]) + ylab(percentage[6])
png(paste0(filename,".PCA_NoForeskin_gender_bodysite_percent_PC56.png"), width=2800, height=1600, res=300)
p1
dev.off()
#######################above For Paper 1 Cluster Figure #######################


df_pca3 <- prcomp(C[Group1 !="InfantOutlier" & Group1 !="Adult20s" & Group1 != "Adult60s",])
df_out3 <-as.data.frame(df_pca3$x)
group=Group1[Group1 !="InfantOutlier" & Group1 != "Adult20s" & Group1 != "Adult60s"]
group[group=="Infant"]="Infant_non-foreskin"
name3=Cname[Group1 !="InfantOutlier"& Group1 != 'Adult20s' & Group1 != "Adult60s"]
p1<-ggplot(df_out3,aes(x=PC1,y=PC2,color=group, label=name3 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()
png(paste0(filename,".PCA_baby_Foreskin.png"), width=800, height=800, res=120)
p1
dev.off()










p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_bw()
png(paste0(filename,".PCA_NoOutlier_BW.png"), width=800, height=800, res=120)
p1
dev.off()






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
