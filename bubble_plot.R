#title: "BabyComplexHeatmap for InfantPaper2"
#author: "Ping hu"
#date: "1/3/2022"

rm(list=ls())


filename="CC_combine.xls"
A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
colnames(A)
library(forcats)
A$Description <- as.factor(A$Description)
A$Description <- fct_inorder(A$Description)
library(ggplot2)
p1<-ggplot(A, aes(compare, Description)) +
  geom_point(aes(color=p.adjust, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))

p2<-ggplot(A, aes(compare, Description)) +
  geom_point(aes(color=p.adjust, size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))+
  theme(legend.direction = "horizontal", legend.position = "top")+#legend??????????????????
  scale_y_discrete(position = "right")#y??????????????????


p3<-ggplot(A, aes(compare, Description)) +
  ####bluecolor
  geom_point(aes(fill=p.adjust, size=Count), color="#999999", shape=21)+theme_bw()+
  scale_color_gradient(low='#99FFFF',high='#0033FF')+
  #scale_fill_manual((values=c("#212c5f", "#3366b1","#42b0e4", "#7bc6ed", "dfe1e0")))+
  geom_point(data=A[which(A$direction=="UP"),], aes(color=p.adjust, size=Count), shape=16)+
  scale_color_gradient(low='#FFFF00',high='#FF3300')+
  theme_bw()+
  #scale_fill_manual((values=c("#212c5f", "#3366b1","#42b0e4", "#7bc6ed", "dfe1e0")))+
  theme(panel.grid.minor.x=element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))+
  theme(legend.direction = "horizontal", legend.position = "top")+#legend??????????????????
  scale_y_discrete(position = "right")#y??????????????????

######https://www.jianshu.com/p/cdc7d9d95fb8
library(tidyverse)
library(ggpubr)
cluster <- A$Description %>% as.data.frame() %>%
  mutate(group=rep(c("cluster1","cluster2","cluster3"),
                   times=c(13,12,13))) %>%
  mutate(p="")%>%
  ggplot(aes(p,.,fill=group))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5))+
  labs(fill = "Pathway")

A$fold <- ''
A$fold <- -log(A$GeneRatio)

p1 <-ggplot(A, aes(Group, Description)) +
  geom_point(aes(color=fold, size=LOGP))+theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_y_discrete(position = "right")
p1%>%
  insert_left(cluster, width = .05)

print(p2)

jpeg(paste0(filename,".jpg"), width = 2800, height = 1800, res=300)

dev.off()
########################################################
