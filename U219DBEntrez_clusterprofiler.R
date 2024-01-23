###################################################################
#########author: Ping Hu
#########Date: May11 2021
#########This script will use the original file's EntrezGeneId to do mapping 
#########also use provided true fold change 
#########need to compare this to the script use the U219.db #######
#########only need to provide the file contain both up and down gene lists
####################################################################
#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
#http://yulab-smu.top/clusterProfiler-book/chapter12.html#upset-plot
#https://f1000research.com/articles/5-1384/v1
rm(list=ls())
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

filename="VPT_FT"
FCCUT=1.2

###GOLEVEL will be "BP", "CC", "MK
GOtest<-function(GENE, NAME, GOLEVEL){
  ego_bpU <- enrichGO(gene          = GENE,
                      universe      = allgenes,
                      OrgDb         = org.Hs.eg.db,
                      ont           = GOLEVEL,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      minGSSize = 5,
                      maxGSSize = 1000,
                      readable      = TRUE)
  ###was able to detect the cholestrol and antigen presenting with max 300 min 5
  jpeg(paste0(filename,".", NAME,".ego",GOLEVEL,".dot.jpg"), width = 800, height = 800)
  p1<-dotplot(ego_bpU,showCategory = 50 )+ ggtitle(paste(filename, NAME, GOLEVEL))
  print(p1)
  dev.off()
  jpeg(paste0(filename,".", NAME,".ego",GOLEVEL,"bar.jpg"), width = 800, height = 800)
  p2<-barplot(ego_bpU,showCategory = 50)+ ggtitle(paste(filename, NAME, GOLEVEL))
  print(p2)
  dev.off()
  write.table(as.data.frame(ego_bpU), file=paste0(filename,".", NAME,".",GOLEVEL), sep="\t", row.name=TRUE
              , col.names = TRUE, eol = "\n", na = "NA")
}

KEGGtest<-function(GENE, NAME){
  kk <- enrichKEGG(gene         = GENE,
                   organism     = 'hsa',
                   minGSSize = 5,
                   maxGSSize = 1000,
                   pvalueCutoff = 0.05)
  
  jpeg(paste0(filename,".",NAME,".kk_dot.jpg"), width = 800, height = 800)
  p3<-dotplot(kk,showCategory = 50 )+ ggtitle(paste(filename,NAME, "kegg"))
  print(p3)
  dev.off()
  jpeg(paste0(filename,".",NAME,".kk_bar.jpg"), width = 800, height = 800)
  p4<-barplot(kk,showCategory = 50)+ ggtitle(paste(filename, NAME, "kegg"))
  print(p4)
  dev.off()
  write.table(as.data.frame(kk), file=paste0(filename,".",NAME,".kk"), sep="\t", row.name=TRUE
              , col.names = TRUE, eol = "\n", na = "NA")
}

A<-read.csv(filename, sep="\t", header=TRUE)
dim(A)

library("hgu219.db")
comparison=colnames(A)[1]
colnames(A)[1]='probe_id'

library(dplyr)
dim(A)
x0<-A %>% filter(abs(P.Value)<=0.05) 
dim(x0)
x1 <- x0  %>% filter(abs(Truefc.from.logFC)>=FCCUT)
dim(x1)
x2 <- x1 %>% group_by(EntrezID) %>% filter(AveExpr == max(AveExpr)) 

#[1] "probe_id"          "GENESYMBOL"        "GENENAME"          "EntrezID"          "logFC"            
#[6] "AveExpr"           "t"                 "P.Value"           "adj.P.Val"         "B"                
#[11] "Truefc.from.logFC" "absTruefc.logFC"

###significant gene list
dim(x2)
###now only 659
#genelist<-x2$Truefc.from.logFC
#names(genelist) <-x2$EntrezID
#genelist <- sort(genelist, decreasing = TRUE)
library("DOSE")
data(geneList, package="DOSE")
foldchanges=x2$Truefc.from.logFC
names(foldchanges)=x2$EntrezID
foldchanges <- sort(foldchanges, decreasing = TRUE)
gene<-names(foldchanges)
geneU<-names(foldchanges [foldchanges >=FCCUT]) ###259
geneD<-names(foldchanges [foldchanges <=-FCCUT]) ###307
allgenes <-unique(unlist(as.list(hgu219ENTREZID))) ##18465 
write.table(x1, file=paste0(filename,".unique_gene"),sep="\t")
write.table(x2, file=paste0(filename,".fc1D2"),sep="\t")

library(magrittr)
library(clusterProfiler)
library(tidyr)
library("hgu219.db")
GOtest(gene,"All","BP")
GOtest(gene,"All","CC")
GOtest(gene,"All","MF")
GOtest(geneU,"Up","BP")
GOtest(geneU,"Up","CC")
GOtest(geneU,"Up","MF")
GOtest(geneD,"Down","BP")
GOtest(geneD,"Down","CC")
GOtest(geneD,"Down","MF")
KEGGtest(gene,"All")
KEGGtest(geneU,"Up")
KEGGtest(geneD,"Down")

