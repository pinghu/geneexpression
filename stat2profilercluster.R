#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
#http://yulab-smu.top/clusterProfiler-book/chapter12.html#upset-plot
#https://f1000research.com/articles/5-1384/v1
rm(list=ls())
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
filename ="comp4_p05_filter"
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
colnames(A)[1]='probe_id'
colnames(A)[2]="pvalue"
colnames(A)[3]="truefc"
egids2=hgu219ENTREZID[A$probe_id]
annots=toTable(egids2)
x<- merge(A,annots,by="probe_id")
rm(annots)
rm(A)
rm(egids2)
library(dplyr)
###only choose single gene
#x1 <- x %>% group_by(gene_id) %>% filter(AveExpr == max(AveExpr))
x1 <- x %>% group_by(gene_id) %>% filter(abs(truefc) == max(abs(truefc)))
#[1] "probe_id"   "GENESYMBOL" "GENENAME"   "EntrezID"   "logFC"      "AveExpr"    "t"          "P.Value"    "adj.P.Val" 
#[10] "B"          "FC"         "TrueFC"     "X.TrueFC."  "truefc"     "gene_id"
###significant gene list
x2 <- x1  %>% filter(abs(truefc)>=FCCUT) 
dim(x2) ###now only 566 genes there are 1358 genes

genelist<-x2$truefc
names(genelist) <-x2$gene_id
genelist <- sort(genelist, decreasing = TRUE)
library("DOSE")
data(geneList, package="DOSE")
foldchanges=x2$truefc
names(foldchanges)=x2$gene_id
gene<-names(foldchanges)
geneU<-names(foldchanges [foldchanges >=FCCUT]) ###259
print("geneU")
length(geneU)
geneD<-names(foldchanges [foldchanges <=-FCCUT]) ###307
print("geneD")
length(geneD)
allgenes <-unique(unlist(as.list(hgu219ENTREZID))) ##18465 
write.table(x1, file=paste0(filename,".unique_gene"),sep="\t")
write.table(x2, file=paste0(filename,".fc1D2"),sep="\t")

library(magrittr)
library(clusterProfiler)
library(tidyr)
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

####cholestrol come out, can I do the mapping
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')
jpeg(paste0(filename,".mkk_dot.jpg"), width = 800, height = 800)
p5<-dotplot(mkk,showCategory = 50 )+ ggtitle(paste(filename, "mkk"))
print(p5)
dev.off()
jpeg(paste0(filename,".mkk_bar.jpg"), width = 800, height = 800)
p6<-barplot(mkk,showCategory = 50)+ ggtitle(paste(filename, "mkk"))
print(p6)
dev.off()
write.table(as.data.frame(mkk), file=paste0(filename,".mkk"), sep="\t", row.name=TRUE
            , col.names = TRUE, eol = "\n", na = "NA")
#######error message from here####for 4 and 5 comp
#Error in graph.data.frame(x, directed = FALSE) : 
#  the data frame should contain at least two columns
mkk1<- setReadable(mkk, 'org.Hs.eg.db', 'ENTREZID')
jpeg(paste0(filename,".mkk_cnet.jpg"), width = 800, height = 800)
p7<-cnetplot(mkk1, foldChange=genelist
         #, circular=TRUE, colorEdge = TRUE
         )
print(p7)
dev.off()
