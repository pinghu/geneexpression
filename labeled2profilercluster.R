#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
#http://yulab-smu.top/clusterProfiler-book/chapter12.html#upset-plot
#https://f1000research.com/articles/5-1384/v1
rm(list=ls())
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#filename="PT3vsFT3.txt"
FCCUT=1.2

log2FC_truefc<-function(VVV){
  XXX=2**VVV
  if(VVV==0){
    XXX=NA
  }else if(VVV<1){
    XXX=-1/VVV
  }
  return(XXX)	
}
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
  write.table(summary(ego_bpU), file=paste0(filename,".", NAME,".",GOLEVEL), sep="\t", row.name=TRUE
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
  write.table(summary(kk), file=paste0(filename,".",NAME,".kk"), sep="\t", row.name=TRUE
              , col.names = TRUE, eol = "\n", na = "NA")
}

A<-read.csv(filename, sep="\t", header=TRUE)
dim(A)

library("hgu219.db")
colnames(A)[1]='probe_id'
egids2=hgu219ENTREZID[A$probe_id]
annots=toTable(egids2)
A$truefc<-log2FC_truefc(A$logFC)
x<- merge(A,annots,by="probe_id")
rm(annots)
rm(A)
rm(egids2)
library(dplyr)
x1 <- x %>% group_by(gene_id) %>% filter(AveExpr == max(AveExpr)) 
#[1] "probe_id"   "GENESYMBOL" "GENENAME"   "EntrezID"   "logFC"      "AveExpr"    "t"          "P.Value"    "adj.P.Val" 
#[10] "B"          "FC"         "TrueFC"     "X.TrueFC."  "truefc"     "gene_id"
###significant gene list
x2 <- x1  %>% filter(abs(TrueFC)>=FCCUT) 
dim(x2) ###now only 566 genes

genelist<-x2$TrueFC
names(genelist) <-x2$gene_id
genelist <- sort(genelist, decreasing = TRUE)
library("DOSE")
data(geneList, package="DOSE")
foldchanges=x2$TrueFC
names(foldchanges)=x2$gene_id
gene<-names(foldchanges)
geneU<-names(foldchanges [foldchanges >=FCCUT]) ###259
geneD<-names(foldchanges [foldchanges <=-FCCUT]) ###307
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
write.table(summary(mkk), file=paste0(filename,".mkk"), sep="\t", row.name=TRUE
            , col.names = TRUE, eol = "\n", na = "NA")

mkk1<- setReadable(mkk, 'org.Hs.eg.db', 'ENTREZID')
jpeg(paste0(filename,".mkk_cnet.jpg"), width = 800, height = 800)
p7<-cnetplot(mkk1, foldChange=genelist
         #, circular=TRUE, colorEdge = TRUE
         )
print(p7)
dev.off()

library(meshes)
mesh <- enrichMeSH(gene, MeSHDb = "MeSH.Hsa.eg.db", database='gendoo', category = 'C')
write.table(summary(mesh), file=paste0(filename,".mesh"), sep="\t", row.name=TRUE
            , col.names = TRUE, eol = "\n", na = "NA")


library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
rst = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
write.table(rbind(rst$greater,rst$less), file=paste0(filename, ".gage_kegg"),sep="\t") 

keggrespathways = data.frame(id=rownames(rst$greater), rst$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()

print("Top Greater Pathway");
keggrespathways 
Lesspathways = data.frame(id=rownames(rst$less), rst$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
print ("Top Less Pathway")
Lesspathways 
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
write.table(gobpres, file=paste0(filename,".gage_gobp"),sep="\t")
goccsets = go.sets.hs[go.subs.hs$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
write.table(gobpres, file=paste0(filename,".gage_gocc"),sep="\t")
gomfsets = go.sets.hs[go.subs.hs$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE)
write.table(gobpres, file=paste0(filename,".gage_gomf"),sep="\t")
library(STRINGdb)
string_db <- STRINGdb$new( version="11", species=9606, 
                           score_threshold=200, input_directory="")
test<-cbind(x2$P.Value, x2$TrueFC, x2$GENESYMBOL)
colnames(test)=c("pvalue", "TrueFC", "gene")
example1_mapped <- string_db$map(test, "gene", removeUnmappedRows = TRUE )
UP<-example1_mapped[example1_mapped$TrueFC>=FCCUT,]
Down<-example1_mapped[example1_mapped$TrueFC<=-FCCUT,]

hitU <- UP$STRING_id
jpeg(paste0(filename,"UP.stringdb.jpg"), width = 800, height = 800)
string_db$plot_network( hitU )
dev.off()
enrichment <- string_db$get_enrichment( hitU )
write.table(enrichment, file=paste0(filename,".stringdb_up"),sep="\t")

hitD <- Down$STRING_id
jpeg(paste0(filename,"Down.stringdb.jpg"), width = 800, height = 800)
string_db$plot_network( hitD )
dev.off()
enrichment <- string_db$get_enrichment( hitD )
write.table(enrichment, file=paste0(filename,".stringdb_down"),sep="\t")
