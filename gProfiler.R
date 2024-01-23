#####################################################################
#Example code from Ping Hu May 8th 2023
#multiple methods for the gene functional analysis  -- need to display 
#results in different tabs, also allow user to select whether to use fold 
#change or p value to change the gene list, if the interactive figures are available, 
#we should switch to use interactive figures for example gprofiler or vocalnoplot.
# multiple databases should be all plotted against such as gage and profilercluster 
#with GO CC BP MF and KEGG, reactome and diseases ... need to display them nicely
##
########################################################################
rm(list=ls())
# Read the CSV file into a data frame
data_file <- read.csv("bj1080.11.bj.txt.stat")
ann_file <- read.csv("HG-U219.na36.ann.csv")
colnames(ann_file)[1]="ProbeID"
colnames(data_file)[1]="ProbeID"
my_data <- merge(data_file, ann_file, by="ProbeID")
dim(my_data)
# Print the first few rows of the data frame
#head(my_data)

library(ggplot2)

p_values <- my_data$p_values
FoldChange <- my_data$FoldChange

sig_level <- 0.05

log2_FC <- log2(FoldChange)

# Calculate the negative log10 p-values
neg_log10_p <- -log10(p_values)

# Create a data frame with the log2 fold change and negative log10 p-values
df <- data.frame(log2_FC, neg_log10_p)
#####################################
# part 1. Generate the volcano plot
#######################################
ggplot(df, aes(x=log2_FC, y=neg_log10_p)) + 
  geom_point(aes(color=ifelse(p_values < sig_level & abs(log2_FC) > 1, "red", "black")), alpha=0.5) +
  geom_hline(yintercept=-log10(sig_level), linetype="dashed") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  scale_color_identity() +
  labs(x="log2(Fold Change)", y="-log10(p-value)", color="Significant") +
  theme_classic()

library(dplyr)
#sig_data <- my_data %>% filter(p_values <= sig_level& abs(log2_FC) > 0.8)
sig_data <- my_data %>% filter(p_values <= sig_level)
sig_probes <- sig_data$GeneSymbol

######################################################################
#part 2. generate the gprofile result table and interactive figures
######################################################################
library(gprofiler2)
#gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), organism = "hsapiens")
gostres <- gost(query = sig_probes, organism = "hsapiens")

# The result is a named list where "result" is a data.frame with the enrichment analysis results
# and "meta" containing a named list with all the metadata for the query.
head(gostres$result)

#p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p <- gostplot(gostres, capped = FALSE, interactive = TRUE)
p

#gconvert enables to map between genes, proteins, microarray probes, common names, various database identifiers, etc, from numerous databases and for many species.
#gconvert(query = c("REAC:R-HSA-3928664", "rs17396340", "NLRP1"), organism = "hsapiens",target="ENSG", mthreshold = Inf, filter_na = TRUE)
#gorth translates gene identifiers between organisms. For example, to convert gene list between mouse (source_organism = mmusculus) and human (target_organism = hsapiens):
#gorth(query = c("Klf4", "Sox2", "71950"), source_organism = "mmusculus",target_organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC")
#gsnpense converts a list of SNP rs-codes (e.g. rs11734132) to chromosomal coordinates, gene names and predicted variant effects. Mapping is only available for variants that overlap with at least one protein coding Ensembl gene.
#gsnpense(query = c("rs11734132", "rs4305276", "rs17396340", "rs3184504"))
#library(enrichR)
#enrich_result <- enrichr(signature="GO_Biological_Process_2018", genes=sig_probes)

#############################
#Part 3 StringDB result
#############################
library(STRINGdb)
string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=200, input_directory="")
STRINGdb$methods()
example1_mapped <- string_db$map( sig_data, "GeneSymbol", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id
par(mar=c(4, 4, 2, 2))
string_db$plot_network( hits )
string_db$get_enrichment( hits )

#############################
#part 4 . profileCluster
#############################
library(magrittr)
library(clusterProfiler)
library(tidyr)
library("hgu219.db")
library("DOSE")
#sig_data <- sort(FoldChange, decreasing = TRUE)
gene<-sig_data$EntrezID
###up-regulated or down-regulated genes
geneU<- sig_data[sig_data$FoldChange >1,]$EntrezID
geneD<-sig_data[sig_data$FoldChange <1,]$EntrezID
allgenes <-unique(unlist(as.list(merged_data$EntrezID)))

###GOLEVEL will be "BP", "CC", "MK", you might need to 
GOtest<-function(GENE, NAME, GOLEVEL, filename){
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

KEGGtest<-function(GENE, NAME, filename){
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
GOtest(gene,"All","BP", "test")
GOtest(gene,"All","CC", "test")
GOtest(gene,"All","MF", "test")
GOtest(geneU,"Up","BP", "test")
GOtest(geneU,"Up","CC", "test")
GOtest(geneU,"Up","MF", "test")
GOtest(geneD,"Down","BP", "test")
GOtest(geneD,"Down","CC", "test")
GOtest(geneD,"Down","MF", "test")

####is the kegg download been blocked by company firewall?
KEGGtest(gene,"All", "test")
KEGGtest(geneU,"Up", "test")
KEGGtest(geneD,"Down", "test")
####Pathway enrichment analysis
x<-enrichWP(gene, organism = "Homo sapiens") 
head(x)
library(ReactomePA)
x <- enrichPathway(gene=gene, pvalueCutoff = 0.05, readable=TRUE)
head(x)
#####disease enrichment analysis
x <- enrichDO(gene          = gene,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = allgenes,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)
#10.1.3 Over-representation analysis for the disease gene network
dgn <- enrichDGN(gene) 
head(dgn)
####MESH analysis
#library(AnnotationHub)
#library(MeSHDbi)
#ah <- AnnotationHub(localHub=TRUE)
#hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
#file_hsa <- hsa[[1]]
#db <- MeSHDbi::MeSHDb(file_hsa)
#x <- enrichMeSH(gene, MeSHDb = db, database='gendoo', category = 'C')
####################
#Part 5 gage method
#####################
data <- read.csv(file = 'bj1080.11.bj.txt.stat')
head(data)
require('AnnotationDbi')
require("BiocGenerics")
require("parallel")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu219.db")

chip='hgu219'
cmd=paste('library(',chip,'.db)',sep='')  
eval(parse(text=cmd))

colnames(data)[1]='probe_id'
cmd=paste('egids2=', chip,'ENTREZID[data$probe_id]',sep='')
eval(parse(text=cmd))
annots=toTable(egids2)
total <- merge(data,annots,by="probe_id")

library(dplyr)
x<- total %>% group_by(gene_id) %>% filter(p_values == min(p_values)) 
remove(total)
remove(data)

library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)


Sig<-x[x$p_values <=0.05,]
foldchanges=Sig$FoldChange
names(foldchanges)=Sig$gene_id

keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)
write.table(keggres, file="test.kegg.txt",sep="\t")

data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
write.table(gobpres, file="test.gobp.txt",sep="\t")

goccsets = go.sets.hs[go.subs.hs$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
write.table(goccres, file="test.gocc.txt",sep="\t")

gomfsets = go.sets.hs[go.subs.hs$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE)
write.table(gomfres, file="test.gomf.txt",sep="\t")
########################
#part 6 fgsea analysis
#######################
library(fgsea)
library(data.table)
library(ggplot2)
data(examplePathways)
data(exampleRanks)
set.seed(42)
fold_all <- x$FoldChange
names(fold_all)=x$gene_id
fold_sorted_desc <- sort(fold_all, decreasing = TRUE)
fold_sorted <- sort(fold_all)
####reactome pathway analysis without any arbituary cut
pathways <- reactomePathways(names(fold_all))
fgseaRes_1 <- fgsea(pathways, fold_sorted, maxSize=500)
fgseaRes_2 <- fgsea(pathways, fold_sorted_desc, maxSize=500)
head(fgseaRes_2)
head(fgseaRes_1)

