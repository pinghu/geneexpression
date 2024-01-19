rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
filename="GSS3049.expected_count.newID2"
A <- read.table(filename, header = TRUE,  sep="\t")

#cnts=read.csv(filename,sep='\t')

d=dim(A)
B=A[1:d[1], 3:d[2]]
#ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
#C=B+ZZ
Cname=colnames(B)
Clen=length(Cname) ##there are 3 annotation columns
CID=rep("NA", Clen)
Chem=rep("NA", Clen)
Dose=rep("NA",Clen)
ChemDose=rep("NA", Clen)
dupID=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
  CID[mm]=splitname[[mm]][1]
  Chem[mm]=splitname[[mm]][2]
  Dose[mm]=splitname[[mm]][3]
  dupID[mm]=splitname[[mm]][4]
}
ChemDose=paste0(CID, ".", Chem, Dose)
SeqID=paste0(CID, ".", dupID)
geneID=A[,1]

samplenames <- Cname

group <- as.factor(ChemDose)
samplenames=paste0(ChemDose, ".", SeqID)
cnts = B
colnames(cnts)=samplenames
rownames(cnts)=geneID
meta<- data.frame(matrix("NA",nrow=length(samplenames), ncol=6))
colnames(meta)=c("ID","Chem", "Dose","CID", "ChemDose", "dup")
meta$ID=samplenames
meta$Chem=Chem
meta$Dose=Dose
meta$ChemDose=ChemDose
meta$dup=SeqID
meta$CID=CID

###########################################
library("DESeq2")
y=round(cnts)
ddsMat <- DESeqDataSetFromMatrix(countData = y,colData = meta, design = ~ ChemDose)

#nrow(ddsMat)###14408
#keep <- rowSums(counts(dds)) > 1
# at least 3 samples with a count of 10 or higher

keep <- rowSums(counts(ddsMat) >= 10) >= 3
dds <- ddsMat[keep,]

###nrow(dds)###11287
#install.packages("hexbin")
###################################################################
#DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function, and the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014).For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. The VST or rlog-transformed data then become approximately homoskedastic (more flat trend in the meanSdPlot), and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and compare the meanSdPlot or PCA plots generated, as described below.
##########################################################
## vst: variance stablizing transformation----------------------------
###In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
##############################################
library("vsn")
vsd <- vst(dds, blind = FALSE)
#head(assay(vsd), 3)
#colData(vsd)
############################################################################
## ----rlog is extremely time consuming, consider to get away from it. 
##########################################################################
#rld <- rlog(dds, blind = FALSE)
#head(assay(rld), 3)
######################################################################
## ----transformplot, fig.width = 6, fig.height = 2.5------------------------
############################################################
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)
new_cnt<-counts(dds, normalized=TRUE)
###################################
#These plot is not as good as the limma data
#I do not like the new count glimma plot, did not seperate the sample well. 
#glMDSPlot(new_cnt, labels=group, 
#          groups=trt, launch=TRUE)

library(Glimma)
glMDSPlot(log2(cnts+1), labels=group, 
          groups=ChemDose, launch=TRUE)
##############################################
write.csv(new_cnt, file = "Deseq.normalized.count.csv")
#write.csv(vsd, file = "Deseq.vsd.csv")
#write.csv(rld, file = "Deseq.rld.csv")
write.csv(cnts, file = "Desq.cleaned.count.csv")

####################################
#df <- bind_rows(
#  as_data_frame(log2(new_cnt[, 1:2]+1)) %>%
#         mutate(transformation = "log2(x + 1)"),
#  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

#colnames(df)[1:2] <- c("x", "y")  
#jpeg('vst-rlog-log2.jpg')
#ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
#  coord_fixed() + facet_grid( . ~ transformation)  
#dev.off()
###########################################################################
## https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
#install.packages("pheatmap")
########################################################
sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")

## ----distheatmap, fig.width = 6.1, fig.height = 4.5------------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
my_sample_col <- data.frame(sample = group)
row.names(my_sample_col) <- colnames(vsd)
jpeg('sampleDistance.jpg', width = 3000, height = 3000,  res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
	 #annotation_row = my_sample_col,
         col = colors)
dev.off()
#######################################################
## Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the counts in dds.

#install.packages("PoiClaClu")
###################################################
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
## ----poisdistheatmap, fig.width = 6.1, fig.height = 4.5--------------------
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- colnames(dds)
colnames(samplePoisDistMatrix) <- NULL

jpeg('samplePoisonDistance.jpg', width = 3000, height = 3000,  res = 300)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,	 
         col = colors)
dev.off()
###################################################################
## ----plotpca, fig.width=6, fig.height=4.5----------------------------------
###############################################################
jpeg('PCA.jpg',width=885, height=600)
plotPCA(vsd, intgroup = "ChemDose")
dev.off()
## --------------------------------------------------------------------------
pcaData <- plotPCA(vsd, intgroup = "ChemDose", returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
## ----ggplotpca, fig.width=6, fig.height=4.5--------------------------------
library(ggpubr)
jpeg('PCA2.jpg',width=1285, height=900)
ggplot(pcaData, aes(x = PC1, y = PC2, color = ChemDose)) +
  geom_text(size =8, label=SeqID) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+theme(legend.text = element_text(size = 18))+
  ggtitle("PCA with VST data")
dev.off()

## ----ggplotpca, fig.width=6, fig.height=4.5--------------------------------
library(ggpubr)
jpeg('PCA3.jpg',width=1285, height=900)
ggplot(pcaData, aes(x = PC1, y = PC2, color = ChemDose)) +
geom_text(size =8, label=SeqID) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+
  stat_chull(aes(color=ChemDose, fill=ChemDose), alpha=0.1, geom="polygon")+
  ggtitle("PCA with VST data")+theme(legend.text = element_text(size = 18))
dev.off()
########################Daw PC1 Data ##########################
# Sort the data frame by PC1 in ascending order
sorted_df <- pcaData[order(pcaData$PC1), ]

# Create a bar plot of PC1
jpeg('PC1_BarPlot.jpg', width = 385, height = 700)

library(ggplot2)

ggplot(sorted_df, aes(y = reorder(row.names(sorted_df), PC1), x = PC1)) +
  geom_bar(stat = "identity", fill = "aquamarine") +
  ylab("") +
  xlab("PC1") +
  ggtitle("PC1 Bar Plot (Sorted)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()


sorted_df2 <- pcaData[order(pcaData$PC2), ]

# Create a bar plot of PC1
jpeg('PC2_BarPlot.jpg', width = 385, height = 700)

library(ggplot2)

ggplot(sorted_df2, aes(y = reorder(row.names(sorted_df2), PC2), x = PC2)) +
  geom_bar(stat = "identity", fill = "aquamarine") +
  ylab("") +
  xlab("PC2") +
  ggtitle("PC2 Bar Plot (Sorted)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()



################################################
##MDS plot
############################################
mds <- as.data.frame(colData(vsd))  %>%
         cbind(cmdscale(sampleDistMatrix))
jpeg('MDS.jpg',width=1285, height=900)
ggplot(mds, aes(x = `1`, y = `2`, color = ChemDose)) +theme_bw()+
geom_text(size = 8, label=SeqID) + coord_fixed()+ ggtitle("MDS with VST data")+theme(legend.text = element_text(size = 18))
dev.off()
####################################################
mdsPois <- as.data.frame(colData(dds)) %>%
   cbind(cmdscale(samplePoisDistMatrix))
jpeg('MDS2.jpg',width=1285, height=900)
ggplot(mdsPois, aes(x = `1`, y = `2`, color = ChemDose)) +theme_bw()+
geom_text(size = 8, label=SeqID) + coord_fixed()+ ggtitle("MDS with PoissonDistances")+theme(legend.text = element_text(size = 18))
dev.off()


########################Analysis##################################
results <- DESeq(dds)
for (i in unique(ChemDose)){
    if(i != "DMSO0_1Percent"){
    	 comparison1 <- results(results, pAdjustMethod = "BH", contrast = c("ChemDose",i,"DMSO0_1Percent"))
	 comparison1a <- cbind(rownames(dds),comparison1)
	 colnames(comparison1a)[1] <- "Genes"
	 XX=2**comparison1a$log2FoldChange
	 comparison1a$TrueFC <- sapply(XX,function(x){ifelse(x>=1, x, -1/x)})
	 YY=cbind(rownames(dds), comparison1a$padj, comparison1a$TrueFC)
	 write.table(YY, file = paste0(i,"_vs_Ctl.stat"), sep="\t")
	 write.table(comparison1a, file = paste0(i, "_vs_Ctl.txt"), sep="\t")
	 print (paste0(i, "_vs_Ctl", " p<=0.05 ",  sum(comparison1a$pvalue<=0.05, na.rm=TRUE), " padj<=0.05 ", sum(comparison1a$padj<=0.05, na.rm=TRUE)))
    }
}

sessionInfo()