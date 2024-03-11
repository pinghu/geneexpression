rm(list=ls())
library("DESeq2")
library("dplyr")
library(geosphere)
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(rstatix)
library("vsn")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <-args[2]
#filename="tmp2.Porphyromonas_gingivalis_ATCC33277.ann.xls"
#outname="test"

drawSortedBarPlot <- function(testData, filename) {
  # Ensure testData is a data frame
  if (!is.data.frame(testData)) {
    testData <- as.data.frame(testData)
  }
  
  # Make sure testData has the expected columns
  if(!("Sample" %in% names(testData)) || !("similarityScore" %in% names(testData))) {
    stop("testData must have 'Sample' and 'similarityScore' columns")
  }
  
  # Convert similarityScore to numeric if it's not already
  testData$similarityScore <- as.numeric(as.character(testData$similarityScore))
  
  # In case of any NA introduced by conversion (e.g., non-numeric values present), handle or warn
  if(any(is.na(testData$similarityScore))) {
    warning("NAs introduced by coercion to numeric in similarityScore")
  }
  
  color_vector<- ifelse(grepl("Ctl",testData$Sample), "blue",
                        ifelse(grepl("SnF|SnCl",testData$Sample), "red", "steelblue"))
  testData$Col =color_vector
  
  sorted_df <- testData[order(testData$similarityScore), ]
  p <- ggplot(sorted_df, aes(y = reorder(Sample, similarityScore), x = similarityScore)) +
    geom_bar(stat = "identity", fill = sorted_df$Col) +
    ylab("Rank") +
    xlab("Distance Score") +
    ggtitle(paste(filename)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 12))
  
  # Save the plot as a JPEG image file
  filename_with_extension <- paste0(filename, ".rank.jpg")
  jpeg(filename_with_extension, width = 285, height = 380)
  print(p)
  dev.off()
  
  # Use rank with ties.method = "min" to handle scores that are the same
  sorted_df$order <- rank(sorted_df$similarityScore, ties.method = "min")
  
  # Specify your desired output file name
  output_tab_delimited_file <- paste0(filename, ".rank")
  
  # Write the specific columns to a tab-delimited text file
  write.table(sorted_df[, c("Sample", "order", "similarityScore")], 
              file = output_tab_delimited_file, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE)
}

calculateSimilarityScores <- function(distance_dict) {
  # Flatten the distance dictionary to a numeric vector
  myDist <- unlist(distance_dict)
  
  # Find minimum and maximum distances
  min_dist <- min(myDist)
  max_dist <- max(myDist)
  
  # Calculate similarity scores based on distances from -1 to 1
  #similarityScore <- 2 * (myDist - min_dist) / (max_dist - min_dist) - 1
  # Calculate similarity scores based on distances from 0 to 1
  similarityScore <- (myDist - min_dist) / (max_dist - min_dist)
  return(similarityScore)
}

cnts=read.csv(filename,sep='\t')
d=dim(cnts)
geneID=cnts[,1]
samplenames <- colnames(cnts)[2:d[2]]
splitname<-strsplit(samplenames, "[.]")
Clen=length(samplenames)
trt=rep("NA", Clen)
id=rep("NA", Clen)
dup=rep("NA", Clen)
for(mm in  1:Clen ){
  trt[mm]=splitname[[mm]][1]
  id[mm]=splitname[[mm]][3]
  dup[mm]=splitname[[mm]][2]
}
group <- as.factor(trt)
samplenames=paste0(trt, ".", id)
cnts = cnts[,2:d[2]]
colnames(cnts)=samplenames
rownames(cnts)=geneID
meta<- data.frame(matrix("NA",nrow=length(samplenames), ncol=4))
colnames(meta)=c("ID","trt", "dup","id")
meta$ID=samplenames
meta$trt=trt
meta$id=id
meta$dup=dup
###########################################
library("DESeq2")
y=round(cnts)
ddsMat <- DESeqDataSetFromMatrix(countData = y,colData = meta, design = ~ trt)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
dds <- ddsMat[keep,]

###################################################################
#DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function, and the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014).For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. The VST or rlog-transformed data then become approximately homoskedastic (more flat trend in the meanSdPlot), and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and compare the meanSdPlot or PCA plots generated, as described below.
##########################################################
## vst: variance stablizing transformation----------------------------
###In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
##############################################

vsd <- vst(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)
new_cnt<-counts(dds, normalized=TRUE)

###################################################################
## ----plotpca, fig.width=6, fig.height=4.5----------------------------------
###############################################################



pcaData <- plotPCA(vsd, intgroup = "trt", returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))



# Calculate the reference point
reference_point <- pcaData %>%
  filter(trt == "Ctl") %>%  # Filter for control group
  summarise(
    avg_PC1 = mean(PC1),  # Calculate average of PC1
    avg_PC2 = mean(PC2)   # Calculate average of PC2
  )

# Calculate the average PC1 and PC2 for each treatment group
avePCAData <- pcaData %>%
  group_by(trt) %>%
  summarise(
    avePC1 = mean(PC1, na.rm = TRUE), # Calculate the average of PC1, remove NA values
    avePC2 = mean(PC2, na.rm = TRUE)  # Calculate the average of PC2, remove NA values
  )

jpeg(paste0(outname, '.PCA2.jpg'), width=600, height=380)
ggplot(pcaData, aes(x = PC1, y = PC2, color = trt)) +
geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+
  stat_chull(aes(color=trt, fill=trt), alpha=0.1, geom="polygon") +
  geom_text(data = avePCAData, aes(label = trt, x = avePC1, y = avePC2), nudge_x =2, nudge_y = 2) + # Label first points
  ggtitle(paste0 (outname, " PCA with VST data"))
dev.off()

jpeg(paste0(outname, '.PCA0.jpg'), width=600, height=380)
ggplot(pcaData, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+
  stat_chull(aes(color=trt, fill=trt), alpha=0.1, geom="polygon") +
  ggtitle(paste(outname, "PCA with VST data"))
dev.off()

distance_dict <- list()
for (i in 1:dim(avePCAData)[1]) {
  myDist<-distm (avePCAData[i,2:3], reference_point, fun = distGeo)
  #myDist <- rdist(as.matrix(result$rotation[i,1:2, drop = FALSE]), as.matrix(reference_point))
  myName=as.character(avePCAData$trt[i])
  distance_dict[myName] <- myDist
}
pca_similarityScore=calculateSimilarityScores(distance_dict)
myResult=cbind(names(distance_dict),unlist(distance_dict), pca_similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
write.table(myResult, file = paste0(outname, ".diffRank.xls"), sep = "\t", row.names = FALSE)
drawSortedBarPlot(myResult, outname)

sessionInfo()