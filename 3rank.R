rm(list=ls())
#install.packages('bigmemory')
#BiocManager::install("NMF", dependencies = TRUE)
library(bigmemory)
library(NMF)
library(geosphere)
library(ggplot2)
library(Rtsne)
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
  
  sorted_df <- testData[order(testData$similarityScore),]
  
  p <- ggplot(sorted_df, aes(y = reorder(Sample, similarityScore), x = similarityScore)) +
    geom_bar(stat = "identity", fill = "blue") +
    ylab(filename) +
    xlab("Similarity Score") +
    ggtitle("Rank") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 12))
  
  # Save the plot as a JPEG image file
  filename_with_extension <- paste0(filename, ".jpg")
  jpeg(filename_with_extension, width = 285, height = 600)
  print(p)
  dev.off()
}

calculateSimilarityScores <- function(distance_dict) {
  # Flatten the distance dictionary to a numeric vector
  myDist <- unlist(distance_dict)
  
  # Find minimum and maximum distances
  min_dist <- min(myDist)
  max_dist <- max(myDist)
  
  # Calculate similarity scores based on distances
  similarityScore <- 2 * (myDist - min_dist) / (max_dist - min_dist) - 1
  
  return(similarityScore)
}



args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
target_row_name <- args[2]
rm(args)
file <- "petro.count1.5.xingtaoLFCDirection"
target_row_name <- "direction.LFC_A1_V1"
#file <- "test_data"
#target_row_name <- "LFC_T_3_vs_V1"
outname=paste0(file, ".", target_row_name, ".nmf_rank")


A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
heatmap_data <- A[, 3:d[2]]
row_names <- A[, 1]
row.names(heatmap_data) <- A[, 1]

######## Part 1 NMF Ranking #########
outname=paste0(file, ".", target_row_name, ".nmf_rank")
ZZ=as.numeric(min(heatmap_data))
data_matrix <- as.matrix(t(heatmap_data)) -ZZ+1
result <- nmf(data_matrix, rank = 2) # Change rank as needed
W_matrix <- basis(result)
row_index <- which(rownames(data_matrix) == target_row_name)
reference_point <- W_matrix[row_index, ]
distance_dict <- list()
for (i in 1:dim(data_matrix)[1]) {
  myDist<-distm (W_matrix[i,], W_matrix[row_index,], fun = distGeo)
  myName=rownames(data_matrix)[i]
  distance_dict[[myName]] <- myDist
}
nmf_similarityScore=calculateSimilarityScores(distance_dict )
myResult=cbind(rownames(data_matrix), unlist(distance_dict), nmf_similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
#write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
drawSortedBarPlot(myResult, paste0(outname))

#######Part2 PCA Ranking ##########################
outname=paste0(file, ".", target_row_name, ".pca_rank")
data_matrix <- as.matrix(heatmap_data)
result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
row_index <- which(colnames(data_matrix) == target_row_name)
reference_point <- result$rotation[row_index,1:2]
distance_dict <- list()
for (i in 1:dim(data_matrix)[2]) {
  myDist<-distm (result$rotation[i,1:2], reference_point, fun = distGeo)
  myName=colnames(data_matrix)[i]
  distance_dict[[myName]] <- myDist
}
pca_similarityScore=calculateSimilarityScores(distance_dict)
myResult=cbind(colnames(data_matrix),unlist(distance_dict), pca_similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
#write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
drawSortedBarPlot(myResult, paste0(outname))

################Part 3 TSNE Ranking ##################################
outname=paste0(file, ".", target_row_name, ".tsne_rank")
data_matrix <- as.matrix(t(heatmap_data))
result <- Rtsne(data_matrix, dims = 2, perplexity = 5)
row_index <- which(rownames(data_matrix) == target_row_name)
reference_point <- result$Y[row_index, ]
distance_dict <- list()
for (i in 1:dim(data_matrix)[1]) {
  myDist<-distm (result$Y[i,], result$Y[row_index,], fun = distGeo)
  myName=rownames(data_matrix)[i]
  distance_dict[[myName]] <- myDist
}
tsne_similarityScore=calculateSimilarityScores(distance_dict)
myResult=cbind(rownames(data_matrix),unlist(distance_dict), tsne_similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
#write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
drawSortedBarPlot(myResult, paste0(outname))
###########################################Combine average the ranking score#####################
scores_matrix <- cbind(pca_similarityScore, nmf_similarityScore, tsne_similarityScore)
scores_mean <- rowMeans(scores_matrix, na.rm = TRUE)  # na.rm=TRUE to remove any NA values in the calculation
combineScore <- cbind(rownames(data_matrix), pca_similarityScore, tsne_similarityScore, nmf_similarityScore, Mean=scores_mean)
colnames(combineScore)=c("Sample", "pca_similarityScore", "tsne_similarityScore", "nmf_similarityScore", "similarityScore")
outname=paste0(file, ".", target_row_name, ".combined_rank")
write.table(combineScore, file = outname, sep = "\t", row.names = FALSE)
drawSortedBarPlot(combineScore, paste0(outname))