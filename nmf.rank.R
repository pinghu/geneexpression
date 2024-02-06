rm(list=ls())
#install.packages('bigmemory')
#BiocManager::install("NMF", dependencies = TRUE)
library(bigmemory)
library(NMF)
library(geosphere)
library(ggplot2)
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
    ggtitle(filename) +
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
similarityScore=calculateSimilarityScores(distance_dict )
#tt<-data.matrix(similarityScore)
myResult=cbind(rownames(data_matrix), unlist(distance_dict), similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
#write.table(similarityScore, file="test.txt", sep="\t", row.names=TRUE)


drawSortedBarPlot(myResult, paste0(outname))






