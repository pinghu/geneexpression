rm(list=ls())
library(geosphere)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
target_row_name <- args[2]
rm(args)
#file <- "petro.count1.5.xingtaoLFCDirection"
#target_row_name <- "direction.LFC_A1_V1"
outname=paste0(file, ".", target_row_name, ".pca_rank")
#file <- "test_data"
#target_row_name <- "LFC_T_3_vs_V1"
A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
heatmap_data <- A[, 3:d[2]]
row_names <- A[, 1]
row.names(heatmap_data) <- A[, 1]
data_matrix <- as.matrix(heatmap_data)
result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
row_index <- which(colnames(data_matrix) == target_row_name)
reference_point <- result$rotation[row_index,1:2]

distance_dict <- list()
for (i in 1:dim(data_matrix)[2]) {
  myDist<-distm (result$rotation[i,1:2], reference_point, fun = distGeo)
  myName=colnames(data_matrix)[i]
  #print (paste(myName, myDist))
  distance_dict[[myName]] <- myDist
}
myDist=unlist(distance_dict)
min_dist <- min(myDist)
max_dist <- max(myDist)
#sorted_distance_dict <- distance_dict[order(unlist(distance_dict))]
#print(sorted_distance_dict)
#####ranking score from -1 to 1, 1 means exactly same, or shortest distant##########
similarityScore=2*(myDist-min_dist)/(min_dist-max_dist) +1
myResult=cbind(colnames(data_matrix),myDist, similarityScore )
write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
