rm(list=ls())
library(Rtsne)
library(geosphere)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
target_row_name <- args[2]
rm(args)
#file <- "petro.count1.5.xingtaoLFCDirection"
#target_row_name <- "direction.LFC_A1_V1"
outname=paste0(file, ".", target_row_name, ".tsne_rank")
A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
heatmap_data <- A[, 3:d[2]]
row_names <- A[, 1]
row.names(heatmap_data) <- A[, 1]
data_matrix <- as.matrix(t(heatmap_data))

result <- Rtsne(data_matrix, dims = 2, perplexity = 5)

row_index <- which(rownames(data_matrix) == target_row_name)
reference_point <- result$Y[row_index, ]

distance_dict <- list()
for (i in 1:dim(data_matrix)[1]) {
  myDist<-distm (result$Y[i,], result$Y[row_index,], fun = distGeo)
  myName=rownames(data_matrix)[i]
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
myResult=cbind(rownames(data_matrix),myDist, similarityScore )
write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
