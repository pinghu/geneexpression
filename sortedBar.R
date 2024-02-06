rm(list=ls())
library(ggplot2)
library(fs)
# Read the data matrix from a tab-delimited file
args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
filename="Irritation3Set1.txt"
A <- read.table(filename, header = TRUE,  sep="\t")
d=dim(A)
data=A[1:d[1], 6:d[2]]

# Get the number of rows and columns in the data matrix
num_rows <- nrow(data)
num_cols <- ncol(data)
# Create a directory for the image files
dir_path <- "./figures"
dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

# Create a barplot for each row
for (i in 1:num_rows) {
  
  fc <- unlist(data[i, ])
  genesymbol=A[i,1]
  genename <- sub("\\[.*", "", A[i,4])  # Remove anything after "["
  # Define colors based on conditions
  color_vector <- ifelse(grepl("etro",colnames(data)), "blue",
                         ifelse(grepl("SLS|TNF|Geraniol",colnames(data) ), "red", "steelblue"))
  testData <- data.frame(Sample = colnames(data), FoldChange = fc, Col=color_vector )
  sorted_df <- testData[order(testData$FoldChange), ]
 

  p<-ggplot(sorted_df, aes(y = reorder(row.names(sorted_df), FoldChange), x = FoldChange)) +
    geom_bar(stat = "identity",  fill = sorted_df$Col) +
    ylab(genename) +
    xlab("") +
    ggtitle(paste(genesymbol, " fold change")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 12))

  # Save the plot as a JPEG image file in the subdirectory
    filename <- paste0(dir_path,"/", genesymbol, ".jpg")
    jpeg(filename, width = 285, height = 600)
    print(p)
    dev.off()
  
}

sessionInfo()