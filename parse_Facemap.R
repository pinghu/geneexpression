library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the Excel file path
#excel_path <- args[1]
#output_file_path <- args[2]
#q_filter <- args[3]
excel_path <- "FMap_TNF4071TK_oldCMP-2024_02_14_1723.xlsx"
output_file_path <- "FMap_TNF4071TK_oldCMP"
data <- read_excel(excel_path, sheet = "TKC")
q_filter <- "CMP424"
# Filter rows where column Q equals the q_filter value and select the specified columns
filtered_data <- data %>%
  filter(cell == "tert keratinocytes") %>%
  select(Score, chem, Report.Name, STID,cell, Conc, batches, chips)
# Write the filtered data to a tab-delimited text file
#write.table(filtered_data, paste0(output_file_path, ".all.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.xlsx(filtered_data, paste0(output_file_path, ".xlsx"))
# The third argument is the text filter for column Q
if(q_filter != "NA"){
  small_data <- filtered_data %>%
  filter(batches == q_filter)%>%
  select(Score, Report.Name, Conc)
  
  color_vector <- ifelse(grepl("etro",small_data$Report.Name), "blue",
                         ifelse(grepl("SLS|TNF|Geraniol",small_data$Report.Name ), "red", "steelblue"))
  testData <- data.frame(Sample = paste(small_data$Report.Name, small_data$Conc), FMAP_Score=small_data$Score, Col=color_vector )
  testData$FMAP_Score <- as.numeric(as.character(testData$FMAP_Score))
  sorted_df <- testData[order(testData$FMAP_Score), ]
  
  
  p<-ggplot(sorted_df, aes(y = reorder(Sample, FMAP_Score), x = FMAP_Score)) +
    geom_bar(stat = "identity",  fill = sorted_df$Col) +
    ylab(paste(output_file_path, q_filter) )+
    xlab(paste("FMAP SCORE")) +
    ggtitle("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 12))
  
  # Save the plot as a JPEG image file in the subdirectory
  filename <- paste0(output_file_path,".", q_filter, ".jpg")
  jpeg(filename, width = 285, height = 380)
  print(p)
  dev.off()
  
}



# Notify the user that the file has been written
cat("The filtered data has been written to", output_file_path, "\n")
