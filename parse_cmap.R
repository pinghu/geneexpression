rm(list=ls())
library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the Excel file path
excel_path <- args[1]
output_file_path <- args[2]
q_filter <- args[3] # if it is NA, then no more filter
excel_path <- "CMap_Irritation1_PCCPerfume-2024_02_21_0908.xlsx"
output_file_path <- "CMap_Irritation_PCCPerfume"
q_filter <- "CMP424"

data <- read_excel(excel_path, sheet = "TKC")
filtered_data <- data %>%
  filter(Cell.Line == "tert keratinocytes") %>%
  select(score, chem, Report.Name, Stock.Tube.ID,Cell.Line, Concentration, batch, chip)
# Write the filtered data to a tab-delimited text file
#write.table(filtered_data, paste0(output_file_path, ".all.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.xlsx(filtered_data, paste0(output_file_path, ".xlsx"))

if(q_filter != "NA"){
  small_data <- filtered_data %>%
  filter(batch == q_filter)%>%
  select(chip, score, Report.Name, Concentration)
  
  color_vector <- ifelse(grepl("etro|yristate",small_data$Report.Name), "blue",
                         ifelse(grepl("SLS|TNF|Geraniol",small_data$Report.Name ), "red", "steelblue"))
  testData <- data.frame(Sample = paste(small_data$Report.Name, small_data$Concentration, small_data$chip), CMAP_Score=small_data$score, Col=color_vector )
  testData$CMAP_Score <- as.numeric(as.character(testData$CMAP_Score))
  sorted_df <- testData[order(testData$CMAP_Score), ]
  
  
  p<-ggplot(sorted_df, aes(y = reorder(Sample, CMAP_Score), x = CMAP_Score)) +
    geom_bar(stat = "identity",  fill = sorted_df$Col) +
    ylab(paste(output_file_path, q_filter) )+
    xlab(paste("CMAP SCORE")) +
    ggtitle("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 18))
  
  # Save the plot as a JPEG image file in the subdirectory
  filename <- paste0(output_file_path,".", q_filter, ".jpg")
  jpeg(filename, width = 585, height = 580)
  print(p)
  dev.off()
  
  # Adding an order column to sorted_df for ranking
  #sorted_df$order <- seq_along(sorted_df$CMAP_Score) ### This is not right
  
  # Assuming sorted_df is already sorted by CMAP_Score in descending order
  # Use rank with ties.method = "min" to handle scores that are the same
  sorted_df$order <- rank(-sorted_df$CMAP_Score, ties.method = "min")
  
  # Specify your desired output file name
  output_tab_delimited_file <- paste0(output_file_path, ".", q_filter, ".rank")
  
  # Write the specific columns to a tab-delimited text file
  write.table(sorted_df[, c("Sample", "order", "CMAP_Score")], 
              file = output_tab_delimited_file, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE)
}

