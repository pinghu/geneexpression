rm(list=ls())
library(openxlsx)

library(org.Hs.eg.db)
#library(DOSE) 
#install.packages("enrichR")
library(enrichR)
#dbs <- listEnrichrDbs()
#write.table(dbs, file = "enrichr.dbs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")

dbs <- c("GO_Biological_Process_2023",
         "GO_Cellular_Component_2023",
         "GO_Molecular_Function_2023",
         "GWAS_Catalog_2023",
         # "Proteomics_Drug_Atlas_2023",
         "WikiPathway_2023_Human",
         "KEGG_2021_Human"#,
         #"GTEx_Tissues_V8_2023", 
         #"Tissue_Protein_Expression_from_Human_Proteome_Map", 
         #"Tabula_Sapiens"
)


######Given a list tell you GO and KEGG matching U133 entrenz gene id
perform_gene_enrichment_analysis <- function(gene_symbols, outname) {
  # Map gene symbols to Entrez Gene IDs using org.Hs.eg.db
  #entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Extract unique Entrez Gene IDs
  #entrez_ids <- unique(entrez_ids$ENTREZID)
  gene_symbols <- unique(gene_symbols)
  
  # Run GO enrichment analysis using clusterProfiler for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)
  #go_BP <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
  #go_MF <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF")
  #go_CC <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC")
  
  # Write the results to files
  # write.table(summary(go_MF), file = paste0(outname, ".goMF"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  # write.table(summary(go_BP), file = paste0(outname, ".goBP"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  # write.table(summary(go_CC), file = paste0(outname, ".goCC"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  # enrichr_result <- enrichr(genes =  gene_symbols , databases  = dbs)
  #  write.table(enrichr_result, file = paste0(outname, ".enrichr"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  # geneU=unique(AsigU[,2])
  combined_results <- data.frame()
  # Initialize an empty data frame to store the filtered results
  combined_results <- data.frame()
  
  for(db in dbs) {
    enrichr_result <- enrichr(genes = gene_symbols, databases = db)
    
    # Access the specific database results and filter for P.value <= 0.05
    # Note: It's assumed enrichr_result is a list with database names as keys
    db_results <- enrichr_result[[db]]  # Access the results for the current database
    filtered_results <- db_results[db_results$P.value <= 0.05, ]
    
    # Add a column to indicate the database for reference in the combined file
    filtered_results$Database <- db
    jpeg(filename = paste0(outname, ".", db, ".jpg"), width = 480, height = 300)  # Adjust size as neede
    p<-plotEnrich(filtered_results, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title= paste(outname, db))
    print(p)
    dev.off()
    # Append the filtered results to the combined_results data frame
    combined_results <- rbind(combined_results, filtered_results)
  }
  
  # Write the combined and filtered results to a single file
  write.table(combined_results, file = paste0(outname, ".enrichr.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  write.xlsx(combined_results, file = paste0(outname, ".enrichr.xlsx"), rowNames = FALSE)
  #plotEnrich(enrichr_result[["GO_Biological_Process_2023"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
outname <- args[2]
rm(args)

#filename="PetroHigh.id"
A<-read.table(filename, sep="\t", header=TRUE)
#perform_gene_enrichment_analysis(A[,1], outname)

Asig=A[as.numeric(A[,6])<=0.05, ]
AsigU = Asig[as.numeric(Asig[,7])>0,]
AsigD = Asig[as.numeric(Asig[,7])<0,]
#5.2 GO classification "MF", "BP", and "CC" subontologies.
#myEZ=AsigU[,3]
perform_gene_enrichment_analysis(AsigU[,2], paste0(outname, ".sigUp"))
perform_gene_enrichment_analysis(AsigD[,2], paste0(outname, ".sigDown"))
perform_gene_enrichment_analysis(Asig[,2], paste0(outname, ".sigAll"))
