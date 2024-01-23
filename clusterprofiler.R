rm(list=ls())
#remotes::install_github('YuLab-SMU/ggtree') 
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE) 
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
          "KEGG_2021_Human",
        "GTEx_Tissues_V8_2023", 
        "Tissue_Protein_Expression_from_Human_Proteome_Map", 
        "Tabula_Sapiens"
)


######Given a list tell you GO and KEGG matching U133 entrenz gene id
perform_gene_enrichment_analysis <- function(gene_symbols, outname) {
  # Map gene symbols to Entrez Gene IDs using org.Hs.eg.db
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Extract unique Entrez Gene IDs
  entrez_ids <- unique(entrez_ids$ENTREZID)
  gene_symbols <- unique(gene_symbols)
  
  # Run GO enrichment analysis using clusterProfiler for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)
  go_BP <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
  go_MF <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF")
  go_CC <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC")
  
  # Write the results to files
  write.table(summary(go_MF), file = paste0(outname, ".goMF"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  write.table(summary(go_BP), file = paste0(outname, ".goBP"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  write.table(summary(go_CC), file = paste0(outname, ".goCC"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
 # enrichr_result <- enrichr(genes =  gene_symbols , databases  = dbs)
#  write.table(enrichr_result, file = paste0(outname, ".enrichr"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
 # geneU=unique(AsigU[,2])
  for(db in dbs) {
    enrichr_result <- enrichr(genes = gene_symbols  , databases  = db)
    write.table(enrichr_result, file = paste0(outname, ".enrichr.",db), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
  } 
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
#COLN <-as.numeric(args[2])
#chip<-args[3]

rm(args)
#filename="agree" ###remove long gene name, there are issure with it
#filename ="GSS3049.T_V2.TNFalpha.xls"
A<-read.table(filename, sep="\t", header=TRUE)
Asig=A[as.numeric(A[,6])<=0.05, ]
AsigU = Asig[as.numeric(Asig[,7])>0,]
AsigD = Asig[as.numeric(Asig[,7])<0,]
#5.2 GO classification "MF", "BP", and "CC" subontologies.
#myEZ=AsigU[,3]
perform_gene_enrichment_analysis(AsigU[,2], paste0(filename, ".sigUp"))
perform_gene_enrichment_analysis(AsigD[,2], paste0(filename, ".sigDown"))
perform_gene_enrichment_analysis(Asig[,2], paste0(filename, ".sigAll"))


#5.3 GO over-representation test
#Gene ID can be mapped to gene Symbol by using paramter readable=TRUE or setReadable function.

#cnetplot(ego, categorySize="geneNum")
####
#kk <- enrichKEGG(gene= mygene,organism='human',pvalueCutoff = 0.05)
#write.table(summary(kk), file=paste0(filename, ".kk"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
