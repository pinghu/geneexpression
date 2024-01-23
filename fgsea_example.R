

# Load libraries ----------------------------------------------------------

#load libraries
library(fgsea)
library(data.table)
library(ggplot2)

# Load data ---------------------------------------------------------------

#load expression data
expression.file <- #path to data
y <- read.delim(expression.file, comment.char = "#",row.names = 1)

#pass a single vector of probenames for gsea (e.g., sample expression, ranked genes, etc.)
#instead of expression, this could be differential expression, correlation, etc.
input.vector <- y[,1] #this is the expression of sample number 1

#path to save results table
save.file <- #myresults.txt
  
# load pathways -----------------------------------------------------------

#Use gmt_files downloaded from MSigDb 
gmt.file <- #path to "go.bp.c5.bp.v7.0.symbols.gmt"
pathways <- gmtPathways(gmt.file)

#Or make your own custom list
#custom.pathways <- list(MyPathway1 = c("Gene1", "Gene2",...),
#                        MyPathway2 = c("Gene3", "Gene4",...))
#pathways <- custom.pathways

# Custom function to prepare for gsea -------------------------------------

#Get gene list for GSEA from diff exp data
prep.for.gsea <- function(tmp.df, method = "max", rem.dup = T) {
  
  organism = "org.Hs.eg.db"
  require(organism, character.only = TRUE)
  require(hgu133plus2.db)
  require(hgu219.db)
  
  #Detect ids
  check.133 <- intersect(rownames(tmp.df),keys(hgu133plus2.db,keytype="PROBEID"))
  check.219 <- intersect(rownames(tmp.df),keys(hgu219.db,keytype="PROBEID"))
  
  #convert ids to gene symbol
  if (length(check.133) > length(check.219)) {
    id.conversion <- AnnotationDbi::select(hgu133plus2.db, keys=check.133, 
                                           columns=c("SYMBOL"), keytype="PROBEID")
  } 
  
  if (length(check.133) < length(check.219)) {
    id.conversion <- AnnotationDbi::select(hgu219.db, keys=check.219, 
                                           columns=c("SYMBOL"), keytype="PROBEID")
  } 
  
  #Index conversions table with unique ids
  ix.id.conversion <- id.conversion
  current.names <- ix.id.conversion$PROBEID
  tmp.names <- strsplit(make.names(ix.id.conversion$PROBEID,unique = T),split="X")
  converted.names <- c()
  for (i in 1:length(tmp.names)) {
    if(!length(tmp.names[i][[1]])>2) {
      converted.names <- c(converted.names,tmp.names[i][[1]][2])
    } else converted.names <- c(converted.names,current.names[i])
  }
  
  #Convert rownames
  rownames(ix.id.conversion) <- converted.names
  
  #return if not removing duplicate gene names
  if(!rem.dup == T) {
    
    #expand list
    avg.logFC <- as.vector(tmp.df[,1])
    gnames <- as.vector(ix.id.conversion[rownames(tmp.df),"SYMBOL"])
    
    #Format gene list
    gene_list <- avg.logFC
    names(gene_list) <- gnames
    gene_list <- gene_list[order(gene_list,decreasing = T)]
    
    return(gene_list)
  }
  
  #Add gene symbol identifier
  tmp.df.ix <- tmp.df
  tmp.df.ix$SYMBOL <- ix.id.conversion[rownames(tmp.df),"SYMBOL"]
  
  #Tabulate gene name mapping
  tmp.tbl <- table(tmp.df.ix$SYMBOL)
  tmp.sngl <- names(tmp.tbl[tmp.tbl==1])
  tmp.dups <- names(tmp.tbl[tmp.tbl>1])
  
  #subset 1:1 map
  tmp.df.sngl <- tmp.df.ix[tmp.df.ix$SYMBOL %in% tmp.sngl,]
  
  #Get logFC values for the multiple maps
  logfc <- c()
  for (gene in tmp.dups) {
    
    tmp.subset <- tmp.df.ix[which(tmp.df.ix$SYMBOL==gene),]
    
    if(method == "mean") {
      mean.logFC <- mean(as.vector(tmp.subset[,1]))
      logfc <- c(logfc,mean.logFC)
    }
    
    if(!method == "mean") {
      
      arg.min <- which.max(abs(tmp.subset[,1]))
      logfc <- c(logfc,tmp.subset[arg.min,1])
    }
  }
  
  #expand list
  avg.logFC <- c(as.vector(tmp.df.sngl[,1]),logfc)
  gnames <- c(tmp.df.sngl$SYMBOL,tmp.dups)
  
  #Format gene list
  gene_list <- avg.logFC
  names(gene_list) <- gnames
  gene_list <- gene_list[order(gene_list,decreasing = T)]
  
  return(gene_list)
  
}

# Perform GSEA ------------------------------------------------------------

#format genes for gsea
tmp.gene.list <- prep.for.gsea(as.data.frame(input.vector))

#perform fgsea
tmp.gsea <- fgsea(pathways = pathways, 
                  stats    = tmp.gene.list,
                  minSize  = 4,
                  maxSize  = 800,
                  nproc=4)

#format results as data frame
tmp.df <- data.frame(NES = tmp.gsea$NES,pval = tmp.gsea$pval, padj = tmp.gsea$padj)
rownames(tmp.df) <- tmp.gsea$pathway

#add leading edge genes
le.genes <- c()
for (i in 1:dim(tmp.df)[1]) {
  tmp.le <- paste(tmp.gsea$leadingEdge[[i]],collapse="/")
  le.genes <- c(le.genes,tmp.le)
}
tmp.df$leadingEdge <- le.genes

#save results to file
write.table(tmp.df, file = save.file, append = FALSE, sep = "\t",
            row.names = F)
