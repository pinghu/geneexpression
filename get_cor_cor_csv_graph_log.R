#options(echo=TRUE) # if you want see commands in output file
rm(list=ls())
#args <- commandArgs(trailingOnly = TRUE)
#print(args)


#filename=args[1]
#outname=args[2]
filename ="GOODBadTestData"
outname ="GoodBadCorrelation"
#rm(args)

A<-read.table(filename, sep="\t", header=TRUE)

d <- dim(A);
B <- data.matrix(A[1:d[1],2:d[2]]);
#ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
#test_data=log10(B+ZZ)
test_data=log10(B)
xx<-t(test_data)
Cname=colnames(B)
Clen=length(Cname)
Irritation=rep("NA", Clen)
CMAP_id=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
  Irritation[mm]=splitname[[mm]][1]
  CMAP_id[mm]=splitname[[mm]][2]
}
#B<-cor(xx)
colnames(xx) =A[,1]
rownames(xx)=colnames(A)[2:d[2]]
#https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
library(tidyverse)
hclust_matrix <- xx
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  #t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
png("output_plot.png", width = 800, height = 800) 
plot(gene_hclust)
abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
dev.off() # Close the plotting device
#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
library("ape")
png("output_plot_circular.png", width = 1080, height = 1080, res=120) 
colors = c("red", "blue", "green", "black")
clus4 = cutree(gene_hclust, 4)
plot(as.phylo(gene_hclust), type = "fan", tip.color = colors[clus4],
     label.offset = 1, cex = 0.7)
dev.off()


cutree(gene_hclust, k = 5)
gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
write.table(gene_cluster,file=paste0(outname,'.genecluster.tsv'),sep="\t")
library(ComplexHeatmap)
png("output_heatmap2.png", width = 800, height = 800) 
#Heatmap(hclust_matrix, show_row_names = TRUE)
HM <- Heatmap(hclust_matrix, km = 4, show_row_names = TRUE, show_column_names = T,  cluster_columns = TRUE, show_heatmap_legend = T)
draw(HM, heatmap_legend_side = "right", show_annotation_legend = F)
dev.off()

#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)
res <- cor(xx)
round(res, 1)
#corrplot(res, type = "upper", order = "hclust", 
#        tl.col = "black", tl.srt = 45)
png("output_correlation.png", width = 1800, height = 1800) 
corrplot(res, method = 'number',p.mat = testRes$p, sig.level = 0.05, order = 'hclust', addrect = 4) # colorful number
#corrplot.mixed(res, order = 'AOE')
dev.off()
png("output_correlation3.png", width = 1250, height = 1250)
testRes = cor.mtest(xx, conf.level = 0.95)
## specialized the insignificant value according to the significant level
corrplot(res, p.mat = testRes$p, sig.level = 0.05, order = 'hclust', addrect = 5)
dev.off()



library("Hmisc")
rCOR=rcorr(as.matrix(xx), type="spearman")
write.table(rCOR$P,file=paste0(outname,'.spearman.corP.csv'),sep=",",row.names=A[,1], col.names = A[,1])
write.table(rCOR$r,file=paste0(outname,'.spearman.corR.csv'),sep=",",row.names=A[,1], col.names = A[,1])
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

ff<-flattenCorrMatrix(rCOR$r, rCOR$P)
write.table(ff,file=paste0(outname,'.spearman.corFlat.tsv'),sep="\t")

rCOR2=rcorr(as.matrix(xx), type="pearson")
write.table(rCOR2$P,file=paste0(outname,'.pearson.corP.csv'),sep=",",row.names=A[,1], col.names = A[,1])
write.table(rCOR2$r,file=paste0(outname,'.pearson.corR.csv'),sep=",",row.names=A[,1], col.names = A[,1])

ff2<-flattenCorrMatrix(rCOR2$r, rCOR2$P)
write.table(ff2,file=paste0(outname,'.pearson.corFlat.tsv'),sep="\t")

#tiff(filename = paste0("Fig4a.tiff"), width = 3660, height = 3660, res=300)
#corrplot.mixed(rCOR2$r, p.mat=rCOR2$P, sig.level=0.05, diag='u', insig = "blank",  order="hclust", upper = "ellipse", lower.col = "black", tl.pos="lt")
#dev.off()

############Now check on the feature selection#############
#https://jtr13.github.io/cc21fall2/feature-selection-in-r.html
library(ggcorrplot)
library(caret)
library(randomForest)
library(gam)

#function to find redundant variables
findCorrelation(cor(xx), cutoff=0.75)
#[1] 47 15 44 17 16  8 28  9 39 56 34 13  2 58 49 66 67 24 50 26 60 65 30 42 55 37 40 29 10 33 20  6 51
#[34] 19 23 22 57 14 48 38 53  3  1 31 46 64 54 21 62 63  4 59 25 41 11
#Variable Importance
#use roc_curve area as score
test=t(B)
colnames(test)=A[,1] 
rownames(test)=colnames(B)
test2 <-data.frame(cbind(test,Irritation))
roc_imp <- filterVarImp(x = test2[,1:68], y = as.numeric(as.factor(test2[,69])))

#sort the score in decreasing order
roc_imp <- data.frame(cbind(variable = rownames(roc_imp), score = roc_imp[,1]))
roc_imp$score <- as.double(roc_imp$score)
roc_imp[order(roc_imp$score,decreasing = TRUE),]
#variable        score
#3      cxcl8 2.720551e+16
#36      mmp3 1.124859e+16
#20      mmp1 1.592263e+15
#11    hspa1b 2.158536e+00
#7      cxcl1 9.293928e-01
#44      btg1 4.407797e-16
#1       olr1          NaN
#train random forest model and calculate feature importance
rf = randomForest(x= test2[,1:68],y=as.numeric(as.factor(test2[,69])) )
var_imp <- varImp(rf, scale = FALSE)

#sort the score in decreasing order
var_imp_df <- data.frame(cbind(variable = rownames(var_imp), score = var_imp[,1]))
var_imp_df$score <- as.double(var_imp_df$score)
var_imp_df[order(var_imp_df$score,decreasing = TRUE),]
png("output.featureImportance.png", width = 600, height = 600) 
ggplot(var_imp_df, aes(x=reorder(variable, score), y=score)) + 
  geom_point() +
  geom_segment(aes(x=variable,xend=variable,y=0,yend=score)) +
  ylab("IncNodePurity") +
  xlab("Variable Name") +
  coord_flip()
dev.off()

#Univariate Feature Selection
filterCtrl <- sbfControl(functions = rfSBF, method = "repeatedcv", repeats = 3)
rfWithFilter <- sbf(x= test2[,1:68],y=as.numeric(as.factor(test2[,69])), sbfControl = filterCtrl)
rfWithFilter

#Using the training set, 5 variables were selected:
#  cxcl1, slc7a11, krtap2.3, ccl20, fth1.

#During resampling, the top 5 selected variables (out of a possible 9):
#  ccl20 (76.7%), fth1 (73.3%), krtap2.3 (60%), cxcl1 (56.7%), slc7a11 (53.3%)
#Recursive Feature elimination
filterCtrl <- rfeControl(functions=rfFuncs, method="cv", number=3)
results <- rfe(x= test2[,1:68],y=as.numeric(as.factor(test2[,69])), sizes=c(1:11), rfeControl=filterCtrl)
results
#The top 5 variables (out of 68):
 # krtap2.3, mt1g, psat1, trib3, il13ra2
# plot the results
plot(results, type=c("g", "o"))
results$fit$importance