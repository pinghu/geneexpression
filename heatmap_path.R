remove(list=ls())
library(ComplexHeatmap)
filename="path_ko00540.xls"
outname="P00540_LPS"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
lfc <- A[, 2:17]
#rownames(lfc) <- paste0(A[, 68],".", A[,71], ".",A[,1])
rownames(lfc) <- paste0( A[,71], ".",A[,1])
# Remove ".lfc" from column names
colnames(lfc) <- sub("\\.lfc$", "", colnames(lfc))
p <- A[,19:34]
adjp<-A[,36:51]
rownames(p)<- rownames(lfc)
rownames(adjp) <- rownames(lfc)
colnames(p) <-sub("\\.p$", "", colnames(p))
colnames(adjp) <-sub("\\.adjp$", "", colnames(adjp))
# Create an empty matrix with the same dimensions as p
p_mat <- matrix("", nrow = nrow(p), ncol = ncol(p))

# Set the row and column names
rownames(p_mat) <- rownames(p)
colnames(p_mat) <- colnames(p)
p_mat[!is.na(p) & p <= 0.05] <- "*"
p_mat[!is.na(adjp) & adjp <= 0.05] <- "**"

# Create row annotation with row names
data <- as.matrix(lfc)
row_avg <- rowMeans(data)
col_avg <- colMeans(data)
#data <- pmin(pmax(data, -3), 3)
library(circlize)
#color_scale = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
#color_scale = colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"))
color_scale = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
#color_scale <- colorRamp2(c(max(data),0, min(data)), c("red","white", "blue"))
row_ha <- rowAnnotation(LFC = anno_barplot(row_avg))
column_ha <- HeatmapAnnotation( LFC= anno_barplot(col_avg))
bacteria <-as.factor(A[, 68])
ht_list<-Heatmap(data, name = "LogFoldChange", top_annotation = column_ha, 
                 left_annotation = row_ha,  row_names_gp = gpar(fontsize = 7),
                 cluster_rows = FALSE, cluster_columns  = FALSE,
                 #rect_gp = gpar(col = "white", lwd = 2),
                 row_split =bacteria,
                 cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                   grid.text(p_mat[i, j], x, y)
                 },
                 col = color_scale)
jpeg("heatmap.jpg", width = 6, height = 8, units = "in", res = 300)
draw(ht_list)#, heatmap_legend_side = "top")
dev.off()

##Other choice
#column_title_gp = gpar(fontsize = 20, fontface = "bold")
#column_title = "I am a column title", column_title_side = "bottom",
#row_title = "I am a row title", show_column_dend = FALSE, row_dend_side = "right", column_dend_side = "bottom",column_dend_height = unit(4, "cm"), 
#row_dend_width = unit(4, "cm"), clustering_distance_rows = "pearson",clustering_distance_rows = function(m) dist(m),clustering_distance_rows = function(x, y) 1 - cor(x, y),
# robust_dist = function(x, y) {
#   qx = quantile(x, c(0.1, 0.9))
#   qy = quantile(y, c(0.1, 0.9))
#   l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
#   x = x[l]
#   y = y[l]
#   sqrt(sum((x - y)^2))
# }
# clustering_distance_rows = robust_dist,
# clustering_distance_columns = robust_dist,
#group = kmeans(t(mat), centers = 3)$cluster
#cluster_columns = cluster_within_group(mat, group)
#cluster_within_group(data, bacteria),
#row_names_gp = gpar(col = )
#row_names_centered = TRUE, column_names_centered = TRUE
#column_names_rot = 45,
# row_names_max_width = max_text_width(
#   rownames(mat2), 
#   gp = gpar(fontsize = 12)
# )
# row_labels = row_labels[rownames(mat)], 
# column_labels = column_labels[colnames(mat)]
# row_labels = expression(alpha, beta, gamma, 
#                         delta, epsilon, zeta, eta, theta, iota, kappa, lambda, mu, nu, xi, 
#                         omicron, pi, rho, sigma)
# row_labels = gt_render(letters[1:18], padding = unit(c(2, 10, 2, 10), "pt")),
# row_names_gp = gpar(box_col = rep(c("red", "green"), times = 9))