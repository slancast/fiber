#This project will perform supervised clustering
#Do pcl data too
#how to plot this with machine learning data? Barplot with stars?
library(supclust)

binary_df <- read.table(paste("/Volumes/My_Book_Duo/B-Cell_Epitope_Prediction_GDrive/Metabolomic-PCL/pcl_LCInulin_UPDATED.csv",sep=""),sep=",",header=TRUE)
binary_df <- read.table(paste("/Users/SLancaster/Desktop/Projects/Fiber/Supervisied_Clustering/metabolomics_LCinulin.csv",sep=""),sep=",",header=TRUE)
binary_df <- data.frame(t(binary_df))
binary_df <- na.omit(binary_df)
binary_df <- data.frame(t(binary_df))
binary_df <- binary_df[order(binary_df$Timepoint),]
binary_df <- as.matrix(binary_df)
class(binary_df) <- "numeric"
df <- binary_df[,2:ncol(binary_df)]
df[df==0.0000000000000355]<-0 #Getting rid of the log2(0.0000001)
df <- t(df)
df <- df[rowSums(df) != 0, ] ####### > vs !=
df <- t(df)
ncol(df)

y <- binary_df[,1]
rownames(df) <- y

fit <- pelora(df, y, noc = 3)

summary(fit)
plot(fit)

genes1.loc <- fit$genes[[1]]
genes1 <- fit$gene.names[genes1.loc]

genes2.loc <- fit$genes[[2]]
genes2 <- fit$gene.names[genes2.loc]

genes3.loc <- fit$genes[[3]]
genes3 <- fit$gene.names[genes3.loc]

genes_total <- unique(c(genes1, genes2, genes3))

sig_genes_heatmap <- df[,genes_total]

h <- heatmap(sig_genes_heatmap, cexCol = 0.2)
# library(gplots)
# heatmap.2(sig_genes_heatmap)
# heatmap.2(heatmap_matrix, dendrogram="row", Colv = FALSE, margins = c(5,20), cexRow=.5)

colnames(sig_genes_heatmap)[h$colInd]


if (FALSE) {
#To make the barplot of the machine learning predictors and compare it to the supervised
#clustering for LCInulin
random_forest_stats <- read.table(paste("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/metabolomics_stats.csv",sep=""),sep=",", header=TRUE, row.names=1)
random_forest_stats <- random_forest_stats[order(random_forest_stats$metabolomics_LCInulin, decreasing = TRUE),]

for_plotting <- sort(random_forest_stats$metabolomics_LCInulin, decreasing=TRUE)[1:15]
overlap <- match(rownames(random_forest_stats)[1:15], colnames(sig_genes_heatmap)[h$colInd])

names(for_plotting) <- rownames(random_forest_stats)[1:15]
par(mar=c(15,4,1,1))
b <- barplot(for_plotting,las=2)
text(b,for_plotting,ifelse(is.na(overlap),"","*"),pos=3,cex=2,xpd=NA)
}

if (FALSE) {
  #To make the bar plot of the machine learning predictors and compare it to the supervised clustering
  #For arabinoxylan
  random_forest_stats <- read.table(paste("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/metabolomics_stats.csv",sep=""),sep=",", header=TRUE, row.names=1)
  random_forest_stats <- random_forest_stats[order(random_forest_stats$metabolomics_arabinoxylan, decreasing = TRUE),]
  
  for_plotting <- sort(random_forest_stats$metabolomics_arabinoxylan, decreasing=TRUE)[1:15]
  overlap <- match(rownames(random_forest_stats)[1:15], colnames(sig_genes_heatmap)[h$colInd])
  
  names(for_plotting) <- rownames(random_forest_stats)[1:15]
  par(mar=c(15,4,1,1))
  b <- barplot(for_plotting,las=2)
  text(b,for_plotting,ifelse(is.na(overlap),"","*"),pos=3,cex=2,xpd=NA)
}

if (FALSE) {
  #The pelora (used above) produces a much better fit than the wilma
  #Saving the code here in case it is useful in the future
  fit <- wilma(df, y, noc = 3, trace = 1)
  
  summary(fit)
  plot(fit)
  fitted(fit)
  
  genes1 <- fit$clist[[1]]
  genes2 <- fit$clist[[2]]
  genes3 <- fit$clist[[3]]
  genes_total <- unique(c(genes1, genes2, genes3))
  
  sig_genes_heatmap <- df[,genes_total]
  
  heatmap(sig_genes_heatmap, Rowv=as.numeric(rownames(df)))
}
