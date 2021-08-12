#Supervised clustering on the rna
#
#

library(supclust)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  drops <- c(paste("X053_",fiber_subset,"_Baseline",sep=""),paste("X053_",fiber_subset,"_Baseline",sep=""),paste("X053_",fiber_subset,"_Baseline",sep=""),paste("X053_",fiber_subset,"_Baseline",sep=""))
  set.seed(333)
  binary_df <- read.table(paste("/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/rna/rna_for_supervised_",fiber_subset,"_binary.txt",sep=""),sep="\t",header=TRUE)
  binary_df <- data.frame(t(binary_df))
  binary_df <- na.omit(binary_df)
  binary_df <- data.frame(t(binary_df))
  drops <- c("binary","participant","week")
  df <- binary_df[,!colnames(binary_df) %in% drops]
  
  
  df <- as.matrix(df)
  class(df) <- "numeric"
  #df[df==0.0000000000000355]<-0 #Getting rid of the log2(0.0000001)
  df <- t(df)
  df <- df[rowSums(df) != 0, ] ####### > vs !=
  df <- t(df)
  ncol(df)
  
  y <- binary_df$binary
  y <- as.matrix(y)
  y <- as.numeric(y)
  names(df) <- y
  
  fit <- pelora(df, y, noc = 3, standardize = FALSE)
  
  summary(fit)
  plot(fit)
  
  genes1.loc <- fit$genes[[1]]
  genes1 <- fit$gene.names[genes1.loc]
  
  genes2.loc <- fit$genes[[2]]
  genes2 <- fit$gene.names[genes2.loc]
  
  genes3.loc <- fit$genes[[3]]
  genes3 <- fit$gene.names[genes3.loc]
  
  genes_total <- unique(c(genes1))
  
  sig_genes_heatmap <- df[,genes_total]
  binary_df <- data.frame(binary_df)
  rownames(sig_genes_heatmap) <- paste(binary_df$participant, binary_df$week)
  
  pdf(paste("/Users/SLancaster/Desktop/Supervised_Clustering/Supervised_RF_Figures/supclust/supclust_heatmap_rna_",fiber_subset,".pdf",sep=""))
  par(oma = c(10,0,0,0))
  h <- heatmap(sig_genes_heatmap, cexCol = 0.2, cexRow = 0.35)
  dev.off()
  colnames(sig_genes_heatmap)[h$colInd]
  
  
  #To make the barplot of the machine learning predictors and compare it to the supervised
  random_forest_stats <- read.table(paste("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/rna_stats.csv",sep=""),sep=",", header=TRUE, row.names=1)
  column <- fiber_subset
  random_forest_stats <- random_forest_stats[order(random_forest_stats[,column], decreasing = TRUE),]

  for_plotting <- sort(random_forest_stats[,column], decreasing=TRUE)[1:25]
  overlap <- match(rownames(random_forest_stats)[1:25], colnames(sig_genes_heatmap)[h$colInd])

  pdf(paste("/Users/SLancaster/Desktop/Supervised_Clustering/Supervised_RF_Figures/supclust/supclust_barplot_rna_",fiber_subset,".pdf",sep=""))
  names(for_plotting) <- rownames(random_forest_stats)[1:25]
  par(mar=c(25,4,1,1))
  b <- barplot(for_plotting,las=2)
  text(b,for_plotting,ifelse(is.na(overlap),"","*"),pos=3,cex=2,xpd=NA)
  dev.off()
}

