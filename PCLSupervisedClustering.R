# The original supervised clustering file combined both the 
# metabolomics and the pcl into a single file. I am splitting
# them up so that each of them has their own supervised clustering
# file. This will be for the the PCL files

library(supclust)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {

set.seed(333)
binary_df <- read.table(paste("/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/pcl/pcl_for_supervised_",fiber_subset,"_binary.txt",sep=""),sep="\t",header=TRUE)
binary_df <- data.frame(t(binary_df))
binary_df <- na.omit(binary_df)
binary_df <- t(binary_df)
df <- as.matrix(binary_df[,5:ncol(binary_df)])
class(df) <- "numeric"
#df[df==0.0000000000000355]<-0 #Getting rid of the log2(0.0000001)
df <- t(df)
df <- df[rowSums(df) != 0, ] ####### > vs !=
df <- df[, !sapply(df, function(x) { sd(x) == 0} )]
df <- t(df)
ncol(df)

y <- as.numeric(binary_df[,4])
rownames(df) <- y

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
rownames(sig_genes_heatmap) <- paste(binary_df$Participant, binary_df$Dose)

pdf(paste("/Users/SLancaster/Desktop/Supervised_Clustering/Supervised_RF_Figures/supclust_heatmap_pcl_",fiber_subset,".pdf",sep=""))
par(oma = c(10,0,0,0))
h <- heatmap(sig_genes_heatmap, cexCol = 0.2, cexRow = 0.35)
dev.off()
colnames(sig_genes_heatmap)[h$colInd]


  #To make the barplot of the machine learning predictors and compare it to the supervised
  random_forest_stats <- read.table(paste("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/pcl_stats.csv",sep=""),sep=",", header=TRUE, row.names=1)
  column <- paste("pcl_",fiber_subset,"_UPDATED",sep="")
  random_forest_stats <- random_forest_stats[order(random_forest_stats[,column], decreasing = TRUE),]
  
  for_plotting <- sort(random_forest_stats[,column], decreasing=TRUE)[1:15]
  overlap <- match(rownames(random_forest_stats)[1:15], colnames(sig_genes_heatmap)[h$colInd])
  
  pdf(paste("/Users/SLancaster/Desktop/Supervised_Clustering/Supervised_RF_Figures/supclust_barplot_pcl_",fiber_subset,".pdf",sep=""))
  names(for_plotting) <- rownames(random_forest_stats)[1:15]
  par(mar=c(25,4,1,1))
  b <- barplot(for_plotting,las=2)
  text(b,for_plotting,ifelse(is.na(overlap),"","*"),pos=3,cex=2,xpd=NA)
  dev.off()
}
