#functions here

#For some of the measurements we have repeats. I think mainly Charles did this in order to see the batch effects.
#However this is probably unnecessary if the data are truly randomized, and then we can perform the best batch
#correction using Combat or something like that
#The input matrix should have the first column as the participant and the second column as the timepoint
averaging_replicates <- function(fiber_subsetted_dataframe) {
df <- data.frame(fiber_subsetted_dataframe)
colnames(df)[1] <- "Participant"
colnames(df)[2] <- "Time_point"
df[,3:ncol(df)] <- data.frame(lapply(df[,3:ncol(df)], as.character), stringsAsFactors=FALSE, row.names = rownames(df)) #To change to a double, this needs to go through character first
df[,3:ncol(df)] <- data.frame(lapply(df[,3:ncol(df)], as.numeric), stringsAsFactors=FALSE, row.names = rownames(df)) #This needs to be done for the transcript data, but causes problems with the pcl data.

df2 <- aggregate(df,list(Participant=df[,1], Visit=df[,2]), mean)
df2 <- data.frame(t(df2)) #popping out the analytes that have NAs
df2 <- na.omit(df2)
df2 <- data.frame(t(df2))
rownames(df2) <- paste(as.character(df2$Participant), as.character(df2$Visit),sep="_")
return(df2)
}

# A standatdized way of determining the correlations
correlations <- function(to_correlate_df) {
  require(Hmisc)
  cor <- rcorr(format(t(to_correlate_df),digits=20), type="spearman")
  cor.data <- cor$r
  cor.data[upper.tri(cor.data, diag = T)]<- 0
  pval.data <- cor$P
  pval.data[upper.tri(pval.data, diag = T)]<- NA
  FDR.data <- apply(pval.data,2,p.adjust,method="BH", n = length(pval.data))
  pdf("./pval_bonferonni_hist.pdf")
  hist(pval.data, breaks = 100, col="darkblue")
  dev.off()
  pdf("./FDR_bonferonni_hist.pdf")
  hist(FDR.data, breaks = 100, col="darkblue")
  dev.off()
  pdf("./cor_bonferonni_hist.pdf")
  hist(cor.data, breaks = 10, col="red")
  dev.off()
  cor.data[FDR.data > 0.05]=0
  return(cor.data)
}

# A standard way of running the standardise funciton
standardise_data <- function(data_frame) {
  require(Mfuzz)
  n <- as.matrix(data_frame)
  class(n) <- "numeric"
  eset1 <- ExpressionSet(n)
  eset1 <- standardise(eset1) #Running standarise 
  o <- exprs(eset1)
  o <- na.omit(o)
  return(o)
}

# This function will perform mean imputation
# For the knn imputation, I can just use the funciton impute.knn
mean_imputation <- function(combined_df){
  imputed_df <- data.frame(matrix(nrow = 0, ncol = 6)) #creating a matrix for the imputed values. I originally ran it on the non-imputed values, but I think I might as well add them
  require(Xmisc)
  combined_df <- as.matrix(combined_df)
  class(combined_df) <- "numeric"
  for (i in 1:nrow(combined_df)) {
    j <- combined_df[i,]
    j[is.na(j)] <- mean(j, na.rm = TRUE)
    imputed_df[i,] <- j
  }
  rownames(imputed_df) <- rownames(combined_df)
}

plotting_fuzzy_clusters <- function(eset, mfuzzcl){
require(Mfuzz)
xaxis_ticks = c("B","10","20","30","D3","D10")
par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), ylim=c(-3,3),time.labels = xaxis_ticks, bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
}

rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.1, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#Must be run with resdata, not res! This is because of the lab=Gene
maplot <- function (res, thresh=0.1, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}


entrezid <- function( resdata ) {
  require(EnsDb.Hsapiens.v79)
  a = resdata$Gene #the column to iterate over will be different if I'm using res vs resdata
  tmp=gsub("\\..*","",a)
  tmp <- as.character(tmp)
  txdb <- EnsDb.Hsapiens.v79 
  df <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "ENTREZID")
  df2 <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "SYMBOL")
  
  ENTREZID <- c() 
  SYMBOL <- c() 
  counter1 <- 0
  for (i in tmp) { 
    counter1 <- counter1 + 1
    j <- match(i,df$GENEID)
    ENTREZID <- c(ENTREZID, toString(df[j,][2]))
    SYMBOL <- c(SYMBOL, toString(df2[j,][2]))}
  resdata$ENTREZID <- ENTREZID
  resdata$SYMBOL <- SYMBOL
  resdata$EnsemblGene <- tmp
  
  resdata
}
