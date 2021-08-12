library(tximport)
require(genefilter)
require(calibrate)
library(readr)
library(EnsDb.Hsapiens.v79)
library(DESeq2)
library(RColorBrewer)
library(gplots)
options(bitmapType='cairo') #work around for X11 missing on cluster

##############################################################
##############################################################
##############################################################
#functions here
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

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
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


maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}


entrezid <- function( resdata ) {
  a = resdata$Gene #the column to iterate over will be different if I'm using res vs resdata
  tmp=gsub("\\..*","",a)
  tmp <- as.character(tmp)
  txdb <- EnsDb.Hsapiens.v79 
  df <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "ENTREZID")
  
  ENTREZID <- c() 
  counter1 <- 0
  for (i in tmp) { 
    counter1 <- counter1 + 1
    j <- match(i,df$GENEID)
    ENTREZID <- c(ENTREZID, toString(df[j,][2]))}
  resdata$ENTREZID <- ENTREZID 
  resdata
}
##############################################################
##############################################################
##############################################################

####locate the directory containing the files####
dir <- "/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/rsem_genes_results_all_trimmed"
#dir <- "/Users/SLancaster/Desktop/rsem_genes_results_all_trimmed"

###create vector of filenames from table
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

##import transcript-level estimates from RSEM files
files <- file.path(dir, samples$run, paste0(samples$run, ".genes.results"))
names(files) <- paste0("sample", 1:nrow(samples))

###
#tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
#head(tx2gene)

###
txi <- tximport(files, type = "rsem")
names(txi)

###remove zero length transcripts
txi$length[txi$length == 0] <- 1

###create dataframe with proper rownames
sampleTable <- data.frame(condition = factor(samples$fiber))
sampleTable[,2] <- data.frame(week = factor(samples$week))
rownames(sampleTable) <- colnames(txi$counts)

###construct a DESeqDataSet from the txi object and sample info in samples.txt file
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition + week)
###prefilter low count genes
dds <- dds[ rowSums(counts(dds)) > 1, ]


###set factor levels
#factor_levels <- unique(sampleTable$condition)
#dds$condition <- factor(dds$condition, levels = factor_levels)
#This loop is slightly clunky, but it works well enough for now. Next time perhaps I will
#create a loop for every factor. That should clean things up nicely.
#The first one will be for "condition" i.e. fiber type
for (i in colnames(sampleTable)[1]) {
  print(i)
  factor_levels <- unique(sampleTable[[i]])
  print("factor levels:")
  print(factor_levels)
  ddsDESeq <- dds
  ddsDESeq[[i]] <- factor(ddsDESeq[[i]], levels = factor_levels)
  for (j in unique(sampleTable[[i]])) { 
    print(j)
    ddsDESeq[[i]] <- relevel(ddsDESeq[[i]], ref = j)   
    ###Differential analysis
    
    ddsDESeq <- DESeq(ddsDESeq)
    
    print("Results Names:")
    print(resultsNames(ddsDESeq))
    for (h in resultsNames(ddsDESeq)[2:4]){
      print(h)
      res <- results(ddsDESeq, contrast = list(h), altHypothesis = "greaterAbs")
      table(res$padj<0.05)
      res <- res[order(res$padj), ]
      head(res)
      res <- subset(res, baseMean > 0.1) #thowing out the extrememly low read count samples.
      resdata <- merge(as.data.frame(res), as.data.frame(counts(ddsDESeq, normalized=TRUE)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene"
      ## Write results
      resdata <- entrezid(resdata)
      write.csv(resdata, file=paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"-diffexpr-results.csv",sep=""))
      
      ###Visualization and data collection
      
      # Plot dispersions
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-dispersions.png",sep=""), 1000, 1000, pointsize=20)
      plotDispEsts(ddsDESeq, main="Dispersion plot")
      dev.off()
      
      # Regularized log transformation for clustering/heatmaps, etc
      print("rld")
      rld <- rlogTransformation(ddsDESeq, blind=FALSE)
      vsd <- varianceStabilizingTransformation(ddsDESeq, blind=FALSE) #blind=false because dds has already been run estimating dispersion values
      head(assay(rld))
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"rld_hist.png",sep=""))
      hist(assay(rld))
      dev.off()
      print("rld finsihed")
      
      # Colors for plots below
      ## Ugly:
      ## (mycols <- 1:length(unique(condition)))
      ## Use RColorBrewer, better looking
      (mycols <- brewer.pal(8, "Dark2")[1:length(unique(factor_levels))])
      
      # Sample distance heatmap
      sampleDists <- as.matrix(dist(t(assay(rld))))
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-heatmap-samples.png",sep=""), w=1000, h=1000, pointsize=20)
      heatmap.2(as.matrix(sampleDists), key=F, trace="none",
                col=colorpanel(100, "green", "red"),
                ColSideColors=mycols[ddsDESeq[[i]]], RowSideColors=mycols[ddsDESeq[[i]]],
                margin=c(10, 10), main="Sample Distance Matrix")
      dev.off()
      
      png(paste("/Users/SLancaster/Desktop/RNAseq_results/",h,"qc-vsd-pca.png",sep=""), 1000, 1000, pointsize=20)
      rld_pca(vsd, colors=mycols, intgroup="condition", xlim=c(-75, 35))
      dev.off()
      
      # Principal components analysis
      ## Could do with built-in DESeq2 function:
      ## DESeq2::plotPCA(rld, intgroup="condition")
      ## I like this better:

      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-pca.png",sep=""), 1000, 1000, pointsize=20)
      rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
      dev.off()
      
      ## Examine plot of p-values
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"res-pvalue.png",sep=""), 1000, 1000, pointsize=20)
      hist(res$pvalue, breaks=50, col="grey")
      dev.off()
      
      ## Examine independent filtering
      #attr(res, "filterThreshold")
      #plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
      
      ## MA plot
      ## Could do with built-in DESeq2 function:
      ## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
      ## Use this instead though:
      print("resLFC")
#      resLFC <- lfcShrink(ddsDESeq, coef=2, res=res)
#      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"built-in-maplot.png",sep=""), 1500, 1000, pointsize=20)
#      plotMA(resLFC, ylim=c(-1,1), cex=1)
#     dev.off()

#      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-maplot.png",sep=""), 1500, 1000, pointsize=20)
#      maplot(resLFC, ylim = c(-.2,.2), main="MA Plot")
#      dev.off()
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-maplot-noLFC.png",sep=""), 1500, 1000, pointsize=20)
      maplot(resdata, ylim = c(-1,1), main="MA Plot")
      dev.off()
      
      ## Volcano plot with "significant" genes labeled
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-volcanoplot.png",sep=""), 1200, 1000, pointsize=20)
      volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-1, 1))
      dev.off()
    }}}

########################################
#The second one is now for column "week"
########################################
  
  for (i in colnames(sampleTable)[2]) {
  print(i)
  factor_levels <- unique(sampleTable[[i]])
  ddsDESeq <- dds
  ddsDESeq[[i]] <- factor(ddsDESeq[[i]], levels = factor_levels)
  for (j in unique(sampleTable[[i]])) { 
  print(j)
  ddsDESeq[[i]] <- relevel(ddsDESeq[[i]], ref = j)   
  ###Differential analysis

  ddsDESeq <- DESeq(ddsDESeq)
  
  resultsNames(ddsDESeq)
  for (h in resultsNames(ddsDESeq)[5:length(resultsNames(ddsDESeq))]){
  print(h)
  res <- results(ddsDESeq, lfcThreshold = 1, contrast = list(h), altHypothesis = "greaterAbs")
  table(res$padj<0.05)
  res <- res[order(res$padj), ]
  head(res)
  res <- subset(res, baseMean > 0.1) #thowing out the extrememly low read count samples
  resdata <- merge(as.data.frame(res), as.data.frame(counts(ddsDESeq, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  head(resdata)
  ## Write results
  write.csv(resdata, file=paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"-diffexpr-results.csv",sep=""))
                 
  ###Visualization and data collection
  
  # Plot dispersions
  png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-dispersions.png",sep=""), 1000, 1000, pointsize=20)
  plotDispEsts(ddsDESeq, main="Dispersion plot")
  dev.off()
  
  # Regularized log transformation for clustering/heatmaps, etc
  print("rld")
  rld <- rlogTransformation(ddsDESeq)
  head(assay(rld))
  hist(assay(rld))
  
  # Colors for plots below
  ## Ugly:
  ## (mycols <- 1:length(unique(condition)))
  ## Use RColorBrewer, better looking
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(factor_levels))])
  
  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(assay(rld))))
  png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-heatmap-samples.png",sep=""), w=1000, h=1000, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[ddsDESeq[[i]]], RowSideColors=mycols[ddsDESeq[[i]]],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  
  # Principal components analysis

  png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-pca.png",sep=""), 1000, 1000, pointsize=20)
  rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
  dev.off()
  
  
  # Get differential expression results
  
  ## Examine plot of p-values
  png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"res-pvalue.png",sep=""), 1000, 1000, pointsize=20)
  hist(res$pvalue, breaks=50, col="grey")
  dev.off()
  
  ## Examine independent filtering
  #attr(res, "filterThreshold")
  #plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
  
  ## MA plot
  png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-maplot-noLFC.png",sep=""), 1500, 1000, pointsize=20)
  maplot(resdata, ylim = c(-1,1), main="MA Plot")
  dev.off()
  
  ## Volcano plot with "significant" genes labeled
  png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-volcanoplot.png",sep=""), 1200, 1000, pointsize=20)
  volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-1, 1))
  dev.off()
  }}}
  
  

  
#dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ week:condition)
sampleTable[,3] <- data.frame(fiber_week = factor(interaction(sampleTable$condition,sampleTable$week)))
#sampleTable$fiber_week <- interaction(samples$fiber, samples$week)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~fiber_week)
###prefilter low count genes
dds <- dds[ rowSums(counts(dds)) > 1, ]

###set factor levels
#factor_levels <- unique(sampleTable$condition)
#dds$condition <- factor(dds$condition, levels = factor_levels)
#This loop is slightly clunky, but it works well enough for now. Next time perhaps I will
#create a loop for every factor. That should clean things up nicely.
#The first one will be for "condition" i.e. fiber type
for (i in colnames(sampleTable)[3]) {
  print(i)
  factor_levels <- unique(sampleTable[[i]])
  ddsDESeq <- dds
  print(factor_levels)
  ddsDESeq[[i]] <- factor(ddsDESeq[[i]], levels = factor_levels)
  for (j in unique(sampleTable[[i]])) { 
    print(j)
    ddsDESeq[[i]] <- relevel(ddsDESeq[[i]], ref = j)   
    ###Differential analysis
    
    ddsDESeq <- DESeq(ddsDESeq)
    
    resultsNames(ddsDESeq)
    for (h in resultsNames(ddsDESeq)[11:length(resultsNames(ddsDESeq))]){
      print(h)
      res <- results(ddsDESeq, contrast = list(h), lfcThreshold = 1, altHypothesis = "greaterAbs")
      table(res$padj<0.05)
      res <- res[order(res$padj), ]
      head(res)
      res <- subset(res, baseMean > 0.1) #thowing out the extrememly low read count samples.
      resdata <- merge(as.data.frame(res), as.data.frame(counts(ddsDESeq, normalized=TRUE)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene"
      head(resdata)
      ## Write results
      write.csv(resdata, file=paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"-diffexpr-results.csv",sep=""))
      
      ###Visualization and data collection
      
      # Plot dispersions
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-dispersions.png",sep=""), 1000, 1000, pointsize=20)
      plotDispEsts(ddsDESeq, main="Dispersion plot")
      dev.off()
      
      # Regularized log transformation for clustering/heatmaps, etc
      print("rld")
      rld <- rlogTransformation(ddsDESeq)
      head(assay(rld))
      hist(assay(rld))
      
      # Colors for plots below
      ## Ugly:
      ## (mycols <- 1:length(unique(condition)))
      ## Use RColorBrewer, better looking
      (mycols <- brewer.pal(8, "Dark2")[1:length(unique(factor_levels))])
      
      # Sample distance heatmap
      sampleDists <- as.matrix(dist(t(assay(rld))))
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-heatmap-samples.png",sep=""), w=1000, h=1000, pointsize=20)
      heatmap.2(as.matrix(sampleDists), key=F, trace="none",
                col=colorpanel(100, "black", "white"),
                ColSideColors=mycols[ddsDESeq[[i]]], RowSideColors=mycols[ddsDESeq[[i]]],
                margin=c(10, 10), main="Sample Distance Matrix")
      dev.off()
      
      # Principal components analysis
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"qc-pca.png",sep=""), 1000, 1000, pointsize=20)
      rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
      dev.off()
      
      
      # Get differential expression results
      ## Examine plot of p-values
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"res-pvalue.png",sep=""), 1000, 1000, pointsize=20)
      hist(res$pvalue, breaks=50, col="grey")
      dev.off()
      
      ## Examine independent filtering
      #attr(res, "filterThreshold")
      #plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
      
      ## MA plot
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-maplot-noLFC.png",sep=""), 1500, 1000, pointsize=20)
      maplot(resdata, ylim = c(-1,1), main="MA Plot")
      dev.off()
      
      ## Volcano plot with "significant" genes labeled
      png(paste("/srv/gsfs0/projects/snyder/slancast/fiber/rnaseq/RNAseq_results/",h,"diffexpr-volcanoplot.png",sep=""), 1200, 1000, pointsize=20)
      volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-1, 1))
      dev.off()
    }}}

