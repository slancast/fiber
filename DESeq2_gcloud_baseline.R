library(gplots)
library(RColorBrewer)
library(DESeq2)
library(readr)
library(tximport)
library(fdrtool)
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
##############################################################
##############################################################
##############################################################

library(stringi)

####### RNA-iPOP3 #######
####locate the directory containing the files####pri
dir <- "~/rnaseq-storage/rsem-results/iPOP-3"

###create vector of filenames from table
samples3 <- read.table(file.path(dir, "RNA_Fiber_ipop_3.csv"), sep=",", header = TRUE)

IDs <- substring(as.character(samples3$ID), 2, 9)
IDs_last_8 <- substring(as.character(samples3$ID), 10)
library(Biostrings)
dna = DNAStringSet(stri_reverse(IDs_last_8))
IDs_last_8_rc <- complement(dna)

##import transcript-level estimates from RSEM files
files3 <- file.path(dir, paste0(IDs, IDs_last_8_rc, "_R1_concatenated_", IDs, IDs_last_8_rc, "_R2_concatenated_rsem.genes.results"))

###
#tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
#head(tx2gene)


####### RNA-iPOP4 #######
####locate the directory containing the files####pri
dir <- "~/rnaseq-storage/rsem-results/iPOP-4"

###create vector of filenames from table
samples4 <- read.table(file.path(dir, "RNA_Fiber_ipop_4.csv"), sep=",", header = TRUE)

IDs <- substring(as.character(samples4$Barcode), 2, 9)
IDs_last_8 <- substring(as.character(samples4$Barcode), 10)
library(Biostrings)
dna = DNAStringSet(stri_reverse(IDs_last_8))
IDs_last_8_rc <- complement(dna)

##import transcript-level estimates from RSEM files
files4 <- file.path(dir, paste0(IDs, IDs_last_8_rc, "_R1_concatenated_", IDs, IDs_last_8_rc, "_R2_concatenated_rsem.genes.results"))

###
#tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
#head(tx2gene)

####### RNA-iPOP1_2 #######
####locate the directory containing the files####pri
dir <- "~/rnaseq-storage/rsem-results/iPOP-1-2"

###create vector of filenames from table
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

##import transcript-level estimates from RSEM files
files1_2 <- file.path(dir, paste0(samples$run, "_rsem.genes.results.txt"))


####### RNA-iPOP5-6-7 #######
####locate the directory containing the files####

dir <- "~/rnaseq-storage/rsem-results/iPOP-5-6-7"

samples5_7 <- read.table(file.path(dir, "ReplaceRNASeq_library_prep.csv"), sep=",", header = TRUE)
##import transcript-level estimates from RSEM files
files5_7 <- file.path(dir, paste0(samples5_7$adapters, ".genes.results"))



samplesa <- samples[,-1]
samples4a <- cbind(samples4[,7],samples4[,8:9],samples4[,6],samples4[,10])
colnames(samples4a) <- colnames(samplesa)
samples3a <- cbind(samples3[,7],samples3[,8:9],samples3[,6],samples3[,10])
colnames(samples3a) <- colnames(samplesa)
###
samples5_7a <- samples5_7[,6:10]
###
#
samplesa <- rbind(samplesa, samples3a, samples4a, samples5_7a)
rownames(samplesa) <- paste0("sample", 1:nrow(samplesa))

files <- c(files1_2,files3, files4, files5_7)
names(files) <- paste0("sample", 1:nrow(samplesa))
txi <- tximport(files, type = "rsem")
names(txi)

library(BiocParallel)
bpparam <- SnowParam(15, log = TRUE, stop.on.error = FALSE)
register(MulticoreParam(15))


###remove zero length transcripts
txi$length[txi$length == 0] <- 1

###create dataframe with proper rownames
sampleTable <- data.frame(condition = factor(samplesa$fiber))
sampleTable[,2] <- data.frame(week = factor(samplesa$week))
rownames(sampleTable) <- colnames(txi$counts)

########################################
#fiber_week
########################################

sampleTable[,3] <- data.frame(fiber_week = factor(interaction(sampleTable$condition,sampleTable$week)))
sampleTable$participant <- factor(samplesa$participant)
sampleTable$batch <- factor(samplesa$batch)
#sampleTable$fiber_week <- interaction(samples$fiber, samples$week)
sampleTable$fiber_week <- gsub(".*\\.Baseline","Baseline",sampleTable$fiber_week)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~participant + batch + fiber_week)
###prefilter low count genes

dds <- dds[ rowMeans(counts(dds)) > 0.15, ] #This cutoff really matters in determining the number of significant values

#Now I will write a function to pull out the raw counts with metadata from dds
metadata <- as.matrix(t(sampleTable))
rna_counts <- as.matrix(counts(dds))
raw_counts <- data.frame(rbind(metadata, rna_counts))
write.table(raw_counts, file="/home/slancast/rnaseq-storage/rsem-results/rna_raw_counts.csv",sep=",")


###set factor levels
#factor_levels <- unique(sampleTable$condition)
#dds$condition <- factor(dds$condition, levels = factor_levels)
#This loop is slightly clunky, but it works well enough for now. Next time perhaps I will
#create a loop for every factor. That should clean things up nicely.
#The first one will be for "condition" i.e. fiber type

i <- colnames(sampleTable)[3]
  print(i)
  factor_levels <- unique(sampleTable[[i]])
  ddsDESeq <- dds
  ddsDESeq[[i]] <- factor(ddsDESeq[[i]], levels = factor_levels)
j <- "Baseline"
    print(j)
    ddsDESeq[[i]] <- relevel(ddsDESeq[[i]], ref = j)  
    ddsDESeq <- DESeq(ddsDESeq, parallel = TRUE, BPPARAM=MulticoreParam(15))
    saveRDS(ddsDESeq,file="/home/slancast/rnaseq-storage/dds_baseline_participant.rds")
    ###Differential analysis
    print("releveled")
    print("names:")
    resultsNames(ddsDESeq)
    for (h in unique(sampleTable[[i]])){
      if (h == j) next
      print(h)
      res <- results(ddsDESeq, contrast = c("fiber_week",j, h), altHypothesis = "greaterAbs", parallel = TRUE, BPPARAM=MulticoreParam(15), cooksCutoff=FALSE, independentFiltering=FALSE)
      
      res <- na.omit(res)
      print("res created")
      res.fixed <- fdrtool(res$stat, statistic = "normal") #adjusting for batch effects (https://support.bioconductor.org/p/99685/; http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html)
      res.fixed <- data.frame(res.fixed)
      padjusted <- p.adjust(res.fixed$pval, method = "BH", n = nrow(res.fixed))
      res$pvalue <- res.fixed$pval
      res$padj <- padjusted
      res <- res[order(res$padj), ]
      head(res)
      
      print("LFCShrink")
      #resLFC <- lfcShrink(ddsDESeq, coef=2, res=res,parallel = TRUE, BPPARAM=MulticoreParam(15))
      print("LFCShrinkDone")
      #      names(resLFC)[1] <- "Gene"
      ## Examine independent filtering
      #attr(res, "filterThreshold")
      #plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
      #png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"resLFC-maplot.png",sep=""), 1500, 1000, pointsize=20)
      #plotMA(resLFC, ylim=c(-1,1), cex=1)
      #dev.off()
      #      res <- subset(res, baseMean > 0.2) #thowing out the extrememly low read count samples.
      
      
      resdata <- merge(as.data.frame(res), as.data.frame(counts(ddsDESeq, normalized=TRUE)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene"
      resdata <- entrezid(resdata)
      ## Write results
      write.csv(resdata, file=paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"-diffexpr-resultsb.csv",sep=""))
      
      ###Visualization and data collection
      
      # Plot dispersions
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"qc-dispersions.png",sep=""), 1000, 1000, pointsize=20)
      plotDispEsts(ddsDESeq, main="Dispersion plot")
      dev.off()
      
      # Regularized log transformation for clustering/heatmaps, etc
      print("rld")
      
      #rld <- bplapply(ddsDESeq, rlogTransformation(ddsDESeq), BPPARAM=MulticoreParam(15)) 
      #head(assay(rld))
      #hist(assay(rld))
      
      # Colors for plots below
      ## Ugly:
      ## (mycols <- 1:length(unique(condition)))
      ## Use RColorBrewer, better looking
      #(mycols <- brewer.pal(8, "Dark2")[1:length(unique(factor_levels))])
      
      # Sample distance heatmap
      #sampleDists <- as.matrix(dist(t(assay(rld))))
      #png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"qc-heatmap-samples.png",sep=""), w=1000, h=1000, pointsize=20)
      #heatmap.2(as.matrix(sampleDists), key=F, trace="none",
                #col=colorpanel(100, "green", "red"),
                #ColSideColors=mycols[ddsDESeq[[i]]], RowSideColors=mycols[ddsDESeq[[i]]],
                #margin=c(10, 10), main="Sample Distance Matrix")
      #dev.off()
      
      # Principal components analysis
      #png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"qc-pca.png",sep=""), 1000, 1000, pointsize=20)
      #rld_pca(rld, colors=mycols, intgroup="condition")
      #dev.off()
      
      
      # Get differential expression results
      ## Examine plot of p-values
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"res-pvalue.png",sep=""), 1000, 1000, pointsize=20)
      hist(res$pvalue, breaks=250, col="grey")
      dev.off()
      
      ## MA plot
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"diffexpr-maplot-noLFC.png",sep=""), 1500, 1000, pointsize=20)
      maplot(resdata, ylim = c(-1,1), main="MA Plot")
      dev.off()
      
      ## Volcano plot with "significant" genes labeled
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"diffexpr-volcanoplot.png",sep=""), 1200, 1000, pointsize=20)
      volcanoplot(resdata,  sigthresh=0.1, textcx=.8, xlim=c(-5, 5))
      dev.off()
    }
    
    system("sudo poweroff")

#####################################################
#Utils
#####################################################
    #The following will troubleshoot why some samples don't work well.
    #the cooks shows that the scinulin samples are outliers.
    filtering_table <- addmargins(table(filtering=(res$padj < .1), noFiltering=(resNoFilt$padj < .1)))
    if (FALSE) {
      write.csv(filtering_table, file=paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"filtering_table.csv",sep=""))
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"cooks_boxplot_1-50.png",sep=""), 1000, 1000, pointsize=20)
      boxplot(log10(assays(ddsDESeq)[["cooks"]][,1:50]), range=0, las=2)
      dev.off()
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"cooks_boxplot_51-100.png",sep=""), 1000, 1000, pointsize=20)
      boxplot(log10(assays(ddsDESeq)[["cooks"]][,51:100]), range=0, las=2)
      dev.off()
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"cooks_boxplot_101-150.png",sep=""), 1000, 1000, pointsize=20)
      boxplot(log10(assays(ddsDESeq)[["cooks"]][,101:150]), range=0, las=2)
      dev.off()
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"cooks_boxplot_151-200.png",sep=""), 1000, 1000, pointsize=20)
      boxplot(log10(assays(ddsDESeq)[["cooks"]][,151:200]), range=0, las=2)
      dev.off()
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"cooks_boxplot_201-250.png",sep=""), 1000, 1000, pointsize=20)
      boxplot(log10(assays(ddsDESeq)[["cooks"]][,201:250]), range=0, las=2)
      dev.off()
      png(paste("~/rnaseq-storage/rsem-results/deseq_results/",j,h,"cooks_boxplot_251-301.png",sep=""), 1000, 1000, pointsize=20)
      boxplot(log10(assays(ddsDESeq)[["cooks"]][,251:301]), range=0, las=2)
      dev.off()}


  