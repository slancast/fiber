##Load libraries
library("tximport")
library("readr")

####locate the directory containing the files####
dir <- "/Users/SLancaster/Desktop/rsem_genes_results_all_trimmed"

###create vector of filenames from table
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

##import transcript-level estimates from RSEM files
files <- file.path(dir, "rsem", samples$run, paste0(samples$run, ".genes.results"))
names(files) <- paste0("sample", 1:5)
txi.rsem <- tximport(files, type = "rsem")
head(txi.rsem$counts)

###
#tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
#head(tx2gene)

###
library(tximport)
library(readr)
txi <- tximport(files, type = "rsem")
names(txi)

###remove zero length transcripts
txi$length[txi$length == 0] <- 1

###DEseq2
library(DESeq2)

###create dataframe with proper rownames
sampleTable <- data.frame(condition = factor(samples$condition))
sampleTable[,2] <- data.frame(sensitivity = factor(samples$sensitivity))
rownames(sampleTable) <- colnames(txi$counts)

###construct a DESeqDataSet from the txi object and sample info in samples.txt file
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~sensitivity)

###prefilter low count genes
dds <- dds[ rowSums(counts(dds)) > 1, ]

###set factor levels
dds$sensitivity <- factor(dds$sensitivity, levels = c("IR","IS"))

###Differential analysis
dds <- DESeq(dds)
res <- results(dds)

###Create a copy of dds to rerun the analysis using a multi-factor design
#ddsMF <- dds
#levels(ddsMF$condition)
#levels(ddsMF$sensitivity)

###Rerun deseq2. Since sensitivity is last, it will be contrasted
#design(ddsMF) <- formula(~ condition + sensitivity)
#ddsMF <- DESeq(ddsMF)

#resMF <- results(ddsMF)
#head(resMF)

#remove decimal and numbers after from ensg id's in res
tmp=gsub("\\..*","",row.names(res))

#Add gene ids from ensembl to dds res file
library(org.Hs.eg.db)

convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

#Convert Ensembl IDs in the rownames of res to gene symbols 
res$hgnc_symbol <- convertIDs( tmp, "ENSEMBL", "SYMBOL", org.Hs.eg.db )
res$entrezid <- convertIDs( row.names(res), "ENSEMBL", "ENTREZID", org.Hs.eg.db )
head(res, 4)


# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

# Get differential expression results for multifactor experiment
resMF <- results(ddsMF)
table(resMF$padj<0.05)
## Order by adjusted p-value
resMF <- resMF[order(resMF$padj), ]
## Merge with normalized count data
resMFdata <- merge(as.data.frame(resMF), as.data.frame(counts(ddsMF, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMFdata)[1] <- "Gene"
head(resMFdata)
## Write results
write.csv(resMFdata, file="diffexpr-results.csv")

##Gene-set enrichment analysis
#source("http://bioconductor.org/biocLite.R")
#biocLite("reactome.db")
library( "reactome.db" )

  #subset the results table to only those genes for which the Reactome database has data
  res2 <- res[ res$entrezid %in% keys( reactome.db, "ENTREZID" ) & !is.na( res$padj ) , ]
  head(res2)
  
  #get a table with the mapping from Entrez IDs to Reactome Path IDs
  reactomeTable <- AnnotationDbi::select( reactome.db, 
                                          keys=as.character(res2$entrezid), keytype="ENTREZID", 
                                          columns=c("ENTREZID","REACTOMEID") )
  head(reactomeTable)
  
  #Transform above table to an incident matrix
  incm <- do.call( rbind, with(reactomeTable, tapply( 
    ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x ) ))
  colnames(incm) <- res2$entrez
  str(incm)
  
  #remove all rows corresponding to Reactome Paths with less than 20 or more than 80 assigned genes
  within <- function(x, lower, upper) (x>=lower & x<=upper)
  incm <- incm[ within(rowSums(incm), lower=20, upper=80), ]
  
  #Statistics for GSEA
  testCategory <- function( reactomeID ) {
    isMember <- incm[ reactomeID, ]
    data.frame( 
      reactomeID  = reactomeID,
      numGenes    = sum( isMember ),
      avgLFC      = mean( res2$log2FoldChange[isMember] ),
      sdLFC       = sd( res2$log2FoldChange[isMember] ),
      zValue      = mean( res2$log2FoldChange[isMember] ) /sd( res2$log2FoldChange[isMember] ),
      strength    = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
      pvalue      = t.test( res2$log2FoldChange[ isMember ] )$p.value,
      reactomeName = reactomePATHID2NAME[[reactomeID]],
      stringsAsFactors = FALSE ) }
  
  
  #Call the function for all Paths in our incidence matrix and collect the results
  reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )
  reactomeResult$padjust <- p.adjust( reactomeResult$pvalue, "BH" )
  
  reactomeResultSignif <- reactomeResult[ reactomeResult$padjust < 0.05, ]
  head( reactomeResultSignif[ order(-reactomeResultSignif$strength), ] )
  
  
  
  
  
  #Tximport level gene id insertion

  k <- keys(txdb, keytype = "GENEID")
  df <- select(txdb, keys = k, keytype = "GENEID", columns = "ENTREZID")
  tx2gene <- df[, 1:2]  # tx ID, then gene ID
  colnames(tx2gene) <- c("TXNAME", "GENEID")
  tx2gene$ENTREZID <- convertIDs( tx2gene[,2], "SYMBOL", "ENTREZID", EnsDb.Hsapiens.v79 )

  
############################################################
#This is to append the entrez IDs into the res files like Charles asked for
  library(EnsDb.Hsapiens.v79)
entrezid <- function( res ) {
  txdb <- EnsDb.Hsapiens.v79 
  k <- row.names(res)
  df <- select(txdb, keys = tmp, keytype = "GENEID", columns = "ENTREZID")
  row.names(res) <- tmp
 
 ENTREZID <- c() 
 counter1 <- 0
 for (i in row.names(res)) {
   counter <- 0
   counter1 <- counter1 + 1
   j <- match(i,df$GENEID)
   ENTREZID <- c(ENTREZID, toString(df[j,][2]))}
  res$ENTREZ_ID <- ENTREZID 
  res
}
############################################################