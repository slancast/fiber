
fiber_subset = "Arabinoxylan"

load(paste("~/NormAggreg_Log_",fiber_subset,"_cytokine_df.RData",sep=""))
load(paste("~/NormAggreg_Log_",fiber_subset,"_pcl_df.RData",sep=""))
load(paste("~/NormAggreg_Log_",fiber_subset,"_RNA_df.RData",sep=""))
load(paste("~/NormAggreg_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))

logaggregnorm_rna <- data.frame(t(logaggregnorm_rna))
logaggregnorm_rna <- cbind(logaggregnorm_rna, Gene = rownames(logaggregnorm_rna))


entrezid <- function( resdata ) {
  require(EnsDb.Hsapiens.v79)
  require(AnnotationDbi)  #the column to iterate over will be different if I'm using res vs resdata
  a = resdata$Gene #the column to iterate over will be different if I'm using res vs resdat
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
  
  resdata
}

logaggregnorm_rna <- entrezid(logaggregnorm_rna)
rownames(logaggregnorm_rna) <- make.names(logaggregnorm_rna$SYMBOL, unique=TRUE)
logaggregnorm_rna <- subset(logaggregnorm_rna, select=-c(SYMBOL,ENTREZID,Gene))
logaggregnorm_rna <- as.matrix(t(logaggregnorm_rna))


m <- rbind(t(logaggregnorm_pcl), t(logaggregnorm_rna), t(logaggregnorm_metabolomics), t(logaggregnorm_cytokine))
rownames(m) <- make.names(rownames(m), unique=TRUE)
lab <- c(rep("Microbiome",ncol(logaggregnorm_pcl)),rep("RNA",ncol(logaggregnorm_rna)),rep("Metabolites",ncol(logaggregnorm_metabolomics)),rep("Cytokine",ncol(logaggregnorm_cytokine)))
color <- c(rep("blue",ncol(logaggregnorm_pcl)),rep("red",ncol(logaggregnorm_rna)),rep("yellow",ncol(logaggregnorm_metabolomics)),rep("orange",ncol(logaggregnorm_cytokine)))
lab <- cbind(lab, color)
rownames(lab) <- rownames(m)

library(Mfuzz)
library(matrixStats)
n <- as.matrix(m)
class(n) <- "numeric"
#n <- n[rowVars(n)>0,] #the mfuzz packages do not like it when you have rows with variaiton of 0.
lab <- lab[row.names(lab) %in% rownames(n),]
rna_eset1 <- ExpressionSet(n)
rna_eset1 <- standardise(rna_eset1) #Running standarise 
o <- exprs(rna_eset1)
o <- na.omit(o)
lab <- lab[row.names(lab) %in% rownames(o),]


library("Hmisc")
# Plot correlation graph
cor <- rcorr(format(t(o),digits=20), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
FDR.data <- apply(pval.data,2,p.adjust,method="bonferroni", n = length(pval.data))
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
save(cor.data, file="bonderroni_corrected_cor.data.RData")
#load("spear_bonferroni__corrected_cor.data.RData")
#cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
library("igraph")
network=graph.adjacency(cor.data, weighted=T, mode="undirected", diag=F)
V(network)$color <- lab[,2]
print(network)

library(icesTAF)
mkdir(paste("./lognormaggreg",fiber_subset,sep=""))
##########################################################################
#Printing subgraphs by name###############################################
##########################################################################
#finding subcluster by vertex name

for (i in colnames(logaggregnorm_cytokine)) {
  print(i)
  vertex <- match(i, V(network)$name)
  cluster_number <- clusters(network)$membership[vertex]#}
  dg <- decompose.graph(network)
  subgraph <- dg[[cluster_number]]
  ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
  
  pdf(paste("./lognormaggreg",fiber_subset,"/",i,"_subgraph.pdf",sep=""))
  par(bg="white", mar=c(0,0,0,0))
  set.seed(4)
  plot(subgraph,
       vertex.size=10,
       vertex.color=V(subgraph)$color,
       vertex.label.cex=0.5,
       vertex.label.color="black",
       vertex.frame.color="black",
       #vertex.label = NA,
       layout = ly
  )
  dev.off()
}

##########################################################################
#Printing All subgraphs###################################################
##########################################################################

if (FALSE) {
  dg <- decompose.graph(network)
  
  for (j in 1:length(dg)) {
    subgraph <- dg[[j]]
    counter <- 0
    ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
    
    pdf(paste("./lognormaggreg",fiber_subset,"/",j,"_subgraph.pdf",sep=""))
    par(bg="white", mar=c(0,0,0,0))
    set.seed(4)
    plot(subgraph,
         vertex.size=10,
         vertex.color=V(subgraph)$color,
         vertex.label.cex=0.5,
         vertex.label.color="black",
         vertex.frame.color="black",
         #vertex.label = NA,
         layout = ly
    )
    dev.off()
    
