
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
  resdata$Gene_woversion <- tmp
  
  resdata
}
options(java.parameters = "-Xmx100g")#This is required otherwise java will run out of memory. The 100g stands for 100 gigs.
library(Mfuzz)
library(matrixStats)
library(RDAVIDWebService)

for (fiber_subset in c("Arabinoxylan", "LCInulin", "Mix")) {
print(paste("fiber_subset: ",fiber_subset,sep="") ) 

load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_RNA_df.RData",sep=""), envir = parent.frame())

combined_df <- t(logaggregnorm_rna)
class(combined_df) <- "numeric"
combined_df <- combined_df[rowVars(combined_df)>0,]
set.seed(1)
eset <- ExpressionSet(combined_df) #Creating the type expression set with the metadata rows as a different argument
m <- exprs(eset)
eset <- standardise(eset) #Running standarise 
m1 = mestimate(eset)

for (z in 2:16) {
  mfuzzcl <- mfuzz(eset, c=z, m=m1)
  
  pdf(paste("/home/slancast/RNA/clusters",z,"multiomics_",fiber_subset,"_parcoord.pdf",sep=""))
  par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
      mar = c(2,1,1,1) + 0.1)
  mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
  dev.off()
  
  pdf(paste("/home/slancast/RNA/clusters",z,"multiomics_",fiber_subset,"_parcoord2.pdf",sep=""))
  par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
      mar = c(2,1,1,1) + 0.1)
  mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
  dev.off()
  
  cluster_membership <- data.frame(mfuzzcl$cluster)
  cluster_membership$Gene <- rownames(cluster_membership)
  cluster_membership <- entrezid(cluster_membership)
  cluster_membership
  
  write.table(cluster_membership, file = paste("/home/slancast/RNA/clusters",z,"multiomics_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")
  
  david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setTimeOut(david, 1000000)
  setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
  background <- addList(david, cluster_membership$Gene_woversion, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")
  
  for (i in unique(cluster_membership$mfuzzcl.cluster)) {
    print(fiber_subset)
    print(z)
    print(i)
    assign(paste("cluster",as.character(i),sep=""), cluster_membership[ which(cluster_membership$mfuzzcl.cluster==i),])
    write.table(get(paste("cluster",as.character(i),sep="")), file = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,".txt",sep=""),sep="\t")
    matrix <- get(paste("cluster",as.character(i),sep=""))
    result<-addList(david, matrix$Gene_woversion, idType="ENSEMBL_GENE_ID",listName=paste("cluster",as.character(i),sep=""), listType="Gene")
    getFunctionalAnnotationChartFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationChart.txt",sep=""))
    getFunctionalAnnotationTableFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationTable.txt",sep=""))
    getClusterReportFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"ClusterReport.txt",sep=""))
  } }

}#Ending fiber subset
