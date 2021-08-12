load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_genef_df.RData",sep=""), envir = parent.frame())
if (FALSE) {
  
print("combined")
class(combined_df) <- "numeric"
library(matrixStats)
combined_df <- combined_df[rowVars(combined_df)>0,]
set.seed(1)
library(Mfuzz)
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

write.table(cluster_membership, file = paste("/home/slancast/RNA/clusters",z,"multiomics_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")

library(RDAVIDWebService)
david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setTimeOut(david, 100000)
background <- addList(david, cluster_membership$EnsemblGene, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")

for (i in unique(cluster_membership$mfuzzcl.cluster)) {
  print(fiber_subset)
  print(z)
  print(i)
  assign(paste("cluster",as.character(i),sep=""), cluster_membership[ which(cluster_membership$mfuzzcl.cluster==i),])
  write.table(get(paste("cluster",as.character(i),sep="")), file = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,".txt",sep=""),sep="\t")
  matrix <- get(paste("cluster",as.character(i),sep=""))
  result<-addList(david, matrix$EnsemblGene, idType="ENSEMBL_GENE_ID",listName=paste("cluster",as.character(i),sep=""), listType="Gene")
  getFunctionalAnnotationChartFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationChart.txt",sep=""))
  getFunctionalAnnotationTableFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationTable.txt",sep=""))
  getClusterReportFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"ClusterReport.txt",sep=""))
} }

}

# cluster_number <- Dmin(eset,m1,crange=seq(2,15,1),repeats=5,visu=TRUE)
# cluster_selection <- cselection(eset,m1,crange=seq(4,20,2),repeats=5,visu=TRUE) }
################################################################################################################
################################################End combined Start RNA##########################################
################################################################################################################

if (FALSE) {
print("rna")
combined_df <- t(logaggregnorm_rna)
class(combined_df) <- "numeric"
library(matrixStats)
combined_df <- combined_df[rowVars(combined_df)>0,]
set.seed(1)
library(Mfuzz)
eset <- ExpressionSet(combined_df) #Creating the type expression set with the metadata rows as a different argument
m <- exprs(eset)
eset <- standardise(eset) #Running standarise 
m <- exprs(eset)
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
  
  write.table(cluster_membership, file = paste("/home/slancast/RNA/clusters",z,"rna_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")
  
  
  library(RDAVIDWebService)
  david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setTimeOut(david, 100000)
  background <- addList(david, cluster_membership$EnsemblGene, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")
  
  for (i in unique(cluster_membership$mfuzzcl.cluster)) {
    print(fiber_subset)
    print(z)
    print(i)
    assign(paste("cluster",as.character(i),sep=""), cluster_membership[ which(cluster_membership$mfuzzcl.cluster==i),])
    write.table(get(paste("cluster",as.character(i),sep="")), file = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,".txt",sep=""),sep="\t")
    matrix <- get(paste("cluster",as.character(i),sep=""))
    result<-addList(david, matrix$EnsemblGene, idType="ENSEMBL_GENE_ID",listName=paste("cluster",as.character(i),sep=""), listType="Gene")
    getFunctionalAnnotationChartFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationChart.txt",sep=""))
    getFunctionalAnnotationTableFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationTable.txt",sep=""))
    getClusterReportFile(david, fileName = paste("/home/slancast/RNA/clusters",z,"cluster",as.character(i),"_",fiber_subset,"ClusterReport.txt",sep=""))
  } }
  
}
################################################################################################################
################################################End RNA Start PCL###############################################
################################################################################################################
  
if (FALSE) {
  print("pcl")
  combined_df <- t(logaggregnorm_pcl)
  class(combined_df) <- "numeric"
  library(matrixStats)
  combined_df <- combined_df[rowVars(combined_df)>0,]
  set.seed(1)
  library(Mfuzz)
  eset <- ExpressionSet(combined_df) #Creating the type expression set with the metadata rows as a different argument
  m <- exprs(eset)
  eset <- standardise(eset) #Running standarise 
  m <- exprs(eset)
  m1 = mestimate(eset)
  
  for (z in 2:16) {
    mfuzzcl <- mfuzz(eset, c=z, m=m1)
    
    pdf(paste("/home/slancast/PCL/clusters",z,"multiomics_",fiber_subset,"_parcoord.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
    dev.off()
    
    pdf(paste("/home/slancast/PCL/clusters",z,"multiomics_",fiber_subset,"_parcoord2.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
    dev.off()
    
    cluster_membership <- data.frame(mfuzzcl$cluster)
    
    write.table(cluster_membership, file = paste("/home/slancast/PCL/clusters",z,"pcl_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")
    }
 
}

################################################################################################################
################################################End PCL Start Genef#############################################
################################################################################################################

print("genef")
logaggregnorm_genef <- t(logaggregnorm_genef)
class(logaggregnorm_genef) <- "numeric"
library(matrixStats)
logaggregnorm_genef <- logaggregnorm_genef[rowVars(logaggregnorm_genef)>0,]
set.seed(1)
library(Mfuzz)
eset <- ExpressionSet(logaggregnorm_genef) #Creating the type expression set with the metadata rows as a different argument
m <- exprs(eset)
eset <- standardise(eset) #Running standarise 
m <- exprs(eset)
m1 = mestimate(eset)

for (z in 2:16) {
  mfuzzcl <- mfuzz(eset, c=z, m=m1)
  
  pdf(paste("/home/slancast/genef/clusters",z,"multiomics_",fiber_subset,"_parcoord.pdf",sep=""))
  par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
      mar = c(2,1,1,1) + 0.1)
  mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
  dev.off()
  
  pdf(paste("/home/slancast/genef/clusters",z,"multiomics_",fiber_subset,"_parcoord2.pdf",sep=""))
  par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
      mar = c(2,1,1,1) + 0.1)
  mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
  dev.off()
  
  cluster_membership <- data.frame(mfuzzcl$cluster)
  
  write.table(cluster_membership, file = paste("/home/slancast/genef/clusters",z,"genef_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")
  }

  ################################################################################################################
  ################################################End PCL Start Metabolome #######################################
  ################################################################################################################
  
  if (FALSE) {
  print("metabolome")
  combined_df <- t(logaggregnorm_metabolomics)
  class(combined_df) <- "numeric"
  library(matrixStats)
  combined_df <- combined_df[rowVars(combined_df)>0,]
  set.seed(1)
  library(Mfuzz)
  eset <- ExpressionSet(combined_df) #Creating the type expression set with the metadata rows as a different argument
  m <- exprs(eset)
  eset <- standardise(eset) #Running standarise 
  m <- exprs(eset)
  m1 = mestimate(eset)
  
  for (z in 2:16) {
    mfuzzcl <- mfuzz(eset, c=z, m=m1)
    
    pdf(paste("/home/slancast/Metabolomics/clusters",z,"multiomics_",fiber_subset,"_parcoord.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
    dev.off()
    
    pdf(paste("/home/slancast/Metabolomics/clusters",z,"multiomics_",fiber_subset,"_parcoord2.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(eset,cl=mfuzzcl,mfrow=c(4,4), bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
    dev.off()
    
    cluster_membership <- data.frame(mfuzzcl$cluster)
    
    write.table(cluster_membership, file = paste("/home/slancast/Metabolomics/clusters",z,"metabolomics_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")
  }
  
  }
  
} 

