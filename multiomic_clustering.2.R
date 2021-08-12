#
#

#
for (fiber_subset in c( "Arabinoxylan", "LCInulin","Mix")) {

print(paste("fiber_subset: ",fiber_subset,sep="") ) 
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_cytokine_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_pathabundance_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_RNA_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_metabolomics_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_lipids_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_proteomics_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_metaphlan_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_Log_",fiber_subset,"_genef_df.RData",sep=""), envir = parent.frame())
load(paste("/home/slancast/NormAggreg_",fiber_subset,"_clinicals_df.RData",sep=""), envir = parent.frame())

#t(logaggregnorm_rna), t(logaggregnorm_pathabundance), t(logaggregnorm_metaphlan), t(logaggregnorm_genef)
combined_df <- rbind(t(logaggregnorm_cytokine), t(logaggregnorm_metabolomics),t(logaggregnorm_lipids),t(logaggregnorm_proteomics),t(aggregnorm_clinicals),t(logaggregnorm_rna), t(logaggregnorm_pathabundance), t(logaggregnorm_metaphlan), t(logaggregnorm_genef))
combined_df <- combined_df[-which(rownames(combined_df) %in% c("UNMAPPED")),] #Uninformative and duplicated so I got rid of them

if (nrow(logaggregnorm_cytokine) == 7) {
  xaxis_ticks = c("B","10","20","30","D3","D10","WF")
} else if (nrow(logaggregnorm_cytokine) == 6) {xaxis_ticks = c("B","10","20","30","D3","D10")}


print("combined")
class(combined_df) <- "numeric"
library(matrixStats)
combined_df <- combined_df[rowVars(combined_df)>0,]
set.seed(1)
library(Mfuzz)
eset <- ExpressionSet(combined_df) #Creating the type expression set with the metadata rows as a different argument
m <- exprs(eset)
eset <- standardise(eset) #Running standarise 
eset <- standardise2(eset) #Running standarise2 

m1 = mestimate(eset)

for (z in 2:16) {
mfuzzcl <- mfuzz(eset, c=z, m=m1)
print(z)

pdf(paste("/home/slancast/multiomics_clusters/clusters",z,"multiomics_",fiber_subset,"_parcoord.pdf",sep=""))
par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
mfuzz.plot2(eset,cl=mfuzzcl,ylim=c(-3,3),mfrow=c(4,4), time.labels = xaxis_ticks,bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
dev.off()

pdf(paste("/home/slancast/multiomics_clusters/clusters",z,"multiomics_",fiber_subset,"_parcoord2.pdf",sep=""))
par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
mfuzz.plot2(eset,cl=mfuzzcl,ylim=c(-3,3),mfrow=c(4,4), time.labels = xaxis_ticks,bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
dev.off()

cluster_membership <- data.frame(mfuzzcl$cluster)
cluster_membership$Gene <- rownames(cluster_membership)

write.table(cluster_membership, file = paste("/home/slancast/multiomics_clusters/clusters",z,"multiomics_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")

}

} #End fiber subset

#The DAVID function is after this:

system("sudo poweroff")


cluster_membership <- entrezid(cluster_membership)



library(RDAVIDWebService)
david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setTimeOut(david, 100000)
background <- addList(david, cluster_membership$EnsemblGene, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")


for (i in unique(cluster_membership$mfuzzcl.cluster)) {
  print(fiber_subset)
  print(z)
  print(i)
  assign(paste("cluster",as.character(i),sep=""), cluster_membership[ which(cluster_membership$mfuzzcl.cluster==i),])
  write.table(get(paste("cluster",as.character(i),sep="")), file = paste("/home/slancast/multiomics_clusters/clusters",z,"cluster",as.character(i),"_",fiber_subset,".txt",sep=""),sep="\t")
  matrix <- get(paste("cluster",as.character(i),sep=""))
  result<-addList(david, matrix$EnsemblGene, idType="ENSEMBL_GENE_ID",listName=paste("cluster",as.character(i),sep=""), listType="Gene")
  getFunctionalAnnotationChartFile(david, fileName = paste("/home/slancast/multiomics_clusters/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationChart.txt",sep=""))
  getFunctionalAnnotationTableFile(david, fileName = paste("/home/slancast/multiomics_clusters/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationTable.txt",sep=""))
  getClusterReportFile(david, fileName = paste("/home/slancast/multiomics_clusters/clusters",z,"cluster",as.character(i),"_",fiber_subset,"ClusterReport.txt",sep=""))
} 


