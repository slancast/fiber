library(Mfuzz)
library(MASS)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  
  load(paste("~/NormAggreg_",fiber_subset,"_clinicals_df.RData",sep=""))
  #load(paste("~/Aggregate_",fiber_subset,"_clinical_df.RData",sep=""))
  
  if (nrow(aggregnorm_clinicals) == 7) {
    xaxis_ticks = c("B","10","20","30","D3","D10","WF")
  } else if (nrow(aggregnorm_clinicals) == 6) {xaxis_ticks = c("B","10","20","30","D3","D10")}
  
  clinical_df <- data.frame(t(aggregnorm_clinicals)) #Must be a matrix where all the values are numeric, but for some fucking reason this dataset won't do that automatically.
  clinical_df <- data.frame(lapply(clinical_df, as.character), stringsAsFactors=FALSE, row.names = rownames(clinical_df)) 
  clinical_df <- data.frame(lapply(clinical_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(clinical_df)) 
  clinical_df <- as.matrix(clinical_df)
  
  clinical_eset1 <- ExpressionSet(clinical_df) #Creating the type expression set with the metadata rows as a different argument
  clinical_eset1 <- standardise(clinical_eset1) #Running standarise 
  clinical_eset1 <- standardise2(clinical_eset1) #Running standarise 
  clinical_m <- exprs(clinical_eset1)
  clinical_m1 = mestimate(clinical_eset1)
  
  clinical_cluster_number <- Dmin(clinical_eset1,clinical_m1,crange=seq(2,15,1),repeats=5,visu=TRUE)
  clinical_mfuzzcl <- mfuzz(clinical_eset1, c=9, m=clinical_m1)
  clinical_standardised_data <- data.matrix(clinical_m)
  clinical_cluster_labels <- clinical_mfuzzcl$cluster

  dir.create(paste("~/clinical/clinical_clustering_",fiber_subset,sep=""))
  
  pdf(paste("~/clinical/clinical_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,"_parcoord.pdf",sep=""))
  par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
      mar = c(2,1,1,1) + 0.1)
  mfuzz.plot2(clinical_eset1,ylim=c(-3,3),cl=clinical_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
  dev.off()
  
  pdf(paste("~/clinical/clinical_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,"_parcoord2.pdf",sep=""))
  par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
      mar = c(2,1,1,1) + 0.1)
  mfuzz.plot2(clinical_eset1,ylim=c(-3,3),cl=clinical_mfuzzcl,mfrow=c(4,4),  bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
  dev.off()
  
  cluster_membership <- data.frame(clinical_mfuzzcl$cluster)
  write.table(cluster_membership, file = paste("~/clinical/clinical_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,".txt",sep=""), sep="\t")
  cluster_matrix <- cbind(clinical_df, cluster_membership)
  write.table(cluster_membership, file = paste("~/clinical/clinical_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,".txt",sep=""), sep="\t")
  
  for (i in unique(cluster_membership$clinical_mfuzzcl.cluster)) {
    testing_cluster <- cluster_matrix[which(cluster_membership$clinical_mfuzzcl.cluster==i),]
    write.table(testing_cluster, file = paste("~/clinical/clinical_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,"cluster",i,".txt",sep=""), sep="\t")
  }
  
}
