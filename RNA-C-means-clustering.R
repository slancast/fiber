library(Mfuzz)
library(MASS)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  
  load(paste("~/NormAggreg_Log_",fiber_subset,"_RNA_df.RData",sep=""))
  
  if (nrow(logaggregnorm_rna) == 7) {
    xaxis_ticks = c("B","10","20","30","D3","D10","WF")
  } else if (nrow(logaggregnorm_rna) == 6) {xaxis_ticks = c("B","10","20","30","D3","D10")}
  
  RNA_df <- data.frame(t(logaggregnorm_rna)) #Must be a matrix where all the values are numeric, but for some fucking reason this dataset won't do that automatically.
  RNA_df <- data.frame(lapply(RNA_df, as.character), stringsAsFactors=FALSE, row.names = rownames(RNA_df)) 
  RNA_df <- data.frame(lapply(RNA_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(RNA_df)) 
  RNA_df <- as.matrix(RNA_df)
  
  RNA_eset1 <- ExpressionSet(RNA_df) #Creating the type expression set with the metadata rows as a different argument
  RNA_eset1 <- standardise(RNA_eset1) #Running standarise 
  RNA_eset1 <- standardise2(RNA_eset1) #Running standarise 
  RNA_m <- exprs(RNA_eset1)
  RNA_m1 = mestimate(RNA_eset1)
  
  dir.create(paste("~/RNA/RNA_clustering_",fiber_subset,sep=""))
  
  for (z in 2:16) {
    print("cluters number:")
    print(z)
    
    RNA_cluster_number <- Dmin(RNA_eset1,RNA_m1,crange=seq(2,15,1),repeats=5,visu=TRUE)
    RNA_mfuzzcl <- mfuzz(RNA_eset1, c=z, m=RNA_m1)
    library(MASS)
    RNA_standardised_data <- data.matrix(RNA_m)
    RNA_cluster_labels <- RNA_mfuzzcl$cluster
    #parcoord(standardised_data, col = cluster_labels, var.label = TRUE)
    
    pdf(paste("~/RNA/RNA_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,"_parcoord.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(RNA_eset1,ylim=c(-3,3),cl=RNA_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
    dev.off()
    
    pdf(paste("~/RNA/RNA_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,"_parcoord2.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(RNA_eset1,ylim=c(-3,3),cl=RNA_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
    dev.off()
    
    cluster_membership <- data.frame(RNA_mfuzzcl$cluster)
    write.table(cluster_membership, file = paste("~/RNA/RNA_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,".txt",sep=""), sep="\t")
    cluster_matrix <- cbind(RNA_df, cluster_membership)
    write.table(cluster_membership, file = paste("~/RNA/RNA_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,".txt",sep=""), sep="\t")
    
    for (i in unique(cluster_membership$RNA_mfuzzcl.cluster)) {
      testing_cluster <- cluster_matrix[which(cluster_membership$RNA_mfuzzcl.cluster==i),]
      write.table(testing_cluster, file = paste("~/RNA/RNA_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,"cluster",i,".txt",sep=""), sep="\t")
    } # End of writing cluster membership
    
  } # End of cluster number loop
  
} # End of fiber loop
