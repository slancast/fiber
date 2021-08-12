library(Mfuzz)
library(MASS)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  
  for (i in c("kingdom","phylum","order","class","family","genus","species","strain")){
    print(i)
    
  
  load(paste("~/taxa/NormAggreg_Log_",fiber_subset,i,"_metaphlan_df.RData",sep=""))
  
  if (nrow(logaggregnorm_metaphlan) == 7) {
    xaxis_ticks = c("B","10","20","30","D3","D10","WF")
  } else if (nrow(logaggregnorm_metaphlan) == 6) {xaxis_ticks = c("B","10","20","30","D3","D10")}
  
  metaphlan_df <- data.frame(t(logaggregnorm_metaphlan)) #Must be a matrix where all the values are numeric, but for some fucking reason this dataset won't do that automatically.
  metaphlan_df <- data.frame(lapply(metaphlan_df, as.character), stringsAsFactors=FALSE, row.names = rownames(metaphlan_df)) 
  metaphlan_df <- data.frame(lapply(metaphlan_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(metaphlan_df)) 
  metaphlan_df <- as.matrix(metaphlan_df)
  
  metaphlan_eset1 <- ExpressionSet(metaphlan_df) #Creating the type expression set with the metadata rows as a different argument
  testing <- metaphlan_eset1 <- standardise(metaphlan_eset1) #Running standarise 
  metaphlan_eset1 <- standardise2(metaphlan_eset1) #Running standarise 
  metaphlan_m <- exprs(metaphlan_eset1)
  metaphlan_m1 = mestimate(metaphlan_eset1)
  
  dir.create(paste("~/metaphlan/metaphlan_clustering_",fiber_subset,i,sep=""))
  a <- tryCatch({ 
  for (z in 2:16) {
    print("cluters number:")
    print(z)
    
    #metaphlan_cluster_number <- Dmin(metaphlan_eset1,metaphlan_m1,crange=seq(2,15,1),repeats=5,visu=TRUE)
    metaphlan_mfuzzcl <- mfuzz(metaphlan_eset1, c=z, m=metaphlan_m1)
    library(MASS)
    metaphlan_standardised_data <- data.matrix(metaphlan_m)
    metaphlan_cluster_labels <- metaphlan_mfuzzcl$cluster
    #parcoord(standardised_data, col = cluster_labels, var.label = TRUE)
    
    pdf(paste("~/metaphlan/metaphlan_clustering_",fiber_subset,i,"/cluster_membership_",z,fiber_subset,i,"_parcoord.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(metaphlan_eset1,ylim=c(-3,3),cl=metaphlan_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
    dev.off()
    
    pdf(paste("~/metaphlan/metaphlan_clustering_",fiber_subset,i,"/cluster_membership_",z,fiber_subset,i,"_parcoord2.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(metaphlan_eset1,ylim=c(-3,3),cl=metaphlan_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
    dev.off()
    
    cluster_membership <- data.frame(metaphlan_mfuzzcl$cluster)
    write.table(cluster_membership, file = paste("~/metaphlan/metaphlan_clustering_",fiber_subset,i,"/cluster_membership_",z,fiber_subset,i,".txt",sep=""), sep="\t")
    cluster_matrix <- cbind(metaphlan_df, cluster_membership)
    write.table(cluster_membership, file = paste("~/metaphlan/metaphlan_clustering_",fiber_subset,i,"/cluster_membership_",z,fiber_subset,i,".txt",sep=""), sep="\t")
    
    for (j in unique(cluster_membership$metaphlan_mfuzzcl.cluster)) {
      testing_cluster <- cluster_matrix[which(cluster_membership$metaphlan_mfuzzcl.cluster==j),]
      write.table(testing_cluster, file = paste("~/metaphlan/metaphlan_clustering_",fiber_subset,i,"/cluster_membership_",z,fiber_subset,i,"cluster",j,".txt",sep=""), sep="\t")
    } # End of writing cluster membership
  }
  },
   error  = function(err){print(err)
     eset = c()})
  
  } # End taxa loop
  
} # End of fiber loop
  
  
