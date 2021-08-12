library(Mfuzz)
library(MASS)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  
  load(paste("~/NormAggreg_Log_",fiber_subset,"_genef_df.RData",sep=""))
  
  if (nrow(logaggregnorm_genef) == 7) {
    xaxis_ticks = c("B","10","20","30","D3","D10","WF")
  } else if (nrow(logaggregnorm_genef) == 6) {xaxis_ticks = c("B","10","20","30","D3","D10")}
  
  genef_df <- data.frame(t(logaggregnorm_genef)) #Must be a matrix where all the values are numeric, but for some fucking reason this dataset won't do that automatically.
  genef_df <- data.frame(lapply(genef_df, as.character), stringsAsFactors=FALSE, row.names = rownames(genef_df)) 
  genef_df <- data.frame(lapply(genef_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(genef_df)) 
  genef_df <- as.matrix(genef_df)
  
  genef_eset1 <- ExpressionSet(genef_df) #Creating the type expression set with the metadata rows as a different argument
  genef_eset1 <- standardise(genef_eset1) #Running standarise 
  genef_eset1 <- standardise2(genef_eset1) #Running standarise 
  genef_m <- exprs(genef_eset1)
  genef_m1 = mestimate(genef_eset1)
  
  dir.create(paste("~/genef/genef_clustering_",fiber_subset,sep=""))
  
  for (z in 2:16) {
    print("cluters number:")
    print(z)
    
    #genef_cluster_number <- Dmin(genef_eset1,genef_m1,crange=seq(2,15,1),repeats=5,visu=TRUE)
    genef_mfuzzcl <- mfuzz(genef_eset1, c=z, m=genef_m1)
    genef_standardised_data <- data.matrix(genef_m)
    genef_cluster_labels <- genef_mfuzzcl$cluster
    #parcoord(standardised_data, col = cluster_labels, var.label = TRUE)
    
    pdf(paste("~/genef/genef_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,"_parcoord.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(genef_eset1,ylim=c(-3,3),cl=genef_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
    dev.off()
    
    pdf(paste("~/genef/genef_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,"_parcoord2.pdf",sep=""))
    par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
        mar = c(2,1,1,1) + 0.1)
    mfuzz.plot2(genef_eset1,ylim=c(-3,3),cl=genef_mfuzzcl,mfrow=c(3,3), time.labels = xaxis_ticks, bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.lab=0.1)
    dev.off()
    
    cluster_membership <- data.frame(genef_mfuzzcl$cluster)
    write.table(cluster_membership, file = paste("~/genef/genef_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,".txt",sep=""), sep="\t")
    cluster_matrix <- cbind(genef_df, cluster_membership)
    write.table(cluster_membership, file = paste("~/genef/genef_clustering_",fiber_subset,"/cluster_membership_",z,fiber_subset,".txt",sep=""), sep="\t")
    
  } # End of cluster number loop
  
} # End of fiber loop

system("sudo poweroff")