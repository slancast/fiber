#In the overfeeder paper, they standardized it using mfuzz::standarsize
#This requires an object of type ExpressionSet. This is a data structure that contains
#the matrix of interest, and then various metadata (arranged by row, column, etc)
library(Mfuzz)
#The data should be in the standard data format. 
pcl_df = read.csv("/Users/SLancaster/Desktop/relab-pathabundance.pcl", sep="\t", header=TRUE, row.names = 1)
bugs_list_df = read.csv("/Users/SLancaster/Desktop/bugs_list.pcl", sep="\t", header=TRUE, row.names = 1)

pcl_metadata_rows = 64
pcl_df <- data.frame(t(pcl_df))
fiber_subset = "Arabinoxylan"
pcl_df <- pcl_df[which(pcl_df$Fiber==fiber_subset),] #subsetting by fiber
pcl_df <- data.frame(t(pcl_df))
pcl_metadata = head(pcl_df,pcl_metadata_rows)
bugs_list_df <- data.frame(t(bugs_list_df))
bugs_list_df<- bugs_list_df[which(bugs_list_df$Fiber==fiber_subset),] #subsetting by fiber
bugs_list_df <- data.frame(t(bugs_list_df))
bugs_list_df2 <- tail(bugs_list_df, -pcl_metadata_rows) #Now getting rid of the metadata
bugs_list_df2 <- as.matrix(bugs_list_df2)
class(bugs_list_df2) <- "numeric"
bugs_list_df2 <- bugs_list_df2/100

pcl_df2 <- tail(pcl_df, -pcl_metadata_rows) #Now getting rid of the metadata
colnames(bugs_list_df2) <- colnames(pcl_df2)
pcl_df2 <- rbind(bugs_list_df2,pcl_df2)
pcl_df2 <- as.matrix(t(pcl_df2))
class(pcl_df2) <- "numeric"

pcl_df3 <- as.data.frame(pcl_metadata["Dose",])
pcl_df3 <- as.data.frame(t(pcl_df3))
pcl_df3 <- cbind(pcl_df3, pcl_df2)
pcl_df3 <- aggregate(pcl_df3,list(Visit=pcl_df3$Dose), mean)
pcl_df3 <- pcl_df3[ , -2]
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
pcl_df3 <- pcl_df3[match(order, pcl_df3$Visit),]
rownames(pcl_df3) <- pcl_df3$Visit

#Rows with 0 variance must be omitted to perform 
#This may not be necessary. Omitting rows with 0 sum sovles this problem
#that I was trying to address with this. I'll keep this funciton in case though.
removing_variance <- function(df3) {
  require(caret)
  pcl_df3 <- t(pcl_df3)
  pcl_df4 <- nearZeroVar(pcl_df3[,2:ncol(pcl_df3)], saveMetrics = TRUE)
  pcl_df3 <- t(pcl_df3)
  pcl_df5 <- cbind(pcl_df3,zeroVar=pcl_df4$zeroVar,nzv=pcl_df4$nzv)
  pcl_df5 <- data.frame(pcl_df5)
  pcl_df5 <- pcl_df5[which(pcl_df5$zeroVar=="FALSE"),] 
  pcl_df5 <- pcl_df5[which(pcl_df5$nzv=="FALSE"),] 
  pcl_df5 <- subset(pcl_df5, select=-c(zeroVar,nzv))
}

pcl_df5 <- pcl_df3
pcl_df5 <- as.matrix(t(pcl_df5[,-1]))#Needs to be a matrix to be an expressionset, and cytokines need to be as rows for standardization
class(pcl_df5) <- "numeric"
pcl_df5 <- pcl_df5[rowMeans(pcl_df5)>0,] #If all the rows are 0
save(pcl_metadata, pcl_df5, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/data/",fiber_subset,"_pcl_df.RData",sep=""))

pcl_eset1 <- ExpressionSet(na.omit(pcl_df5)) #Creating the type expression set with the metadata rows as a different argument
pcl_eset1 <- standardise(pcl_eset1) #Running standarise 
pcl_m <- exprs(pcl_eset1)
pcl_m1 = mestimate(pcl_eset1)

pcl_cluster_number <- Dmin(pcl_eset1,pcl_m1,crange=seq(2,15,1),repeats=5,visu=TRUE)
pcl_mfuzzcl <- mfuzz(pcl_eset1, c=10, m=pcl_m1)
library(MASS)
pcl_standardised_data <- data.matrix(pcl_m)
pcl_cluster_labels <- pcl_mfuzzcl$cluster
#parcoord(standardised_data, col = cluster_labels, var.label = TRUE)
dir.create(paste("/Users/SLancaster/Desktop/metagenome_clustering_",fiber_subset,sep=""))

pdf(paste("/Users/SLancaster/Desktop/metagenome_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,"_parcoord.pdf",sep=""))
par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
mfuzz.plot2(pcl_eset1,cl=pcl_mfuzzcl,mfrow=c(4,4), bg="black",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
dev.off()

pdf(paste("/Users/SLancaster/Desktop/metagenome_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,"_parcoord2.pdf",sep=""))
par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
mfuzz.plot2(pcl_eset1,cl=pcl_mfuzzcl,mfrow=c(4,4), bg="white",ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
dev.off()

cluster_membership <- data.frame(pcl_mfuzzcl$cluster)
write.table(cluster_membership, file = paste("/Users/SLancaster/Desktop/metagenome_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,".txt",sep=""), sep="\t")
cluster_matrix <- cbind(pcl_df5, cluster_membership)
write.table(cluster_membership, file = paste("/Users/SLancaster/Desktop/metagenome_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,".txt",sep=""), sep="\t")

for (i in unique(cluster_membership$pcl_mfuzzcl.cluster)) {
  testing_cluster <- cluster_matrix[which(cluster_membership$pcl_mfuzzcl.cluster==i),]
  write.table(testing_cluster, file = paste("/Users/SLancaster/Desktop/metagenome_clustering_",fiber_subset,"/cluster_membership_",fiber_subset,"cluster",i,".txt",sep=""), sep="\t")
  }

