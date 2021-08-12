# I want some personalized responses to fiber supplementation
# this will run the multiomics dataset subseting by individual
#First attach participant number
#Then add in missing values in all potential timepoints
#Then merge the datasets
#Lastly run the clustering.
#
#For some reason the cluster had problems loading the LCInulin genefamiles file the
#The error was: "Error: error reading from connection"
#Online solutions were to change permissions on the file, didn't solve anything
#Or that it mysteriously was fixed after reboot, didn't solve either
#It loads fine on the laptop, so I will try to run the LCInulin on the laptop
#
#

source("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Software/Utils.r")


for (participant in c("1","23","28","53","63","64","69","74","100","1005","1008","1010","1015","104","107","111","114","123")){
print("Participant:")
print(participant)
time_points <- c("Baseline","10","20","30","WashoutD3","WashoutD10")

for (fiber_subset in c("LCInulin")) {

print(participant)
print(paste("fiber_subset: ",fiber_subset,sep="") )

#Loading in the datasets. The normalized dataset does not contain the metadata
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_clinical_df.RData",sep=""), envir = parent.frame())
  
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_cytokine_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_pcl_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_rna_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_metabolome_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_genef_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_lipids_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_metaphlan_df.RData",sep=""), envir = parent.frame())

load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_cytokine_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_rna_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metabolomics_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_genef_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_lipids_df.RData",sep=""), envir = parent.frame())
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metaphlan_df.RData",sep=""), envir = parent.frame())

#Appending relevant metadata to the normalized dataframes
cytokine_df <-  cbind(t(cytokine_metadata["Participant",]),t(cytokine_metadata["Week",]),normalized_cytokine_df)
pcl_df <- cbind(t(pcl_metadata["Participant",]),t(pcl_metadata["Dose",]),normalized_pcl_df)
rna_df <-  cbind(t(rna_metadata["participant",]),t(rna_metadata["week",]),normalized_rna_df)
metabolomics_df <- cbind(t(metabolomics_metadata["subject_id",]),t(metabolomics_metadata["timepoint",]),normalized_metabolomics_df)
genef_df <- cbind(t(genef_metadata["Participant",]),t(genef_metadata["Dose",]),normalized_genef_df)
metaphlan_df <- cbind(t(metaphlan_metadata["Participant",]),t(metaphlan_metadata["Dose",]),normalized_metaphlan_df)
lipids_df <- cbind(t(lipids_metadata["subject_id",]),t(lipids_metadata["timepoint",]),normalized_lipids_df)

clinical_metadata <- data.frame(t(clinical_metadata))
clinical_df <- cbind(clinical_metadata$participant,clinical_metadata$week,clinical_df2)

#Subsetting by the participant and lining up the correct time points
print("Clinicals")
clinical_df <- data.frame(clinical_df)
participants <- gsub(" ","",as.character(clinical_df$clinical_metadata.participant)) #There were extra spaces in the cytokine participant labels.
clinical_df2 <- clinical_df[which(participants==participant),] 
clinical_df2 <- clinical_df2[which(clinical_df2$clinical_metadata.week %in% time_points),] 
missing_samples <- is.na(match(time_points, clinical_df2$clinical_metadata.week))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(clinical_df2)-2))
    clinical_df2 <- rbind(clinical_df2,new_row) 
  }#Ending if statement
}#Ending loop
clinical_df2 <- clinical_df2[match(time_points, clinical_df2$clinical_metadata.week),]
clinical_df2 <- data.frame(t(clinical_df2))
colnames(clinical_df2) <- time_points

#Subsetting by the participant and lining up the correct time points
print("cytokines")
cytokine_df <- data.frame(cytokine_df)
participants <- gsub(" ","",as.character(cytokine_df$Participant)) #There were extra spaces in the cytokine participant labels.
cytokine_df2 <- cytokine_df[which(participants==participant),] 
cytokine_df2 <- cytokine_df2[which(cytokine_df2$Week %in% time_points),] 
missing_samples <- is.na(match(time_points, cytokine_df2$Week))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(cytokine_df2)-2))
    cytokine_df2 <- rbind(cytokine_df2,new_row) 
  }#Ending if statement
}#Ending loop
cytokine_df2 <- cytokine_df2[match(time_points, cytokine_df2$Week),]
cytokine_df2 <- data.frame(t(cytokine_df2))
colnames(cytokine_df2) <- time_points

print("pcl")
pcl_df <- data.frame(pcl_df)
pcl_df2 <- pcl_df[which(as.character(pcl_df$Participant)==participant),] 
pcl_df2 <- pcl_df2[which(pcl_df2$Dose %in% time_points),] 
missing_samples <- is.na(match(time_points, pcl_df2$Dose))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
  new_row <- c(participant,time_points[counter],rep(NA,ncol(pcl_df2)-2))
  pcl_df2 <- rbind(pcl_df2,new_row) 
  }#Ending if statement
}#Ending loop
pcl_df2 <- pcl_df2[match(time_points, pcl_df2$Dose),]
pcl_df2 <- data.frame(t(pcl_df2))
colnames(pcl_df2) <- time_points

print("metaphlan")
metaphlan_df <- data.frame(metaphlan_df)
metaphlan_df2 <- metaphlan_df[which(as.character(metaphlan_df$Participant)==participant),] 
metaphlan_df2 <- metaphlan_df2[which(metaphlan_df2$Dose %in% time_points),] 
missing_samples <- is.na(match(time_points, metaphlan_df2$Dose))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(metaphlan_df2)-2))
    metaphlan_df2 <- rbind(metaphlan_df2,new_row) 
  }#Ending if statement
}#Ending loop
metaphlan_df2 <- metaphlan_df2[match(time_points, metaphlan_df2$Dose),]
metaphlan_df2 <- data.frame(t(metaphlan_df2))
colnames(metaphlan_df2) <- time_points

print("lipids")
lipids_df <- data.frame(lipids_df)
lipids_df2 <- lipids_df[which(as.character(lipids_df$subject_id)==participant),] 
lipids_df2 <- lipids_df2[which(lipids_df2$timepoint %in% time_points),] 
missing_samples <- is.na(match(time_points, lipids_df2$timepoint))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(lipids_df2)-2))
    lipids_df2 <- rbind(lipids_df2,new_row) 
  }#Ending if statement
}#Ending loop
lipids_df2 <- lipids_df2[match(time_points, lipids_df2$subject_id),]
lipids_df2 <- data.frame(t(lipids_df2))
colnames(lipids_df2) <- time_points

print("rna")
rna_df <- data.frame(rna_df)
rna_df2 <- rna_df[which(as.character(rna_df$participant)==participant),] 
rna_df2 <- rna_df2[which(rna_df2$week %in% time_points),] 
missing_samples <- is.na(match(time_points, rna_df2$week))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(rna_df2)-2))
    rna_df2 <- rbind(rna_df2,new_row) 
  }#Ending if statement
}#Ending loop
rna_df2 <- rna_df2[match(time_points, rna_df2$week),]
rna_df2 <- data.frame(t(rna_df2))
colnames(rna_df2) <- time_points

print("genef")
genef_df <- data.frame(genef_df)
genef_df2 <- genef_df[which(as.character(genef_df$Participant)==participant),] 
genef_df2 <- genef_df2[which(genef_df2$Dose %in% time_points),] 
missing_samples <- is.na(match(time_points, genef_df2$Dose))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(genef_df2)-2))
    genef_df2 <- rbind(genef_df2,new_row) 
  }#Ending if statement
}#Ending loop
genef_df2 <- genef_df2[match(time_points, genef_df2$Dose),]
genef_df2 <- data.frame(t(genef_df2))
colnames(genef_df2) <- time_points

print("metabolomcis")
metabolomics_df <- data.frame(metabolomics_df)
metabolomics_df2 <- metabolomics_df[which(metabolomics_df$subject_id==participant),] 
metabolomics_df2 <- metabolomics_df2[which(metabolomics_df2$timepoint %in% time_points),] 
missing_samples <- is.na(match(time_points, metabolomics_df2$timepoint))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(metabolomics_df2)-2))
    metabolomics_df2 <- rbind(metabolomics_df2,new_row) 
  }#Ending if statement
}#Ending loop
metabolomics_df2 <- metabolomics_df2[match(time_points, metabolomics_df2$timepoint),]
metabolomics_df2 <- data.frame(t(metabolomics_df2))
colnames(metabolomics_df2) <- time_points

metabolomics_df <- data.frame(metabolomics_df)
metabolomics_df2 <- metabolomics_df[which(metabolomics_df$subject_id==participant),] 
metabolomics_df2 <- metabolomics_df2[which(metabolomics_df2$timepoint %in% time_points),] 
missing_samples <- is.na(match(time_points, metabolomics_df2$timepoint))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(metabolomics_df2)-2))
    metabolomics_df2 <- rbind(metabolomics_df2,new_row) 
  }#Ending if statement
}#Ending loop
metabolomics_df2 <- metabolomics_df2[match(time_points, metabolomics_df2$timepoint),]
metabolomics_df2 <- data.frame(t(metabolomics_df2))
colnames(metabolomics_df2) <- time_points

metabolomics_df <- data.frame(metabolomics_df)
metabolomics_df2 <- metabolomics_df[which(metabolomics_df$subject_id==participant),] 
metabolomics_df2 <- metabolomics_df2[which(metabolomics_df2$timepoint %in% time_points),] 
missing_samples <- is.na(match(time_points, metabolomics_df2$timepoint))
counter = 0 
for (i in missing_samples) {
  counter = counter + 1
  print(i)
  if (isTRUE(i)) {
    new_row <- c(participant,time_points[counter],rep(NA,ncol(metabolomics_df2)-2))
    metabolomics_df2 <- rbind(metabolomics_df2,new_row) 
  }#Ending if statement
}#Ending loop
metabolomics_df2 <- metabolomics_df2[match(time_points, metabolomics_df2$timepoint),]
metabolomics_df2 <- data.frame(t(metabolomics_df2))
colnames(metabolomics_df2) <- time_points

combined_df <- rbind(clinical_df2[3:nrow(clinical_df2),], cytokine_df2[3:nrow(cytokine_df2),], metaphlan_df2[3:nrow(metaphlan_df2),], pcl_df2[3:nrow(pcl_df2),], genef_df2[3:nrow(genef_df2),], rna_df2[3:nrow(rna_df2),], metabolomics_df2[3:nrow(metabolomics_df2),], genef_df2[3:nrow(genef_df2),])
combined_df <- na.omit(combined_df)

print("combined")
combined_df <- as.matrix(combined_df)
class(combined_df) <- "numeric"
set.seed(1)
library(Mfuzz)
eset <- ExpressionSet(combined_df) #Creating the type expression set with the metadata rows as a different argument
testing <- tryCatch({ eset <- standardise(eset) 
m <- exprs(eset) 
# m <- na.omit(m) #Somehow running standardise introduces NAs into the matrix
# eset <- ExpressionSet(m) #Recreating the standardised expression set without the NAs
m1 = mestimate(eset)

for (z in 2:16) {
print("cluters number:")
print(z)
mfuzzcl <- mfuzz(eset, c=z, m=m1)

library(icesTAF)
mkdir(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,sep=""))

pdf(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,"/clusters",z,"personalized_",fiber_subset,"_parcoord2.pdf",sep=""))
plotting_fuzzy_clusters(eset, mfuzzcl)
dev.off()

cluster_membership <- data.frame(mfuzzcl$cluster)
cluster_membership$Gene <- rownames(cluster_membership)
#cluster_membership <- entrezid(cluster_membership)

write.table(cluster_membership, file = paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,"/clusters",z,"personalized_",fiber_subset,"_cluster_membership.txt",sep=""), sep="\t")



# library(RDAVIDWebService)
# david<-DAVIDWebService(email="slancast@stanford.edu", url="https:/david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
# setTimeOut(david, 1000000)
# background <- addList(david, cluster_membership$EnsemblGene, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")
# 
# for (i in unique(cluster_membership$mfuzzcl.cluster)) {
#   print(fiber_subset)
#   print(i)
#   assign(paste("cluster",as.character(i),sep=""), cluster_membership[ which(cluster_membership$mfuzzcl.cluster==i),])
#   write.table(get(paste("cluster",as.character(i),sep="")), file = paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,"/clusters",z,"cluster",as.character(i),"_",fiber_subset,".txt",sep=""),sep="\t")
#   matrix <- get(paste("cluster",as.character(i),sep=""))
#   result<-addList(david, matrix$EnsemblGene, idType="ENSEMBL_GENE_ID",listName=paste("cluster",as.character(i),sep=""), listType="Gene")
#   getFunctionalAnnotationChartFile(david, fileName = paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,"/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationChart.txt",sep=""))
#   getFunctionalAnnotationTableFile(david, fileName = paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,"/clusters",z,"cluster",as.character(i),"_",fiber_subset,"AnnotationTable.txt",sep=""))
#   getClusterReportFile(david, fileName = paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/",participant,fiber_subset,"/clusters",z,"cluster",as.character(i),"_",fiber_subset,"ClusterReport.txt",sep=""))
# } } 
} #ending cluster loop
},
error  = function(err){print(err)
  eset = c()})#Running standarise 
  
} #Ending fiber subset loop

} #Ending the personalized loop
