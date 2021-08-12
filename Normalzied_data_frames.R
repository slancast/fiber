#This will take the baseline means and subtract them from the data files in order to find the adjusted 
#values.
#

for (fiber_subset in c("Arabinoxylan","LCInulin", "Mix")) {
  print(fiber_subset)
if (fiber_subset == "Mix") {
  order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
  order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","Washout_Final")
}

########################################
#############Cytokine###################
########################################
if (FALSE) {
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/cytokine_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_cytokine_df.RData",sep=""))

cytokine_metadata <- data.frame(t(cytokine_metadata))
cytokine_metadata <- cytokine_metadata[which(cytokine_metadata$Fiber==fiber_subset),] #subsetting by fiber

participant <- cytokine_metadata$Participant
cytokine_df2 <- cbind(participant,cytokine_df2)

if (length(colnames(cytokine_df2[2:length(cytokine_df2)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
  print("Warning column names are not the same size")
}

baseline_mean_row <- apply(cytokine_df2, 1, function(x) match(x[1],baseline_means[,1]))

cytokine_df2 <- data.frame(lapply(cytokine_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(cytokine_df2)) #To change to a double, this needs to go through character first
cytokine_df2 <- data.frame(lapply(cytokine_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(cytokine_df2)) #then numeric

normalized_cytokine_df <- c()
for (i in 1:length(baseline_mean_row)){
  normalized_row <- as.numeric(cytokine_df2[i,2:ncol(cytokine_df2)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
  normalized_cytokine_df <- rbind(normalized_cytokine_df, normalized_row)
}

rownames(normalized_cytokine_df) <- rownames(cytokine_df2)
colnames(normalized_cytokine_df) <- colnames(cytokine_df2[2:ncol(cytokine_df2)])

save(normalized_cytokine_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_cytokine_df.RData",sep=""))

###################### log_transorming ##########################
#AggregNorm#

metadata <- data.frame(t(cytokine_metadata))
datamatrix3 <- as.data.frame(metadata["Week",])
datamatrix3 <- as.data.frame(t(datamatrix3))
datamatrix3 <- cbind(datamatrix3, normalized_cytokine_df)
datamatrix3 <- aggregate(datamatrix3,list(Visit=datamatrix3$Week), mean)
datamatrix3 <- datamatrix3[ , -2]
datamatrix3 <- datamatrix3[match(order, datamatrix3$Visit),]
datamatrix3 <- na.omit(datamatrix3)
rownames(datamatrix3) <- datamatrix3$Visit
datamatrix3 <- datamatrix3[,-1]
aggregnorm_cytokine <- datamatrix3

save(aggregnorm_cytokine, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_cytokine_df.RData",sep=""))
}

  
  ##########################################
  #################LIPIDS###################
  ##########################################
  
  load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/lipids_baselines.RData")
  load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_lipids_df.RData",sep=""))
  
  lipids_metadata <- data.frame(t(lipids_metadata))
  participant <- lipids_metadata$subject_id
  lipids_df2 <- cbind(participant,lipids_df2)
  
  if (length(colnames(lipids_df2[,2:ncol(lipids_df2)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
    print("Warning column names are not the same size")
  }
  
  baseline_mean_row <- apply(lipids_df2, 1, function(x) match(x[1],baseline_means[,1]))
  
  lipids_df2 <- data.frame(lapply(lipids_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(lipids_df2)) #To change to a double, this needs to go through character first
  lipids_df2 <- data.frame(lapply(lipids_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(lipids_df2)) #then numeric
  
  normalized_lipids_df <- c()
  for (i in 1:length(baseline_mean_row)){
    normalized_row <- as.numeric(lipids_df2[i,2:ncol(lipids_df2)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
    normalized_lipids_df <- rbind(normalized_lipids_df, normalized_row)
  }
  
  rownames(normalized_lipids_df) <- rownames(lipids_df2)
  colnames(normalized_lipids_df) <- colnames(lipids_df2[2:ncol(lipids_df2)])
  
  save(normalized_lipids_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_lipids_df.RData",sep=""))
  
  ###################### log_transorming ##########################
  #AggregNorm#
  
  lipids_metadata <- data.frame(t(lipids_metadata))
  lipids_df3 <- as.data.frame(lipids_metadata["timepoint",])
  lipids_df3 <- as.data.frame(t(lipids_df3))
  lipids_df3 <- cbind(lipids_df3, normalized_lipids_df)
  lipids_df3 <- aggregate(lipids_df3,list(Visit=lipids_df3$timepoint), mean)
  lipids_df3 <- lipids_df3[ , -2]
  lipids_df3 <- lipids_df3[match(order, lipids_df3$Visit),]
  lipids_df3 <- na.omit(lipids_df3)
  rownames(lipids_df3) <- lipids_df3$Visit
  lipids_df3 <- lipids_df3[,-1]
  aggregnorm_lipids <- lipids_df3
  
  save(aggregnorm_lipids, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_lipids_df.RData",sep=""))

  
#######################################
#################RNA###################
#######################################
  
if (FALSE) {

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/rna_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_rna_df.RData",sep=""))

rna_metadata <- data.frame(t(rna_metadata))
participant <- rna_metadata$participant
rna_df2 <- cbind(participant,rna_df2)

if (length(colnames(rna_df2[,2:ncol(rna_df2)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
  print("Warning column names are not the same size")
}

baseline_mean_row <- apply(rna_df2, 1, function(x) match(x[1],baseline_means[,1]))

rna_df2 <- data.frame(lapply(rna_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(rna_df2)) #To change to a double, this needs to go through character first
rna_df2 <- data.frame(lapply(rna_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(rna_df2)) #then numeric

normalized_rna_df <- c()
for (i in 1:length(baseline_mean_row)){
  normalized_row <- as.numeric(rna_df2[i,2:ncol(rna_df2)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
  normalized_rna_df <- rbind(normalized_rna_df, normalized_row)
}

rownames(normalized_rna_df) <- rownames(rna_df2)
colnames(normalized_rna_df) <- colnames(rna_df2[2:ncol(rna_df2)])

save(normalized_rna_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_rna_df.RData",sep=""))

###################### log_transorming ##########################
#AggregNorm#

rna_metadata <- data.frame(t(rna_metadata))
rna_df3 <- as.data.frame(rna_metadata["week",])
rna_df3 <- as.data.frame(t(rna_df3))
rna_df3 <- cbind(rna_df3, normalized_rna_df)
rna_df3 <- aggregate(rna_df3,list(Visit=rna_df3$week), mean)
rna_df3 <- rna_df3[ , -2]
rna_df3 <- rna_df3[match(order, rna_df3$Visit),]
rna_df3 <- na.omit(rna_df3)
rownames(rna_df3) <- rna_df3$Visit
rna_df3 <- rna_df3[,-1]
aggregnorm_rna <- rna_df3

save(aggregnorm_rna, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_RNA_df.RData",sep=""))
}

#######################################
#################PCL###################
#######################################
  if (FALSE) {
#I think these are different because there are different bugs list files that I used
#the baseline ones are with just the genus, and the pcl_df has the entire bugs list.
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/pcl_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_pcl_df.RData",sep=""))

pcl_metadata <- data.frame(t(pcl_metadata))
pcl_metadata <- pcl_metadata[which(pcl_metadata$Fiber==fiber_subset),] #subsetting by fiber

participant <- pcl_metadata$Participant
pcl_df2 <- data.frame(cbind(as.character(participant),pcl_df2))

if (length(colnames(pcl_df2[,2:ncol(pcl_df2)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
  print("Warning column names are not the same size")
}

baseline_mean_row <- apply(pcl_df2, 1, function(x) match(x[1],baseline_means[,1]))

pcl_df2 <- data.frame(lapply(pcl_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(pcl_df2)) #To change to a double, this needs to go through character first
pcl_df2 <- data.frame(lapply(pcl_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(pcl_df2)) #then numeric

normalized_pcl_df <- c()
for (i in 1:length(baseline_mean_row)){
  normalized_row <- as.numeric(pcl_df2[i,2:ncol(pcl_df2)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
  normalized_pcl_df <- rbind(normalized_pcl_df, normalized_row)
}

rownames(normalized_pcl_df) <- rownames(pcl_df2)
colnames(normalized_pcl_df) <- colnames(pcl_df2[2:ncol(pcl_df2)])

save(normalized_pcl_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_pcl_df.RData",sep=""))

###################### log_transorming ##########################
#AggregNorm#

pcl_metadata <- data.frame(t(pcl_metadata))
pcl_df3 <- as.data.frame(pcl_metadata["Dose",])
pcl_df3 <- as.data.frame(t(pcl_df3))
pcl_df3 <- cbind(pcl_df3, normalized_pcl_df)
pcl_df3 <- aggregate(pcl_df3,list(Visit=pcl_df3$Dose), mean)
pcl_df3 <- pcl_df3[ , -2]
pcl_df3 <- pcl_df3[match(order, pcl_df3$Visit),]
pcl_df3 <- na.omit(pcl_df3)
rownames(pcl_df3) <- pcl_df3$Visit
pcl_df3 <- pcl_df3[,-1]
aggregnorm_pcl <- pcl_df3

save(aggregnorm_pcl, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_pcl_df.RData",sep=""))
  }
  
#############################################
#################Metabolome##################
#############################################
if (FALSE) {

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/metabolomic_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_metabolomics_df.RData",sep=""))

metabolomics_metadata <- data.frame(t(metabolomics_metadata))
metabolomics_metadata <- metabolomics_metadata[which(metabolomics_metadata$fiber==fiber_subset),] #subsetting by fiber

participant <- metabolomics_metadata$subject_id
metabolomics_df2 <- cbind(participant,metabolomics_df2)

if (length(colnames(metabolomics_df2[,2:ncol(metabolomics_df2)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
  print("Warning column names are not the same size")
}

baseline_mean_row <- apply(metabolomics_df2, 1, function(x) match(as.character(x[1]),trimws(as.character(baseline_means[,1]))))

metabolomics_df2 <- data.frame(lapply(metabolomics_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df2)) #To change to a double, this needs to go through character first
metabolomics_df2 <- data.frame(lapply(metabolomics_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df2)) #then numeric

normalized_metabolomics_df <- c()
for (i in 1:length(baseline_mean_row)){
  normalized_row <- as.numeric(metabolomics_df2[i,2:ncol(metabolomics_df2)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
  normalized_metabolomics_df <- rbind(normalized_metabolomics_df, normalized_row)
}

rownames(normalized_metabolomics_df) <- rownames(metabolomics_df2)
colnames(normalized_metabolomics_df) <- colnames(metabolomics_df2[2:ncol(metabolomics_df2)])

save(normalized_metabolomics_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_metabolome_df.RData",sep=""))

###################### log_transorming ##########################
#AggregNorm#

metabolomics_metadata <- data.frame(t(metabolomics_metadata))
metabolomics_df3 <- as.data.frame(metabolomics_metadata["timepoint",])
metabolomics_df3 <- as.data.frame(t(metabolomics_df3))
metabolomics_df3 <- cbind(metabolomics_df3, normalized_metabolomics_df)
metabolomics_df3 <- aggregate(metabolomics_df3,list(Visit=metabolomics_df3$timepoint), mean)
metabolomics_df3 <- metabolomics_df3[ , -2]
metabolomics_df3 <- metabolomics_df3[match(order, metabolomics_df3$Visit),]
rownames(metabolomics_df3) <- metabolomics_df3$Visit
metabolomics_df3  <- metabolomics_df3[,-2]
metabolomics_df3 <- t(metabolomics_df3)
metabolomics_df3 <- na.omit(metabolomics_df3)
metabolomics_df3 <- t(metabolomics_df3)
metabolomics_df3  <- metabolomics_df3[,-1]
aggregnorm_metabolomics <- metabolomics_df3

save(aggregnorm_metabolomics, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_metabolomics_df.RData",sep=""))
}

}

