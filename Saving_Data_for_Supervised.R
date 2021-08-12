#Creating the metabolomics supervised clustering files
#

#################################################################################
#                                   Metabolomics                                #
#################################################################################

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Mix_metabolome_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Mix_metabolomics_df.RData")

normalized_metabolomics_df <- data.frame(normalized_metabolomics_df)
normalized_metabolomics_df <- cbind(t(metabolomics_metadata[c(1,4),]),normalized_metabolomics_df)

write.table(normalized_metabolomics_df, file="/Users/SLancaster/Desktop/metabolomics_for_supervised_Mix.txt",sep="\t")
            
#Creating the pcl supervised clustering files
#

#################################################################################
#                                      pcl                                      #
#################################################################################

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Arabinoxylan_pcl_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Arabinoxylan_pcl_df.RData")

normalized_pcl_df <- data.frame(normalized_pcl_df)
normalized_pcl_df <- cbind(t(pcl_metadata[c("Participant","Dose"),]),normalized_pcl_df)

write.table(normalized_pcl_df, file="/Users/SLancaster/Desktop/pcl_for_supervised_Arabinoxylan.txt",sep="\t")

#Creating the pcl supervised clustering files
#

#################################################################################
#                                      rna                                      #
#################################################################################

baseline <- c()

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {

load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_rna_df.RData",sep=""))
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_rna_df.RData",sep=""))

normalized_rna_df <- data.frame(normalized_rna_df)
normalized_rna_df <- cbind(t(rna_metadata[c("participant","week"),]),normalized_rna_df)
baseline_rna_df <- normalized_rna_df[normalized_rna_df$week %in% c("Baseline"),]
baseline <- rbind(baseline,baseline_rna_df )
normalized_rna_df <- normalized_rna_df[normalized_rna_df$week %in% c("10","20","30"),]

write.table(normalized_rna_df, file=paste("/Users/SLancaster/Desktop/rna_for_supervised_",fiber_subset,"_fiberonly.txt",sep=""),sep="\t")
}

baseline_drops <- c("X028_LCInulin_Baseline.1","X053_LCInulin_Baseline.1","X1015_lcInulin_Baseline","X1008_Mix_Baseline") #For some reason these rows don't have any reads in them
baseline <- baseline[!rownames(baseline) %in% baseline_drops, ]
write.table(baseline, file="/Users/SLancaster/Desktop/rna_for_supervised_baseline.txt", sep="\t")
baseline$binary <- rep(0,nrow(baseline))

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  rna_df <- read.csv(paste("/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/rna/rna_for_supervised_",fiber_subset,"_fiberonly.txt",sep=""),header=TRUE, sep="\t",row.names = 1)
  
  rna_df$binary <- rep(1,nrow(rna_df))
  rna_df <- rbind(rna_df, baseline)
  
  write.table(rna_df, file=paste("/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/rna/rna_for_supervised_",fiber_subset,"_binary.txt",sep=""),sep="\t")
}

#################################################################################
#                                   lipids                                      #
#################################################################################

baseline <- c()

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  
  load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_lipids_df.RData",sep=""))
  load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_lipids_df.RData",sep=""))
  
  normalized_lipids_df <- data.frame(normalized_lipids_df)
  normalized_lipids_df <- cbind(t(lipids_metadata[c("subject_id","timepoint"),]),normalized_lipids_df)
  baseline_lipids_df <- normalized_lipids_df[normalized_lipids_df$timepoint %in% c("Baseline"),]
  baseline <- rbind(baseline,baseline_lipids_df )
  normalized_lipids_df <- normalized_lipids_df[normalized_lipids_df$timepoint %in% c("10","20","30"),]
  
  write.table(normalized_lipids_df, file=paste("/Users/SLancaster/Desktop/lipids_for_supervised_",fiber_subset,"_fiberonly.txt",sep=""),sep="\t")
}

baseline <- baseline[!rownames(baseline) %in% baseline_drops, ]
write.table(baseline, file="/Users/SLancaster/Desktop/lipids_for_supervised_baseline.txt", sep="\t")
baseline$binary <- rep(0,nrow(baseline))

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  lipids_df <- read.csv(paste("/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/lipids/lipids_for_supervised_",fiber_subset,"_fiberonly.txt",sep=""),header=TRUE, sep="\t",row.names = 1)
  
  lipids_df$binary <- rep(1,nrow(lipids_df))
  lipids_df <- rbind(lipids_df, baseline)
  
  write.table(lipids_df, file=paste("/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/lipids/lipids_for_supervised_",fiber_subset,"_binary.txt",sep=""),sep="\t")
}

#################################################################################
#                                   lipids                                      #
#################################################################################




            