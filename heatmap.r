#clutering the omics data
#First I attach the timepoint to the data

fiber_subset = "Arabinoxylan"
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_metabolome_df.RData",sep=""))

#rename the rows (samples), so they are more informative
metabolomics_clustring <- normalized_metabolomics_df
metabolomics_clustring <- cbind(paste(t(metabolomics_metadata["timepoint",]),t(metabolomics_metadata["subject_id",])),metabolomics_clustring)
rownames(metabolomics_clustring) <- metabolomics_clustring[,1]
metabolomics_clustring <- metabolomics_clustring[,-1]

# The make it numeric
class(metabolomics_clustring) <- "numeric"

#plotting
pdf(paste("/Users/SLancaster/Desktop/metabolomicsheatmap",fiber_subset,".pdf",sep=""), width = 15, height = 15 )
par(xaxt="n",mar=c(5,6,4,1), oma=c(0,0,0,3))
heatmap(metabolomics_clustring, xaxt='n', ann=FALSE, xlab="")
dev.off()

#Create new dataframe to write for the supervised clustering
normalized_metabolomics_df <- cbind(t(metabolomics_metadata["timepoint",]), normalized_metabolomics_df)
write.table(normalized_metabolomics_df, file=paste("/Users/SLancaster/Desktop/metabolomics_",fiber_subset,".txt",sep=""), sep="\t")

##################################################
##################### PCL ########################
##################################################

fiber_subset = "Mix"
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""))
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_pcl_df.RData",sep=""))

#rename the rows (samples), so they are more informative
pcl_clustring <- normalized_pcl_df
pcl_clustring <- cbind(paste(t(pcl_metadata["Dose",]),t(pcl_metadata["Participant",])),pcl_clustring)
rownames(pcl_clustring) <- pcl_clustring[,1]
pcl_clustring <- pcl_clustring[,-1]

# The make it numeric
class(pcl_clustring) <- "numeric"

#plotting
pdf(paste("/Users/SLancaster/Desktop/pclheatmap",fiber_subset,".pdf",sep=""), width = 15, height = 15 )
par(xaxt="n",mar=c(5,6,4,1), oma=c(0,0,0,3))
heatmap(pcl_clustring, xaxt='n', ann=FALSE, xlab="")
dev.off()

#Create new dataframe to write for the supervised clustering
normalized_pcl_df <- cbind(t(pcl_metadata["Dose",]), normalized_pcl_df)
write.table(normalized_pcl_df, file=paste("/Users/SLancaster/Desktop/pcl_",fiber_subset,".txt",sep=""), sep="\t")

##################################################
################## Cytokine ######################
##################################################

fiber_subset = "Arabinoxylan"
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_cytokine_df.RData",sep=""))
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_cytokine_df.RData",sep=""))

#rename the rows (samples), so they are more informative
cytokine_clustring <- normalized_cytokine_df
cytokine_clustring <- cbind(paste(t(cytokine_metadata["Week",]),t(cytokine_metadata["Participant",])),cytokine_clustring)
rownames(cytokine_clustring) <- cytokine_clustring[,1]
cytokine_clustring <- cytokine_clustring[,-1]

# The make it numeric
class(cytokine_clustring) <- "numeric"

#plotting
pdf(paste("/Users/SLancaster/Desktop/cytokineheatmap",fiber_subset,".pdf",sep=""), width = 15, height = 15 )
par(mar=c(5,6,4,1), oma=c(0,0,0,3))
heatmap(cytokine_clustring)
dev.off()

##################################################
##################### RNA ########################
##################################################
fiber_subset = "Mix"
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_rna_df.RData",sep=""))
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_rna_df.RData",sep=""))

normalized_rna_df <- cbind(t(rna_metadata["week",]), normalized_rna_df)
write.table(normalized_rna_df, file=paste("/Users/SLancaster/Desktop/iPOP_machine_learning/rna_",fiber_subset,".txt",sep=""), sep="\t")

