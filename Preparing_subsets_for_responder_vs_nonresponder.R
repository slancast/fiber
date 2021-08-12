#!/usr/bin/R
#
#Dalia is very interested at the starting levels for the interesting genefamilies
#I think the best way to do this would be to make a subest of the genes and then work from there
#

#genes_of_interest <- read.csv("/Users/SLancaster/Desktop/kegg_bile_acids.csv")
genes_of_interest <- read.csv("/Users/SLancaster/Desktop/KeggAXuptotalCluster.csv")

genes <- as.character(genes_of_interest[,1])
genes <- gsub(" ","",genes)

interesting_gene_baselines <- baseline_means[,which(colnames(baseline_means) %in% c("Participant",genes))]

baseline_means <- interesting_gene_baselines

#I was able to then pipe this into the responders vs non-responders pipeline
#I wasn't able to find anything with the fist pass done on 11.18.19

#I am curious about how these also look during the course of fiber treatment
#See below

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Arabinoxylan_genef_df.RData")

interesting_genes <- normalized_genef_df[,which(colnames(normalized_genef_df) %in% c(genes))]

data_frame2 <- interesting_genes

#I will do the same thing here with the bile acids and phenolics sent by Brittany

phenolics <-  read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/Fiber_phenolic_compounds_quantifications.csv", header=TRUE, row.names = NULL)
bile_acids <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/Bile_acids_microbial_metabolite_list.csv", header=TRUE, row.names = NULL)

colnames(normalized_metabolomics_df)

phenolic_names <- make.names(phenolics$Metabolite_val,unique = TRUE)
#Not all phenolics from this dataset are in the larger dataset, so the following snipped of code will select the ones present.
known_phenolcis <- colnames(normalized_metabolomics_df)[which(colnames(normalized_metabolomics_df) %in% phenolic_names)] 

bile_acid_names <- make.names(bile_acids$Metabolite_val,unique = TRUE)
which(colnames(normalized_metabolomics_df) %in% bile_acid_names)

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/metabolomic_baselines.RData")

interesting_metabolomcis_baselines <- baseline_means[,which(colnames(baseline_means) %in% c("Participant",bile_acid_names))]

baseline_means <- interesting_metabolomcis_baselines 

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Arabinoxylan_metabolomics_df.RData")

interesting_metabolites <- normalized_metabolomics_df[,which(colnames(normalized_metabolomics_df) %in% c(bile_acid_names))]

data_frame2 <- interesting_metabolites


