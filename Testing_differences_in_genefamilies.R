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


