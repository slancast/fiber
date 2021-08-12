#############################################################################
#                                 clinicals                                 #
#############################################################################
#To find the baseline for the clinicals the preprocessing of the dataset
#needs to be a little different. This is mainly because of the NA values
#and because the data are much more heterogenous than a typical dataset
#where every value is genearted in the same way.
library(openxlsx)
clinical_df <- read.xlsx("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Data/total_metadata_cleaned.xlsx", colNames = TRUE)
clinical_df <- as.matrix(t( clinical_df))
class(clinical_df) <- "character"
for (i in 1:nrow(clinical_df)){ #going in and taking out some of the unnecessary charaters int the clinical dataset
  j <- clinical_df[i,]
  j <- gsub(" \\(L\\)","", j)
  j <- gsub(" \\(H\\)","", j)
  j <- gsub("<","", j)
  j <- gsub(">","", j)
  clinical_df[i,] <- j
}


clinical_metadata_rows = 5
week <- clinical_df["week",]
participant <- clinical_df["participant",]
clinical_metadata = head(clinical_df,clinical_metadata_rows)

clinical_df2 = tail(clinical_df, -clinical_metadata_rows) #Now getting rid of the metadata
clinical_df2 <- data.frame(t(clinical_df2))
clinical_df2 <- data.frame(lapply(clinical_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(clinical_df2)) #To change to a double, this needs to go through character first
clinical_df2 <- data.frame(lapply(clinical_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(clinical_df2)) #then numeric
clinical_df2 <- clinical_df2[ , -which(colnames(clinical_df2) %in% c("Alb.Creat.Interp","Sex","Ethnicity","IR_IS_from_SSPG","DOB.month.year","date_time_blood_draw_clinicals"))]
rownames(clinical_df2) <- make.names(as.character(unlist(clinical_metadata[1,])), unique=TRUE)

#This in the point where I begin imputing
clinical_toimpute <- clinical_df2[-which(rowMeans(is.na(clinical_df2)) > 0.5),] #eliminating rows missing most of the values
to_remove <- as.vector(which(rowMeans(is.na(clinical_df2)) > 0.5))
week <- week[-to_remove]
participant <- participant[-to_remove]

library(impute) 
clinical_toimpute.i <- impute.knn(t(clinical_toimpute)) #Imputing the individual entries missing 
clinical_imputed <- data.frame(t(clinical_toimpute.i$data)) #pulling out the data from the imputed object
clinical_imputed <- cbind(week,participant,clinical_imputed)

baselines <- clinical_imputed[ which(clinical_imputed$week=='Baseline'),]
baseline_means <- aggregate(baselines,list(Visit=baselines$participant), mean)
baseline_means <- baseline_means[,-2]
baseline_means[,3:ncol(baseline_means)] <- as.numeric(as.matrix(baseline_means[,3:ncol(baseline_means)]))

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/clinical_baselines.RData",sep=""))

#############################################################################
#                                     RNA                                   #
#############################################################################

rna_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/rna-final-combat.txt", sep="\t", header=TRUE, row.names = 1)

rna_metadata_rows = 5
rna_metadata = head(rna_df,rna_metadata_rows)
rna_df <- tail(rna_df, -rna_metadata_rows )
rna_df <- data.frame(t(rna_df))
rna_df <- data.frame(lapply(rna_df, as.character), stringsAsFactors=FALSE, row.names = rownames(rna_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
rna_df <- data.frame(lapply(rna_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(rna_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.

rna_df3 <- as.data.frame(rna_metadata["week",])
rna_df3 <- as.data.frame(t(rna_df3))
partic <- as.data.frame(t(rna_metadata["participant",]))
rna_df3 <- cbind(partic, rna_df3)
rna_df3 <- data.frame(rna_df3)
rna_df3 <- cbind(rna_df3, rna_df)

rna_df4 <- aggregate(rna_df3,list(Participant=rna_df3$participant, Dose=rna_df3$week),  mean) #Changing the week column to Dose
rna_df5 <- rna_df4[,-3]
rna_df5 <- rna_df5[,-3]
baseline_means <- rna_df5[ which(rna_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/rna_baselines.RData",sep=""))

###############################################################################
#                                     Lipids                                   #
###############################################################################

lipids_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Lipids/Data/lipidomics_metadata_transposed_imputed.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header = TRUE)

lipids_metadata_rows = 4
lipids_metadata = head(lipids_df,lipids_metadata_rows)
lipids_df <- tail(lipids_df, -lipids_metadata_rows )
lipids_df <- data.frame(t(lipids_df))
lipids_df <- data.frame(lapply(lipids_df, as.character), stringsAsFactors=FALSE, row.names = rownames(lipids_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
lipids_df <- data.frame(lapply(lipids_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(lipids_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
lipids_df <- log(lipids_df, 2)

lipids_df3 <- as.data.frame(lipids_metadata["timepoint",])
lipids_df3 <- as.data.frame(t(lipids_df3))
partic <- as.data.frame(t(lipids_metadata["subject_id",]))
lipids_df3 <- cbind(partic, lipids_df3)
lipids_df3 <- data.frame(lipids_df3)
lipids_df3 <- cbind(lipids_df3, lipids_df)

lipids_df4 <- aggregate(lipids_df3,list(Participant=lipids_df3$subject_id, Dose=lipids_df3$timepoint),  mean) #Changing the week column to Dose
lipids_df5 <- lipids_df4[,-3]
lipids_df5 <- lipids_df5[,-3]
baseline_means <- lipids_df5[ which(lipids_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/lipids_baselines.RData",sep=""))

#######################################################################################
#                                     pathcoverage                                   #
#######################################################################################

pathcoverage_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/pathcoverage-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)

pathcoverage_metadata_rows = 5
bugs_list_df <- tail(bugs_list_df, -pathcoverage_metadata_rows) #Now getting rid of the metadata
bugs_list_df <- as.matrix(bugs_list_df)
class(bugs_list_df) <- "numeric"
pathcoverage_metadata = data.frame(head(pathcoverage_df,pathcoverage_metadata_rows))

pathcoverage_df2 <- tail(pathcoverage_df, -pathcoverage_metadata_rows) #Now getting rid of the metadata
pathcoverage_df2 <- as.matrix(pathcoverage_df2)
class(pathcoverage_df2) <- "numeric"

pathcoverage_df3 <- as.data.frame(pathcoverage_metadata["Dose",])
pathcoverage_df3 <- as.data.frame(t(pathcoverage_df3))
partic <- as.data.frame(t(pathcoverage_metadata["Participant",]))
pathcoverage_df3 <- cbind(partic, pathcoverage_df3)
pathcoverage_df3 <- cbind(pathcoverage_df3, data.frame(t(pathcoverage_df2)))
pathcoverage_df3 <- cbind(pathcoverage_df3)
pathcoverage_df4 <- aggregate(pathcoverage_df3,list(Participant=pathcoverage_df3$Participant, Dose=pathcoverage_df3$Dose),  mean)
pathcoverage_df5 <- pathcoverage_df4[,-3]
pathcoverage_df5 <- pathcoverage_df5[,-3]
baseline_means <- pathcoverage_df5[ which(pathcoverage_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/pathcoverage_baselines.RData",sep=""))

#######################################################################################
#                                     pathabundance                                   #
#######################################################################################

pathabundance_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/pathabundance-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)

pathabundance_metadata_rows = 5
bugs_list_df <- tail(bugs_list_df, -pathabundance_metadata_rows) #Now getting rid of the metadata
bugs_list_df <- as.matrix(bugs_list_df)
class(bugs_list_df) <- "numeric"
pathabundance_metadata = data.frame(head(pathabundance_df,pathabundance_metadata_rows))

pathabundance_df2 <- tail(pathabundance_df, -pathabundance_metadata_rows) #Now getting rid of the metadata
pathabundance_df2 <- as.matrix(pathabundance_df2)
class(pathabundance_df2) <- "numeric"

pathabundance_df3 <- as.data.frame(pathabundance_metadata["Dose",])
pathabundance_df3 <- as.data.frame(t(pathabundance_df3))
partic <- as.data.frame(t(pathabundance_metadata["Participant",]))
pathabundance_df3 <- cbind(partic, pathabundance_df3)
pathabundance_df3 <- cbind(pathabundance_df3, data.frame(t(pathabundance_df2)))
pathabundance_df3 <- cbind(pathabundance_df3)
pathabundance_df4 <- aggregate(pathabundance_df3,list(Participant=pathabundance_df3$Participant, Dose=pathabundance_df3$Dose),  mean)
pathabundance_df5 <- pathabundance_df4[,-3]
pathabundance_df5 <- pathabundance_df5[,-3]
baseline_means <- pathabundance_df5[ which(pathabundance_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/pathabundance_baselines.RData",sep=""))

###############################################################################
#                                 Metaphlan                                   #
###############################################################################

bugs_list_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)

metaphlan_metadata_rows = 5
metaphlan_metadata = data.frame(head(bugs_list_df,metaphlan_metadata_rows))
bugs_list_df <- tail(bugs_list_df, -metaphlan_metadata_rows) #Now getting rid of the metadata
bugs_list_df <- as.matrix(bugs_list_df)
class(bugs_list_df) <- "numeric"
#bugs_list_df <- bugs_list_df[!rownames(bugs_list_df) %in% "k__Archaea|p__Euryarchaeota|c__Methanobacteria",]
bugs_list_df <- bugs_list_df/100

metaphlan_df3 <- as.data.frame(metaphlan_metadata["Dose",])
metaphlan_df3 <- as.data.frame(t(metaphlan_df3))
partic <- as.data.frame(t(metaphlan_metadata["Participant",]))
metaphlan_df3 <- cbind(partic, metaphlan_df3)
metaphlan_df3 <- cbind(metaphlan_df3, data.frame(t(bugs_list_df)))
metaphlan_df3 <- cbind(metaphlan_df3)
metaphlan_df4 <- aggregate(metaphlan_df3,list(Participant=metaphlan_df3$Participant, Dose=metaphlan_df3$Dose),  mean)
metaphlan_df5 <- metaphlan_df4[,-3]
metaphlan_df5 <- metaphlan_df5[,-3]
baseline_means <- metaphlan_df5[ which(metaphlan_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/metaphlan_baselines.RData",sep=""))

##################################################################################
#                                 Genefamilies                                   #
##################################################################################

genef_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/genefamilies-90ko-groupedonly-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)

genef_metadata_rows = 5
genef_metadata = head(genef_df,genef_metadata_rows)
genef_df <- tail(genef_df, -genef_metadata_rows )
genef_df <- data.frame(t(genef_df))
genef_df <- data.frame(lapply(genef_df, as.character), stringsAsFactors=FALSE, row.names = rownames(genef_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
genef_df <- data.frame(lapply(genef_df, as.numeric), stringsAsFactors=FALSE, row.names = rownames(genef_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.

genef_df3 <- as.data.frame(genef_metadata["Dose",])
genef_df3 <- as.data.frame(t(genef_df3))
partic <- as.data.frame(t(genef_metadata["Participant",]))
genef_df3 <- cbind(partic, genef_df3)
genef_df3 <- data.frame(genef_df3)
genef_df3 <- cbind(genef_df3, genef_df)

genef_df4 <- aggregate(genef_df3,list(Participant=genef_df3$Participant, Dose=genef_df3$Dose),  mean) #Changint the week column to Dose
genef_df5 <- genef_df4[,-3]
genef_df5 <- genef_df5[,-3]
baseline_means <- genef_df5[ which(genef_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/genef_baselines.RData",sep=""))

#############################################################################
#                              Cytokines                                    #
#############################################################################

datamatrix = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Cytokine/Data/CytokineWashoutFinal.txt", sep="\t", header=TRUE)
#Getting rid of the CHEX
datamatrix <- datamatrix[,-ncol(datamatrix)]
datamatrix <- datamatrix[,-ncol(datamatrix)]
datamatrix <- datamatrix[,-ncol(datamatrix)]
datamatrix <- datamatrix[,-ncol(datamatrix)]
#continuing with standard analysis
metadata_rows = 5
datamatrix <- data.frame(t(datamatrix))
cytokine_metadata <- head(datamatrix,metadata_rows)
metadata <- cytokine_metadata
datamatrix2 <- tail(datamatrix, -metadata_rows) #Now getting rid of the metadata
datamatrix2 <- as.matrix(t(datamatrix2))
class(datamatrix2) <- "numeric"
cytokine_df <- log(datamatrix2, 2)

cytokine_df3 <- as.data.frame(cytokine_metadata["Week",])
cytokine_df3 <- as.data.frame(t(cytokine_df3))
partic <- as.data.frame(t(cytokine_metadata["Participant",]))
cytokine_df3 <- cbind(partic, cytokine_df3)
cytokine_df3 <- cbind(cytokine_df3, cytokine_df)

cytokine_df4 <- aggregate(cytokine_df3,list(Participant=cytokine_df3$Participant, Dose=cytokine_df3$Week),  mean)
cytokine_df5 <- cytokine_df4[,-3]
cytokine_df5 <- cytokine_df5[,-3]
cytokine_df5 <- cytokine_df5[,-3]
baseline_means <- cytokine_df5[which(cytokine_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/cytokine_baselines.RData",sep=""))

#############################################################################
#                              Metabolome                                   #
#############################################################################

metabolomics_metadata_rows = 4

metabolomics_df <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/fiber_metabolites_named.txt", header = TRUE, sep="\t")
rownames(metabolomics_df) <- make.names(metabolomics_df[,1], unique=TRUE)
metabolomics_df <- metabolomics_df[,-1]

metabolomics_df <- data.frame(t(metabolomics_df))
meta_metadata <- head(metabolomics_df, metabolomics_metadata_rows)
metabolomics_df2 = tail(metabolomics_df, -metabolomics_metadata_rows) #Now getting rid of the metadata
metabolomics_df2 <- data.frame(t(metabolomics_df2))
metabolomics_df2 <- data.frame(lapply(metabolomics_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df2)) #To change to a double, this needs to go through character first
metabolomics_df2 <- data.frame(lapply(metabolomics_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df2)) #then numeric
metabolomics_df2 <- log(metabolomics_df2, 2)

meta_df3 <- as.data.frame(meta_metadata["timepoint",])
meta_df3 <- as.data.frame(t(meta_df3))
partic <- as.data.frame(t(meta_metadata["subject_id",]))
meta_df3 <- cbind(partic, meta_df3)
meta_df3 <- cbind(meta_df3, metabolomics_df2)

meta_df4 <- aggregate(meta_df3,list(Participant=meta_df3$subject_id, Dose=meta_df3$timepoint),  mean)
meta_df5 <- meta_df4[,-3]
meta_df5 <- meta_df5[,-3]
baseline_means <- meta_df5[ which(meta_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/metabolomic_baselines.RData",sep=""))

#############################################################################
#                              Proteomics                                   #
#############################################################################

load("/Users/SLancaster/Desktop/Projects/Fiber/Proteomics/Data/proteomics_with_metadata.RData")
proteomics_metadata_rows = 64

prot_metadata <- head(proteomics_with_metadata, proteomics_metadata_rows)
proteomics_df2 <- tail(proteomics_with_metadata, -proteomics_metadata_rows) #Now getting rid of the metadata
proteomics_df2 <- data.frame(t(proteomics_df2))
proteomics_df2 <- data.frame(lapply(proteomics_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(proteomics_df2)) #To change to a double, this needs to go through character first
proteomics_df2 <- data.frame(lapply(proteomics_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(proteomics_df2)) #then numeric

prot_df3 <- as.data.frame(prot_metadata["week",])
prot_df3 <- as.data.frame(t(prot_df3))
partic <- as.data.frame(t(prot_metadata["participant",]))
prot_df3 <- cbind(partic, prot_df3)
prot_df3 <- cbind(prot_df3, proteomics_df2)
prot_df4 <- aggregate(prot_df3,list(Participant=prot_df3$participant, Dose=prot_df3$week),  mean)
prot_df5 <- prot_df4[,-3]
prot_df5 <- prot_df5[,-3]
baseline_means <- prot_df5[ which(prot_df5$Dose=='Baseline'),]

save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/proteomic_baselines.RData",sep=""))

