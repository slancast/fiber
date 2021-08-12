for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
print(fiber_subset) #

if (fiber_subset == "Mix") {
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}

source("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Software/Utils.r")
#A lot of these could probably be turned into funcitons. For now it seems organized enough to use them like this.
########################
########Clinical########
########################
if (FALSE) {
  
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
clinical_df <- data.frame(t(clinical_df)) #Turning the matrix so that fiber is a column making subsetting easier
clinical_df <- data.frame(lapply(clinical_df, as.character), stringsAsFactors=FALSE, row.names = rownames(clinical_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
clinical_df <- clinical_df[which(clinical_df$fiber==fiber_subset),] #subsetting by fiber
clinical_df <- data.frame(t(clinical_df))
clinical_metadata = head(clinical_df,clinical_metadata_rows) #snipping of the metadata, that is the top rows

clinical_df2 = tail(clinical_df, -clinical_metadata_rows) #Now getting rid of the metadata
clinical_df2 <- data.frame(t(clinical_df2))
clinical_df2 <- data.frame(lapply(clinical_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(clinical_df2)) #To change to a double, this needs to go through character first
clinical_df2 <- data.frame(lapply(clinical_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(clinical_df2)) #then numeric
clinical_df2 <- clinical_df2[ , -which(colnames(clinical_df2) %in% c("Alb.Creat.Interp","Sex","Ethnicity","IR_IS_from_SSPG","DOB.month.year","date_time_blood_draw_clinicals"))]
rownames(clinical_df2) <- make.names(as.character(unlist(clinical_metadata[1,])), unique=TRUE)

#This in the point where I begin imputing
clinical_toimpute <- clinical_df2[-which(rowMeans(is.na(clinical_df2)) > 0.5),] #eliminating rows missing most of the values
library(impute) 
clinical_toimpute.i <- impute.knn(t(clinical_toimpute)) #Imputing the individual entries missing 
clinical_imputed <- data.frame(t(clinical_toimpute.i$data)) #pulling out the data from the imputed object

save(clinical_metadata, clinical_imputed, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_clinical_df.RData",sep=""))

clinical_df3 <- as.data.frame(clinical_metadata["week",]) #pulling out the metadata about the week
clinical_df3 <- as.data.frame(t(clinical_df3))
#The clinical dataset is missing so many values and this leads
#to many NAs after aggregating by the mean
clinical_df3 <- clinical_df3[-which(rowMeans(is.na(clinical_df2)) > 0.5),] #Throwing out the same rows with mostly missing data from the metadata
clinical_df3 <- cbind(clinical_df3, clinical_imputed )
colnames(clinical_df3)[1] <- "week" #Using this order of commands the name of the column containing the week information has changed
clinical_df4 <- aggregate(clinical_df3,list(Visit=clinical_df3$week), mean) #Starting the classic order for aggregation.
clinical_df4 <- clinical_df4[ , -2]
clinical_df4 <- clinical_df4[match(order, clinical_df4$Visit),]
clinical_df4 <- data.frame(t(na.omit(t(clinical_df4))))
rownames(clinical_df4) <- clinical_df4$Visit
aggregated_clinical_df <- clinical_df4[,-1]

save(clinical_metadata, aggregated_clinical_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_clinical_df.RData",sep=""))

clinical_metadata <- clinical_metadata[,-which(rowMeans(is.na(clinical_df2)) > 0.5)] #Throwing out the same rows with mostly missing data from the metadata
clinical_df2 <- clinical_imputed
save(clinical_metadata, clinical_df2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_clinical_df.RData",sep=""))

log_clinicals <- clinical_df2 + 1 #I add one to every data frame to avoid the probem of 0s.
log_clinicals <- log(clinical_df2, 2)

save(clinical_metadata, log_clinicals, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_clinical_df.RData",sep=""))

}
######################
########Lipids########
######################

if (FALSE) {

lipids_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Lipids/Data/lipidomics_metadata_transposed.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header = TRUE)

#Full#

lipids_metadata_rows = 4 #Assigning the number of rows first
#Turning the matrix so that fiber is a column making subsetting easier
lipids_metadata <- head(lipids_df,lipids_metadata_rows) 
lipids_df2 <- tail(lipids_df,-lipids_metadata_rows)

lipids_toimpute <- lipids_df2[-which(rowMeans(is.na(lipids_df2)) > 0.5),] #eliminating rows missing most of the values
lipids_toimpute <- as.matrix(lipids_toimpute)
class(lipids_toimpute) <- "double"
library(impute) 
lipids_toimpute.i <- impute.knn(lipids_toimpute) #Imputing the individual entries missing 
lipids_imputed <- data.frame(t(lipids_toimpute.i$data)) #pulling out the data from the imputed object
lipids_df2 <- data.frame(t(lipids_imputed))
lipids_df2 <- rbind(lipids_metadata, lipids_df2 )

write.table(lipids_df2 , "/Users/SLancaster/Desktop/Projects/Fiber/Lipids/Data/lipidomics_metadata_transposed_imputed.csv", sep=",")

lipids_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Lipids/Data/lipidomics_metadata_transposed_imputed.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header = TRUE)

lipids_df <- data.frame(t(lipids_df))
lipids_df <- data.frame(lapply(lipids_df, as.character), stringsAsFactors=FALSE, row.names = rownames(lipids_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
lipids_df <- lipids_df[which(lipids_df$fiber==fiber_subset),] #subsetting by fiber
lipids_df <- data.frame(t(lipids_df))

lipids_metadata <- head(lipids_df,lipids_metadata_rows) 
lipids_df2 <- tail(lipids_df,-lipids_metadata_rows)

#Eliminating the sample name and fiber columns
#lipids_df <- lipids_df[,-which(names(lipids_df) %in% c("sample_name","fiber"))] #At this point with the new function averaging the replicates this command should eliminate all but the two columns -- participant and time point columns

#lipids_df <- averaging_replicates(lipids_df)

lipids_df2 <- data.frame(lapply(lipids_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(lipids_df2)) #To change to a double, this needs to go through character first
lipids_df2 <- data.frame(lapply(lipids_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(lipids_df2)) #then numeric

save(lipids_metadata, lipids_df2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_lipids_df.RData",sep=""))

#Just like the metabolomics, the lipidomics had additional files that have additional details
#about the about the compounds found in the mass spectra. Again like the metabolomics I will
#attach these files to the log data file for future accessing.

load("/Users/SLancaster/Desktop/Projects/Fiber/Lipids/Data/Lipids.RData")

lipid_class <- gsub("\\..*","",colnames(lipids_df2))
lipid_class <- gsub("[0-9].*","",lipid_class)
names(lipid_class) <- colnames(lipids_df2)

log_lipids <- log(lipids_df2+1, 2)
save(lipids_metadata, log_lipids, lipid_class, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_lipids_df.RData",sep=""))

#Aggregate#

lipids_df3 <- as.data.frame(lipids_metadata["timepoint",])
lipids_df3 <- as.data.frame(t(lipids_df3))
lipids_df3 <- cbind(lipids_df3, t(lipids_df2))
lipids_df3 <- aggregate(lipids_df3,list(Visit=lipids_df3$timepoint), mean)
lipids_df3 <- lipids_df3[ , -2]
lipids_df3 <- lipids_df3[match(order, lipids_df3$Visit),]
lipids_df3 <- data.frame(t(na.omit(t(lipids_df3))))
rownames(lipids_df3) <- lipids_df3$Visit
lipids_df3 <- lipids_df3[,-1]

save(lipids_metadata, lipids_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_lipids_df.RData",sep=""))

}


###################
########RNA########
################### 

if (FALSE) {

rna_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/rna-final-combat.txt", sep="\t", row.names = 1, stringsAsFactors=FALSE)

#Full#

rna_metadata_rows = 5
rna_df <- data.frame(t(rna_df))
rna_df <- data.frame(lapply(rna_df, as.character), stringsAsFactors=FALSE, row.names = rownames(rna_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
rna_df <- rna_df[which(rna_df$condition==fiber_subset),] #subsetting by fiber
rna_df <- data.frame(t(rna_df))
rna_metadata = head(rna_df,rna_metadata_rows)
rna_df2 = tail(rna_df, -rna_metadata_rows) #Now getting rid of the metadata

for_averaging <- data.frame(rbind(rna_metadata["participant",],rna_metadata["week",]))
for_averaging <- data.frame(cbind(t(for_averaging),t(rna_df2 )))
rownames(for_averaging) <- colnames(rna_df2)
rna_df2 <- averaging_replicates(for_averaging)
rna_metadata <- data.frame(t(rna_df2[,1:2]))
rna_df2 <- rna_df2[,-1]
rna_df2 <- rna_df2[,-1]

log_rna <- rna_df2
save(rna_metadata, log_rna, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_rna_df.RData",sep=""))

#Aggregate#

rna_df3 <- as.data.frame(rna_metadata["Visit",])
rna_df3 <- as.data.frame(t(rna_df3))
rna_df3 <- cbind(rna_df3, rna_df2)
rna_df3 <- aggregate(rna_df3,list(Visit=rna_df3$Visit), mean)
rna_df3 <- rna_df3[ , -2]
rna_df3 <- rna_df3[match(order, rna_df3$Visit),]
rna_df3 <- na.omit(rna_df3)
rownames(rna_df3) <- rna_df3$Visit
rna_df3 <- rna_df3[,-1]

save(rna_metadata, rna_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_RNA_df.RData",sep=""))
}


###################################
###########pathcoverage############
###################################
if (FALSE) {
  pathcoverage_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/pathcoverage-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)
  
  #Full#
  
  pathcoverage_metadata_rows = 5
  pathcoverage_df <- data.frame(t(pathcoverage_df))
  pathcoverage_df <- pathcoverage_df[which(pathcoverage_df$Fiber==fiber_subset),] #subsetting by fiber
  pathcoverage_df <- data.frame(t(pathcoverage_df))
  pathcoverage_metadata = head(pathcoverage_df,pathcoverage_metadata_rows)
  pathcoverage_df2 <- tail(pathcoverage_df, -pathcoverage_metadata_rows) #Now getting rid of the metadata
  
  for_averaging <- data.frame(rbind(pathcoverage_metadata["Participant",],pathcoverage_metadata["Dose",]))
  for_averaging <- data.frame(cbind(t(for_averaging),t(pathcoverage_df2 )))
  rownames(for_averaging) <- colnames(pathcoverage_df2)
  pathcoverage_df2 <- averaging_replicates(for_averaging)
  pathcoverage_metadata <- data.frame(t(pathcoverage_df2[,1:2]))
  pathcoverage_df2 <- pathcoverage_df2[,-1]
  pathcoverage_df2 <- pathcoverage_df2[,-1]
  
  log_pathcoverage <- pathcoverage_df2

  save(pathcoverage_metadata, log_pathcoverage, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_pathcoverage_df.RData",sep=""))
  
  #Aggregate#
  
  pathcoverage_df2 <- as.matrix(pathcoverage_df2)
  class(pathcoverage_df2) <- "numeric"
  pathcoverage_df3 <- as.data.frame(pathcoverage_metadata["Visit",])
  pathcoverage_df3 <- as.data.frame(t(pathcoverage_df3))
  pathcoverage_df3 <- cbind(pathcoverage_df3, pathcoverage_df2)
  pathcoverage_df3 <- aggregate(pathcoverage_df3,list(Visit=pathcoverage_df3$Visit), mean)
  pathcoverage_df3 <- pathcoverage_df3[ , -2]
  pathcoverage_df3 <- pathcoverage_df3[match(order, pathcoverage_df3$Visit),]
  pathcoverage_df3 <- na.omit(pathcoverage_df3)
  rownames(pathcoverage_df3) <- pathcoverage_df3$Visit
  pathcoverage_df3 <- pathcoverage_df3[,-1]
  
  save(pathcoverage_metadata, pathcoverage_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_pathcoverage_df.RData",sep=""))
}



####################################
###########pathabundance############
####################################
if (FALSE) {
pathabundance_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/pathabundance-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)

#Full#

pathabundance_metadata_rows = 5
pathabundance_df <- data.frame(t(pathabundance_df))
pathabundance_df <- pathabundance_df[which(pathabundance_df$Fiber==fiber_subset),] #subsetting by fiber
pathabundance_df <- data.frame(t(pathabundance_df))
pathabundance_metadata = head(pathabundance_df,pathabundance_metadata_rows)
pathabundance_df2 <- tail(pathabundance_df, -pathabundance_metadata_rows) #Now getting rid of the metadata

for_averaging <- data.frame(rbind(pathabundance_metadata["Participant",],pathabundance_metadata["Dose",]))
for_averaging <- data.frame(cbind(t(for_averaging),t(pathabundance_df2)))
rownames(for_averaging) <- colnames(pathabundance_df2)
pathabundance_df2 <- averaging_replicates(for_averaging)
pathabundance_metadata <- data.frame(t(pathabundance_df2[,1:2]))
pathabundance_df2 <- pathabundance_df2[,-1]
pathabundance_df2 <- pathabundance_df2[,-1]

log_pathabundance <- pathabundance_df2

save(pathabundance_metadata, log_pathabundance, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_pathabundance_df.RData",sep=""))

#Aggregate#

pathabundance_df2 <- as.matrix(pathabundance_df2)
class(pathabundance_df2) <- "numeric"
pathabundance_df3 <- as.data.frame(pathabundance_metadata["Visit",])
pathabundance_df3 <- as.data.frame(t(pathabundance_df3))
pathabundance_df3 <- cbind(pathabundance_df3, pathabundance_df2)
pathabundance_df3 <- aggregate(pathabundance_df3,list(Visit=pathabundance_df3$Visit), mean)
pathabundance_df3 <- pathabundance_df3[ , -2]
pathabundance_df3 <- pathabundance_df3[match(order, pathabundance_df3$Visit),]
pathabundance_df3 <- na.omit(pathabundance_df3)
rownames(pathabundance_df3) <- pathabundance_df3$Visit
pathabundance_df3 <- pathabundance_df3[,-1]

save(pathabundance_metadata, pathabundance_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_pathabundance_df.RData",sep=""))
}

################################
###########Metaphlan############
################################

#if (FALSE) {

bugs_list_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)
metaphlan_metadata_rows = 5

bugs_list_df <- data.frame(t(bugs_list_df))
bugs_list_df<- bugs_list_df[which(bugs_list_df$Fiber==fiber_subset),] #subsetting by fiber
bugs_list_df <- data.frame(t(bugs_list_df))
metaphlan_metadata = head(bugs_list_df,metaphlan_metadata_rows)
bugs_list_df2 <- tail(bugs_list_df, -metaphlan_metadata_rows) #Now getting rid of the metadata
bugs_list_df2 <- as.matrix(bugs_list_df2)
#bugs_list_df2 <- bugs_list_df2[!rownames(bugs_list_df2) %in% "k__Archaea.p__Euryarchaeota.c__Methanobacteria", ]
class(bugs_list_df2) <- "numeric"
#bugs_list_df2 <- data.frame(bugs_list_df2/100)

for_averaging <- data.frame(rbind(metaphlan_metadata["Participant",],metaphlan_metadata["Dose",]))
for_averaging <- data.frame(cbind(t(for_averaging),t(bugs_list_df2 )))
rownames(for_averaging) <- colnames(bugs_list_df2)
metaphlan_df2 <- averaging_replicates(for_averaging)
metaphlan_metadata <- data.frame(t(metaphlan_df2[,1:2]))
metaphlan_df2 <- metaphlan_df2[,-1]
metaphlan_df2 <- metaphlan_df2[,-1]

log_metaphlan <- metaphlan_df2

save(metaphlan_metadata, log_metaphlan, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))

#Aggregate#
metaphlan_df2 <- as.matrix(metaphlan_df2)
class(metaphlan_df2) <- "numeric"
metaphlan_df3 <- as.data.frame(metaphlan_metadata["Visit",])
metaphlan_df3 <- as.data.frame(t(metaphlan_df3))
metaphlan_df3 <- cbind(metaphlan_df3, metaphlan_df2)
metaphlan_df3 <- aggregate(metaphlan_df3,list(Visit=metaphlan_df3$Visit), mean)
metaphlan_df3 <- metaphlan_df3[ , -2]
metaphlan_df3 <- metaphlan_df3[match(order, metaphlan_df3$Visit),]
metaphlan_df3 <- na.omit(metaphlan_df3)
rownames(metaphlan_df3) <- metaphlan_df3$Visit
metaphlan_df3 <- metaphlan_df3[,-1]

save(metaphlan_metadata, metaphlan_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_metaphlan_df.RData",sep=""))

#}

###################################
###########Genefamilies############
###################################

if (FALSE) {
genef_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/genefamilies-90ko-groupedonly-final-combat.pcl", sep="\t", row.names = 1, stringsAsFactors=FALSE)

#Full#

genef_metadata_rows = 5
genef_df <- data.frame(t(genef_df))
genef_df <- data.frame(lapply(genef_df, as.character), stringsAsFactors=FALSE, row.names = rownames(genef_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
genef_df <- genef_df[which(genef_df$Fiber==fiber_subset),] #subsetting by fiber
genef_df <- data.frame(t(genef_df))
genef_metadata = head(genef_df,genef_metadata_rows)
genef_df2 = tail(genef_df, -genef_metadata_rows) #Now getting rid of the metadata

for_averaging <- data.frame(rbind(genef_metadata["Participant",],genef_metadata["Dose",]))
for_averaging <- data.frame(cbind(t(for_averaging),t(genef_df2 )))
rownames(for_averaging) <- colnames(genef_df2)
genef_df2 <- averaging_replicates(for_averaging)
genef_metadata <- data.frame(t(genef_df2[,1:2]))
genef_df2 <- genef_df2[,-1]
genef_df2 <- genef_df2[,-1]

log_genef <- genef_df2

save(genef_metadata, log_genef, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_genef_df.RData",sep=""))

genef_df2 <- as.matrix(genef_df2)
class(genef_df2) <- "numeric"
genef_df3 <- as.data.frame(genef_metadata["Visit",])
genef_df3 <- as.data.frame(t(genef_df3))
genef_df3 <- cbind(genef_df3, genef_df2)
genef_df3 <- aggregate(genef_df3,list(Visit=genef_df3$Visit), mean)
genef_df3 <- genef_df3[ , -2]
genef_df3 <- genef_df3[match(order, genef_df3$Visit),]
genef_df3 <- na.omit(genef_df3)
rownames(genef_df3) <- genef_df3$Visit
genef_df3 <- genef_df3[,-1]

save(genef_metadata, genef_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_genef_df.RData",sep=""))
}

##########################
########Cytokine##########
##########################
if (FALSE) {

cytokine_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Cytokine/Data/CytokineWashoutFinal.txt", sep="\t", header=TRUE)
cytokine_df = data.frame(t(cytokine_df))

#Full#

metadata_rows = 6

cytokine_df <- data.frame(t(cytokine_df))
#Getting rid of the CHEX
cytokine_df <- cytokine_df[,-ncol(cytokine_df)]
cytokine_df <- cytokine_df[,-ncol(cytokine_df)]
cytokine_df <- cytokine_df[,-ncol(cytokine_df)]
cytokine_df <- cytokine_df[,-ncol(cytokine_df)]
#continuing with standard analysis
cytokine_df <- data.frame(lapply(cytokine_df, as.character), stringsAsFactors=FALSE, row.names = rownames(cytokine_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
cytokine_df <- cytokine_df[which(cytokine_df$Fiber==fiber_subset),] #subsetting by fiber
cytokine_df <- data.frame(t(cytokine_df))
cytokine_metadata <- head(cytokine_df, metadata_rows) #Putting metadata into its own matrix


#Since the cytokine data doesn't have any sample names, I'll just create them here
participant <- as.character(unlist(cytokine_metadata["Participant",])) #to paste the lists must be character vectors
fiber <- as.character(unlist(cytokine_metadata["Fiber",]))
week <- as.character(unlist(cytokine_metadata["Week",]))
sample_names <- paste(participant, fiber, week, sep="_")
colnames(cytokine_metadata) <- sample_names
colnames(cytokine_df) <- sample_names
#ending creation of sample names

#The second df has the metadata removed
cytokine_df2 = data.frame(tail(cytokine_df, -metadata_rows)) #Now getting rid of the metadata
cytokine_df2 <- data.frame(lapply(cytokine_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(cytokine_df2)) #This needs to be done for the transcript data, but causes problems with the pcl data.
cytokine_df2 <- data.frame(lapply(cytokine_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(cytokine_df2)) #This needs to be done for the transcript data, but causes problems with the pcl data.
cytokine_df2 <- data.frame(t(cytokine_df2))

save(cytokine_metadata, cytokine_df2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_cytokine_df.RData",sep=""))

log_cytokine <- cytokine_df2
log_cytokine <- log_cytokine + 1
log_cytokine <- log(log_cytokine, 2)

save(cytokine_metadata, log_cytokine, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_cytokine_df.RData",sep=""))

#Aggregate#

datamatrix3 <- as.data.frame(cytokine_metadata["Week",])
datamatrix3 <- as.data.frame(t(datamatrix3))
datamatrix3 <- cbind(datamatrix3, cytokine_df2)
datamatrix3 <- aggregate(datamatrix3,list(Visit=datamatrix3$Week), mean)
datamatrix3 <- datamatrix3[ , -2]
datamatrix3 <- datamatrix3[match(order, datamatrix3$Visit),]
datamatrix3 <- na.omit(datamatrix3)
rownames(datamatrix3) <- datamatrix3$Visit
datamatrix3 <- datamatrix3[,-1]

save(cytokine_metadata, cytokine_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_cytokine_df.RData",sep=""))
}

###########################
########Metabolome#########
###########################
if (FALSE) {

if (fiber_subset == "Mix") {
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","Washout_Final")
}

#Full#

metabolomics_metadata_rows = 4
metabolomics_df <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/fiber_metabolites_named.txt", header = TRUE,sep="\t")
rownames(metabolomics_df) <- make.names(metabolomics_df[,1], unique=TRUE)
metabolomics_df <- metabolomics_df[,-1]

metabolomics_df <- data.frame(lapply(metabolomics_df, as.character), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
metabolomics_df <- metabolomics_df[which(metabolomics_df$fiber==fiber_subset),] #subsetting by fiber
metabolomics_df <- data.frame(t(metabolomics_df))
metabolomics_metadata = head(metabolomics_df,metabolomics_metadata_rows)

metabolomics_df2 = tail(metabolomics_df, -metabolomics_metadata_rows) #Now getting rid of the metadata
metabolomics_df2 <- data.frame(t(metabolomics_df2))
metabolomics_df2 <- data.frame(lapply(metabolomics_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df2)) #To change to a double, this needs to go through character first
metabolomics_df2 <- data.frame(lapply(metabolomics_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(metabolomics_df2)) #then numeric

save(metabolomics_metadata, metabolomics_df2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_metabolomics_df.RData",sep=""))

log_metabolomics <- metabolomics_df2
log_metabolomics <- log_metabolomics + 1
log_metabolomics <- log(log_metabolomics, 2)#negative values generate NaNs


annotated_metabolomics_data <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Metabolomics/Data/data_validated_transposed.csv",row.names=1,header=TRUE)
metabolome_superpathway <- annotated_metabolomics_data["Super.pathway",]
metabolome_superpathway[1] <- NA
metabolome_superpathway <- as.character(as.matrix(metabolome_superpathway))
metabolome_superpathway[is.na(metabolome_superpathway)] <- "Other"
metabolome_superpathway <- data.frame(t(metabolome_superpathway))
colnames(metabolome_superpathway) <- colnames(log_metabolomics)
log_metabolomics <- rbind(log_metabolomics, metabolome_superpathway)
log_metabolomics <- data.frame(t(na.omit(t(log_metabolomics)))) #Tossing out the metabolites that contain NaNs
metabolome_superpathway <- log_metabolomics["1",]
log_metabolomics <- log_metabolomics[-nrow(log_metabolomics),]

save(metabolomics_metadata, log_metabolomics, metabolome_superpathway, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))

#Aggregate#

metabolomics_df3 <- as.data.frame(metabolomics_metadata["timepoint",])
metabolomics_df3 <- as.data.frame(t(metabolomics_df3))
metabolomics_df3 <- cbind(metabolomics_df3, metabolomics_df2)
metabolomics_df3 <- aggregate(metabolomics_df3,list(Visit=metabolomics_df3$timepoint), mean)
metabolomics_df3 <- metabolomics_df3[ , -2]
metabolomics_df3 <- metabolomics_df3[match(order, metabolomics_df3$Visit),]
metabolomics_df3 <- na.omit(metabolomics_df3)
rownames(metabolomics_df3) <- metabolomics_df3$Visit
metabolomics_df3<- metabolomics_df3[,-1]

save(metabolomics_metadata, metabolomics_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_metabolomics_df.RData",sep=""))
}


###########################
########Proteome###########
###########################
if (FALSE) {

if (fiber_subset == "Mix") {
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","Washout_Final")
}

#Full#

load("/Users/SLancaster/Desktop/Projects/Fiber/Proteomics/Data/proteomics_with_metadata.RData")
proteomics_metadata_rows = 64

proteomics_df <- data.frame(t(proteomics_with_metadata))
proteomics_df <- data.frame(lapply(proteomics_df, as.character), stringsAsFactors=FALSE, row.names = rownames(proteomics_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
proteomics_df <- proteomics_df[which(proteomics_df$fiber==fiber_subset),] #subsetting by fiber
proteomics_df <- data.frame(t(proteomics_df))
proteomics_metadata <- head(proteomics_df,proteomics_metadata_rows)

proteomics_df2 <- tail(proteomics_df, -proteomics_metadata_rows) #Now getting rid of the metadata
proteomics_df2 <- data.frame(t(proteomics_df2))
proteomics_df2 <- data.frame(lapply(proteomics_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(proteomics_df2)) #To change to a double, this needs to go through character first
proteomics_df2 <- data.frame(lapply(proteomics_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(proteomics_df2)) #then numeric
proteomics_df2 <- 2^proteomics_df2 #undoing the log

save(proteomics_metadata, proteomics_df2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_proteomics_df.RData",sep=""))

#The proteomics data are already log transformed.

log_proteomics <- proteomics_df2
log_proteomics <- log_proteomics + 1 #Standard if there are any 0s
log_proteomics <- log(log_proteomics, 2)#negative values generate NaNs

load("/Users/SLancaster/Desktop/Projects/Fiber/Proteomics/Data/Proteomics.RData")
protein_annotations <- cbind(as.character(Proteomics_data_1pept_Log2_filt80perc_ImptND_Combat$Protein.name), as.character(Proteomics_data_1pept_Log2_filt80perc_ImptND_Combat$Gene.name))
rownames(protein_annotations) <- rownames(Proteomics_data_1pept_Log2_filt80perc_ImptND_Combat)

save(proteomics_metadata, log_proteomics, protein_annotations, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_proteomics_df.RData",sep=""))

#Aggregate#

proteomics_df3 <- as.data.frame(proteomics_metadata["week",])
proteomics_df3 <- as.data.frame(t(proteomics_df3))
proteomics_df3 <- cbind(proteomics_df3, proteomics_df2)
proteomics_df3 <- aggregate(proteomics_df3,list(Visit=proteomics_df3$week), mean)
proteomics_df3 <- proteomics_df3[ , -2]
proteomics_df3 <- proteomics_df3[match(order, proteomics_df3$Visit),]
proteomics_df3 <- na.omit(proteomics_df3)
rownames(proteomics_df3) <- proteomics_df3$Visit
proteomics_df3<- proteomics_df3[,-1]


save(proteomics_metadata, proteomics_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggregate_",fiber_subset,"_proteomics_df.RData",sep=""))
}


} # end of fiber_subset loop

#################
####Utils########
#################
#finding wthere the are problem values:
# a <- apply(o, 1, function(x) any(is.na(x) | is.infinite(x)))