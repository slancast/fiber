#This will take the baseline means and subtract them from the data files in sorted_order to find the adjusted 
#values.
#
# 
for (fiber_subset in c("Arabinoxylan", "LCInulin", "Mix")) {
print(fiber_subset)
if (fiber_subset == "Mix") {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}

##########################################
#################LIPIDS###################
##########################################
#if (FALSE) {
#For the lipids and metabolites they used a slightly different format for WashoutFinal vs Washout_Final
if (fiber_subset == "Mix") {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","Washout_Final")
}

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/lipids_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_lipids_df.RData",sep=""))

lipids_metadata <- data.frame(t(lipids_metadata))
subject_id <- lipids_metadata$subject_id
log_lipids <- cbind(subject_id,data.frame(t(log_lipids)))

if (length(colnames(log_lipids[,2:ncol(log_lipids)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_lipids, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

log_lipids <- data.frame(lapply(log_lipids, as.character), stringsAsFactors=FALSE, row.names = rownames(log_lipids)) #To change to a double, this needs to go through character first
log_lipids <- data.frame(lapply(log_lipids, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_lipids)) #then numeric

normalized_lipids_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_lipids[i,2:ncol(log_lipids)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_lipids_df <- rbind(normalized_lipids_df, normalized_row)
}

normalized_lipids_df <- data.frame(normalized_lipids_df)
rownames(normalized_lipids_df) <- rownames(log_lipids)
colnames(normalized_lipids_df) <- colnames(log_lipids[2:ncol(log_lipids)])

save(normalized_lipids_df, lipids_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_lipids_df.RData",sep=""))

#logaggregnorm#

lipids_df3<- as.data.frame(lipids_metadata$timepoint)
lipids_df3<- cbind(lipids_df3, normalized_lipids_df)
lipids_df3<- aggregate(lipids_df3,list(Visit=lipids_df3[,1]), mean)
lipids_df3<- lipids_df3[ , -2]
lipids_df3<- lipids_df3[match(sorted_order, lipids_df3$Visit),]
rownames(lipids_df3) <- lipids_df3$Visit
lipids_df3<- t(lipids_df3)
lipids_df3<- na.omit(lipids_df3)
lipids_df3<- data.frame(t(lipids_df3))
lipids_df3 <- lipids_df3[,-1]
logaggregnorm_lipids <- lipids_df3

save(logaggregnorm_lipids, lipids_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_lipids_df.RData",sep=""))

#}

#######################################
#################RNA###################
#######################################

#if (FALSE) {
if (fiber_subset == "Mix") {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}


load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/rna_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_rna_df.RData",sep=""))

rna_metadata <- data.frame(t(rna_metadata))
participant <- rna_metadata$Participant
log_rna <- cbind(participant,log_rna)

if (length(colnames(log_rna[,2:ncol(log_rna)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_rna, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

log_rna <- data.frame(lapply(log_rna, as.character), stringsAsFactors=FALSE, row.names = rownames(log_rna)) #To change to a double, this needs to go through character first
log_rna <- data.frame(lapply(log_rna, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_rna)) #then numeric

normalized_rna_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_rna[i,2:ncol(log_rna)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_rna_df <- rbind(normalized_rna_df, normalized_row)
}

rownames(normalized_rna_df) <- rownames(log_rna)
colnames(normalized_rna_df) <- colnames(log_rna[2:ncol(log_rna)])

save(normalized_rna_df, rna_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_rna_df.RData",sep=""))

#logaggregnorm#

normalized_rna_df<- as.data.frame(rna_metadata$Visit)
normalized_rna_df<- cbind(normalized_rna_df, log_rna[,-1])
normalized_rna_df<- aggregate(normalized_rna_df,list(Visit=normalized_rna_df[,1]), mean)
normalized_rna_df<- normalized_rna_df[ , -2]
normalized_rna_df<- normalized_rna_df[match(sorted_order, normalized_rna_df$Visit),]
rownames(normalized_rna_df) <- as.character(normalized_rna_df$Visit)
normalized_rna_df <- normalized_rna_df[,-2]
normalized_rna_df<- t(normalized_rna_df)
normalized_rna_df<- na.omit(normalized_rna_df)
normalized_rna_df<- data.frame(t(normalized_rna_df))
normalized_rna_df <- normalized_rna_df[,-1]
logaggregnorm_rna <- normalized_rna_df

save(logaggregnorm_rna, rna_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_rna_df.RData",sep=""))
#}
  
  #################################################
  #################pathcoverage###################
  #################################################
if (fiber_subset == "Mix") {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}


  if (FALSE) {
  #I think these are different because there are different bugs list files that I used
  #the baseline ones are with just the genus, and the pcl_df has the entire bugs list.
  load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/pathcoverage_baselines.RData")
  load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_pathcoverage_df.RData",sep=""))
  
  pathcoverage_metadata <- data.frame(t(pathcoverage_metadata))
  #pathcoverage_metadata <- pathcoverage_metadata[which(pathcoverage_metadata$Fiber==fiber_subset),] #subsetting by fiber
  
  participant <- pathcoverage_metadata$Participant
  log_pathcoverage <- data.frame(cbind(as.character(participant),log_pathcoverage))
  
  if (length(colnames(log_pathcoverage[,2:ncol(log_pathcoverage)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
    print("Warning column names are not the same size")
  }
  
  baseline_mean_row <- apply(log_pathcoverage, 1, function(x) match(x[1],baseline_means[,1]))
 #baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)
  
  log_pathcoverage <- data.frame(lapply(log_pathcoverage, as.character), stringsAsFactors=FALSE, row.names = rownames(log_pathcoverage)) #To change to a double, this needs to go through character first
  log_pathcoverage <- data.frame(lapply(log_pathcoverage, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_pathcoverage)) #then numeric
  
  normalized_pathcoverage_df <- c()
  for (i in 1:length(baseline_mean_row)){
    normalized_row <- as.numeric(log_pathcoverage[i,2:ncol(log_pathcoverage)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
    normalized_pathcoverage_df <- rbind(normalized_pathcoverage_df, normalized_row)
  }
  
  rownames(normalized_pathcoverage_df) <- rownames(log_pathcoverage)
  colnames(normalized_pathcoverage_df) <- colnames(log_pathcoverage[2:ncol(log_pathcoverage)])
  
  save(normalized_pathcoverage_df, pathcoverage_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_pathcoverage_df.RData",sep=""))
  
  #logaggregnorm#
  
  pathcoverage_df3<- as.data.frame(pathcoverage_metadata$Visit)
  pathcoverage_df3<- cbind(pathcoverage_df3, normalized_pathcoverage_df)
  pathcoverage_df3<- aggregate(pathcoverage_df3,list(Visit=pathcoverage_df3[,1]), mean)
  pathcoverage_df3<- pathcoverage_df3
  pathcoverage_df3<- pathcoverage_df3[match(sorted_order, pathcoverage_df3$Visit),]
  rownames(pathcoverage_df3) <- pathcoverage_df3$Visit
  pathcoverage_df3 <- pathcoverage_df3[,-2]
  pathcoverage_df3<- t(pathcoverage_df3)
  pathcoverage_df3<- na.omit(pathcoverage_df3)
  pathcoverage_df3<- data.frame(t(pathcoverage_df3))
  pathcoverage_df3 <- pathcoverage_df3[,-1]
  logaggregnorm_pathcoverage <- pathcoverage_df3
  
  save(logaggregnorm_pathcoverage, pathcoverage_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_pathcoverage_df.RData",sep=""))
  }
  
#################################################
#################pathabundance###################
#################################################
if (FALSE) {
#I think these are different because there are different bugs list files that I used
#the baseline ones are with just the genus, and the pcl_df has the entire bugs list.
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/pathabundance_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_pathabundance_df.RData",sep=""))

pathabundance_metadata <- data.frame(t(pathabundance_metadata))
#pathabundance_metadata <- pathabundance_metadata[which(pathabundance_metadata$Fiber==fiber_subset),] #subsetting by fiber

participant <- pathabundance_metadata$Participant
log_pathabundance <- data.frame(cbind(as.character(participant),log_pathabundance))

if (length(colnames(log_pathabundance[,2:ncol(log_pathabundance)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_pathabundance, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

log_pathabundance <- data.frame(lapply(log_pathabundance, as.character), stringsAsFactors=FALSE, row.names = rownames(log_pathabundance)) #To change to a double, this needs to go through character first
log_pathabundance <- data.frame(lapply(log_pathabundance, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_pathabundance)) #then numeric

normalized_pathabundance_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_pathabundance[i,2:ncol(log_pathabundance)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_pathabundance_df <- rbind(normalized_pathabundance_df, normalized_row)
}

rownames(normalized_pathabundance_df) <- rownames(log_pathabundance)
colnames(normalized_pathabundance_df) <- colnames(log_pathabundance[2:ncol(log_pathabundance)])

save(normalized_pathabundance_df, pathabundance_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_pathabundance_df.RData",sep=""))

#logaggregnorm#

pathabundance_df3<- as.data.frame(pathabundance_metadata$Visit)
pathabundance_df3<- cbind(pathabundance_df3, normalized_pathabundance_df)
pathabundance_df3<- aggregate(pathabundance_df3,list(Visit=pathabundance_df3[,1]), mean)
pathabundance_df3<- pathabundance_df3[ , -2]
pathabundance_df3<- pathabundance_df3[match(sorted_order, pathabundance_df3$Visit),]
rownames(pathabundance_df3) <- pathabundance_df3$Visit
pathabundance_df3 <- pathabundance_df3
pathabundance_df3<- t(pathabundance_df3)
pathabundance_df3<- na.omit(pathabundance_df3)
pathabundance_df3<- data.frame(t(pathabundance_df3))
pathabundance_df3 <- pathabundance_df3[,-1]
logaggregnorm_pathabundance <- pathabundance_df3

save(logaggregnorm_pathabundance, pathabundance_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_pathabundance_df.RData",sep=""))
}

#############################################
#################Metaphlan###################
#############################################
if (FALSE) {
#I think these are different because there are different bugs list files that I used
#the baseline ones are with just the genus, and the metaphlan_df has the entire bugs list.
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/metaphlan_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))

metaphlan_metadata <- data.frame(t(metaphlan_metadata))
#metaphlan_metadata <- metaphlan_metadata[which(metaphlan_metadata$Fiber==fiber_subset),] #subsetting by fiber

participant <- metaphlan_metadata$Participant
log_metaphlan <- data.frame(cbind(as.character(participant),log_metaphlan))

if (length(colnames(log_metaphlan[,2:ncol(log_metaphlan)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_metaphlan, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

log_metaphlan <- data.frame(lapply(log_metaphlan, as.character), stringsAsFactors=FALSE, row.names = rownames(log_metaphlan)) #To change to a double, this needs to go through character first
log_metaphlan <- data.frame(lapply(log_metaphlan, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_metaphlan)) #then numeric

normalized_metaphlan_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_metaphlan[i,2:ncol(log_metaphlan)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_metaphlan_df <- rbind(normalized_metaphlan_df, normalized_row)
}

rownames(normalized_metaphlan_df) <- rownames(log_metaphlan)
colnames(normalized_metaphlan_df) <- colnames(log_metaphlan[2:ncol(log_metaphlan)])

save(normalized_metaphlan_df, metaphlan_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))

#logaggregnorm#

metaphlan_df3<- as.data.frame(metaphlan_metadata$Visit)
metaphlan_df3<- cbind(metaphlan_df3, normalized_metaphlan_df)
metaphlan_df3<- aggregate(metaphlan_df3,list(Visit=metaphlan_df3[,1]), mean)
metaphlan_df3<- metaphlan_df3[ , -2]
metaphlan_df3<- metaphlan_df3[match(sorted_order, metaphlan_df3$Visit),]
rownames(metaphlan_df3) <- metaphlan_df3$Visit
metaphlan_df3 <- metaphlan_df3
metaphlan_df3<- t(metaphlan_df3)
metaphlan_df3<- na.omit(metaphlan_df3)
metaphlan_df3<- data.frame(t(metaphlan_df3))
metaphlan_df3 <- metaphlan_df3[,-1]
logaggregnorm_metaphlan <- metaphlan_df3

save(logaggregnorm_metaphlan, metaphlan_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))
}

################################################
#################Genefamilies###################
################################################
if (FALSE) {
#I think these are different because there are different bugs list files that I used
#the baseline ones are with just the genus, and the genef_df has the entire bugs list.

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/genef_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_genef_df.RData",sep=""))

genef_metadata <- data.frame(t(genef_metadata))
participant <- as.character(genef_metadata$Participant)
log_genef <- data.frame(cbind(as.character(participant),log_genef))

if (length(colnames(log_genef[,2:ncol(log_genef)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_genef, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

log_genef <- data.frame(lapply(log_genef, as.character), stringsAsFactors=FALSE, row.names = rownames(log_genef)) #To change to a double, this needs to go through character first
log_genef <- data.frame(lapply(log_genef, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_genef)) #then numeric

normalized_genef_df <- data.frame(
  for (i in 2:ncol(log_genef)){
    i=double()
  }
)

for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_genef[i,2:ncol(log_genef)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_genef_df <- rbind(normalized_genef_df, normalized_row)
}

colnames(normalized_genef_df) <- colnames(log_genef)[2:ncol(log_genef)]
rownames(normalized_genef_df) <- rownames(log_genef)


save(normalized_genef_df, genef_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_genef_df.RData",sep=""))

#logaggregnorm#

genef_df3<- as.data.frame(genef_metadata$Visit)
genef_df3<- cbind(genef_df3, normalized_genef_df)
genef_df3<- aggregate(genef_df3,list(Visit=genef_df3[,1]), mean)
genef_df3<- genef_df3[ , -2]
genef_df3<- genef_df3[match(sorted_order, genef_df3$Visit),]
rownames(genef_df3) <- genef_df3$Visit
genef_df3 <- genef_df3
genef_df3<- t(genef_df3)
genef_df3<- na.omit(genef_df3)
genef_df3<- data.frame(t(genef_df3))
genef_df3 <- genef_df3[,-1]
logaggregnorm_genef <- genef_df3


save(logaggregnorm_genef, genef_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_genef_df.RData",sep=""))
}

########################################
#############Cytokine###################
########################################
#if (FALSE) {
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/cytokine_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_cytokine_df.RData",sep=""))

cytokine_metadata <- data.frame(t(cytokine_metadata))

participant <- cytokine_metadata$Participant
log_cytokine <- cbind(participant,log_cytokine)

if (length(colnames(log_cytokine[2:length(log_cytokine)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_cytokine, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

log_cytokine <- data.frame(lapply(log_cytokine, as.character), stringsAsFactors=FALSE, row.names = rownames(log_cytokine)) #To change to a double, this needs to go through character first
log_cytokine <- data.frame(lapply(log_cytokine, as.numeric), stringsAsFactors=FALSE, row.names = rownames(log_cytokine)) #then numeric

normalized_cytokine_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_cytokine[i,2:ncol(log_cytokine)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_cytokine_df <- rbind(normalized_cytokine_df, normalized_row)
}

rownames(normalized_cytokine_df) <- rownames(log_cytokine)
colnames(normalized_cytokine_df) <- colnames(log_cytokine[2:ncol(log_cytokine)])

save(normalized_cytokine_df, cytokine_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_cytokine_df.RData",sep=""))

#logaggregnorm#

cytokine_df3<- as.data.frame(cytokine_metadata$Week)
cytokine_df3<- cbind(cytokine_df3, normalized_cytokine_df)
cytokine_df3<- aggregate(cytokine_df3,list(Visit=cytokine_df3[,1]), mean)
cytokine_df3<- cytokine_df3[ , -2]
cytokine_df3<- cytokine_df3[match(sorted_order, cytokine_df3$Visit),]
rownames(cytokine_df3) <- cytokine_df3$Visit
cytokine_df3 <- cytokine_df3
cytokine_df3<- t(cytokine_df3)
cytokine_df3<- na.omit(cytokine_df3)
cytokine_df3<- data.frame(t(cytokine_df3))
cytokine_df3 <- cytokine_df3[,-1]
logaggregnorm_cytokine <- cytokine_df3

save(logaggregnorm_cytokine, cytokine_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_cytokine_df.RData",sep=""))
#}

#############################################
#################Metabolome##################
#############################################
#if (FALSE) {
if (fiber_subset == "Mix") {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","Washout_Final")
}

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/metabolomic_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))

metabolomics_metadata <- data.frame(t(metabolomics_metadata))

participant <- data.frame(metabolomics_metadata$subject_id)
log_metabolomics <- data.frame(log_metabolomics)
log_metabolomics<- cbind(participant,log_metabolomics)

metabolomic_rownames <- rownames(log_metabolomics)
log_metabolomics <- data.frame(lapply(log_metabolomics, as.character), stringsAsFactors=FALSE) #To change to a double, this needs to go through character first
log_metabolomics <- data.frame(lapply(log_metabolomics, as.numeric), stringsAsFactors=FALSE) #then numeric
rownames(log_metabolomics) <- metabolomic_rownames

baseline_mean_header <- baseline_means[,1:2]
baseline_means <- baseline_means[,colnames(baseline_means) %in% colnames(log_metabolomics)]
baseline_means <- cbind(baseline_mean_header, baseline_means)
if (length(colnames(log_metabolomics[,2:ncol(log_metabolomics)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_metabolomics, 1, function(x) match(as.character(x[1]),trimws(as.character(baseline_means[,1]))))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

normalized_metabolomics_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_metabolomics[i,2:ncol(log_metabolomics)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_metabolomics_df <- rbind(normalized_metabolomics_df, normalized_row)
}

rownames(normalized_metabolomics_df) <- rownames(log_metabolomics)
colnames(normalized_metabolomics_df) <- colnames(log_metabolomics[2:ncol(log_metabolomics)])

save(normalized_metabolomics_df, metabolomics_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))

#logaggregnorm#

metabolomics_df3<- as.data.frame(metabolomics_metadata$timepoint)
metabolomics_df3<- cbind(metabolomics_df3, normalized_metabolomics_df)
metabolomics_df3<- aggregate(metabolomics_df3,list(Visit=metabolomics_df3[,1]), mean)
metabolomics_df3<- metabolomics_df3[ , -2]
metabolomics_df3<- metabolomics_df3[match(sorted_order, metabolomics_df3$Visit),]
rownames(metabolomics_df3) <- metabolomics_df3$Visit
metabolomics_df3 <- metabolomics_df3
metabolomics_df3<- t(metabolomics_df3)
metabolomics_df3<- na.omit(metabolomics_df3)
metabolomics_df3<- data.frame(t(metabolomics_df3))
metabolomics_df3 <- metabolomics_df3[,-1]
logaggregnorm_metabolomics <- metabolomics_df3

save(logaggregnorm_metabolomics, metabolomics_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))
#}

#############################################
##################Proteome###################
#############################################
if (FALSE) {
if (fiber_subset == "Mix") {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/proteomic_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_proteomics_df.RData",sep=""))

proteomics_metadata <- data.frame(t(proteomics_metadata))
participant <- data.frame(proteomics_metadata$participant)
log_proteomics <- data.frame(log_proteomics)
log_proteomics<- cbind(participant,log_proteomics)

proteomic_rownames <- rownames(log_proteomics)
log_proteomics <- data.frame(lapply(log_proteomics, as.character), stringsAsFactors=FALSE) #To change to a double, this needs to go through character first
log_proteomics <- data.frame(lapply(log_proteomics, as.numeric), stringsAsFactors=FALSE) #then numeric
rownames(log_proteomics) <- proteomic_rownames

baseline_mean_header <- baseline_means[,1:2]
baseline_means <- baseline_means[,colnames(baseline_means) %in% colnames(log_proteomics)]
baseline_means <- cbind(baseline_mean_header, baseline_means)
if (length(colnames(log_proteomics[,2:ncol(log_proteomics)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
print("Warning column names are not the same size")
}

baseline_mean_row <- apply(log_proteomics, 1, function(x) match(as.character(x[1]),trimws(as.character(baseline_means[,1]))))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)

normalized_proteomics_df <- c()
for (i in 1:length(baseline_mean_row)){
normalized_row <- as.numeric(log_proteomics[i,2:ncol(log_proteomics)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
normalized_proteomics_df <- rbind(normalized_proteomics_df, normalized_row)
}

rownames(normalized_proteomics_df) <- rownames(log_proteomics)
colnames(normalized_proteomics_df) <- colnames(log_proteomics[2:ncol(log_proteomics)])

save(normalized_proteomics_df, proteomics_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_proteomics_df.RData",sep=""))

#logaggregnorm#

proteomics_df3<- as.data.frame(proteomics_metadata$week)
proteomics_df3<- cbind(proteomics_df3, normalized_proteomics_df)
proteomics_df3<- aggregate(proteomics_df3,list(Visit=proteomics_df3[,1]), mean)
proteomics_df3<- proteomics_df3[ , -2]
proteomics_df3<- proteomics_df3[match(sorted_order, proteomics_df3$Visit),]
rownames(proteomics_df3) <- proteomics_df3$Visit
proteomics_df3 <- proteomics_df3
proteomics_df3 <- proteomics_df3[,-1]
logaggregnorm_proteomics <- proteomics_df3

save(logaggregnorm_proteomics, proteomics_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_",fiber_subset,"_proteomics_df.RData",sep=""))
#}


#############################################
##################Clinicals##################
#############################################
#All the NAs in the clinicals makes this part difficult
#I'm going to back burner this until people start clamoring for it

if (FALSE) {
if (fiber_subset == "Mix") {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/clinical_baselines.RData")
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_clinical_df.RData",sep=""))

baseline_colnames <- colnames(baseline_means)

clinical_metadata <- data.frame(t(clinical_metadata))
participant <- data.frame(clinical_metadata$participant)
clinical_df2 <- data.frame(clinical_df2)
clinical_df2<- cbind(participant,clinical_df2)

clinical_rownames <- rownames(clinical_df2)
clinical_df2 <- data.frame(lapply(clinical_df2, as.character), stringsAsFactors=FALSE) #To change to a double, this needs to go through character first
clinical_df2 <- data.frame(lapply(clinical_df2, as.numeric), stringsAsFactors=FALSE) #then numeric
rownames(clinical_df2) <- clinical_rownames

baseline_mean_header <- baseline_means[,1:2]
baseline_means <- baseline_means[,colnames(baseline_means) %in% colnames(clinical_df2)]
baseline_means <- cbind(baseline_mean_header, baseline_means)

columns <- match(colnames(baseline_means[3:length(baseline_means)]),colnames(clinical_df2))
columns <- c(1,columns)
clinical_df2.a <- clinical_df2[,columns]

if (length(colnames(clinical_df2.a[,2:ncol(clinical_df2.a)])) != length(colnames(baseline_means[3:length(baseline_means)]))) {
  print("Warning column names are not the same size")
}

baseline_mean_row <- apply(clinical_df2.a, 1, function(x) match(as.character(x[1]),trimws(as.character(baseline_means[,1]))))

baseline_means <- as.matrix(baseline_means)
normalized_clinicals_df <- c()
for (i in 1:length(baseline_mean_row)){
  normalized_row <- as.numeric(clinical_df2.a[i,2:ncol(clinical_df2.a)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
  normalized_clinicals_df <- rbind(normalized_clinicals_df, normalized_row)
}

rownames(normalized_clinicals_df) <- clinical_rownames
colnames(normalized_clinicals_df) <- colnames(clinical_df2.a[2:ncol(clinical_df2.a)])
normalized_clinicals_df <- data.frame(normalized_clinicals_df)

save(normalized_clinicals_df, clinical_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_clinicals_df.RData",sep=""))
# 
# load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_clinical_df.RData",sep=""))
# 
# class(baseline_means) <- "numeric"
# colnames(baseline_means) <- baseline_colnames
# #baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)]+1,2)
# 
# clinical_metadata <- data.frame(t(clinical_metadata))
# participant <- data.frame(clinical_metadata$participant)
# log_clinicals <- data.frame(log_clinicals)
# log_clinicals<- cbind(participant,log_clinicals)
# 
# clinical_rownames <- rownames(log_clinicals)
# log_clinicals <- data.frame(lapply(log_clinicals, as.character), stringsAsFactors=FALSE) #To change to a double, this needs to go through character first
# log_clinicals <- data.frame(lapply(log_clinicals, as.numeric), stringsAsFactors=FALSE) #then numeric
# rownames(log_clinicals) <- clinical_rownames
# 
# columns <- match(colnames(baseline_means[,3:ncol(baseline_means)]),colnames(log_clinicals))
# columns <- c(1,columns)
# log_clinicals.a <- log_clinicals[,columns]
# colnames(log_clinicals.a) <- colnames(log_clinicals)[columns]
#   
# if (length(colnames(log_clinicals.a[,2:ncol(log_clinicals.a)])) != length(colnames(baseline_means[,3:ncol(baseline_means)]))) {
#   print("Warning column names are not the same size")
# }
# 
# baseline_mean_row <- apply(log_clinicals.a, 1, function(x) match(as.character(x[1]),trimws(as.character(baseline_means[,1]))))
# 
# baseline_means <- as.matrix(baseline_means)
# normalized_clinicals_df <- c()
# for (i in 1:length(baseline_mean_row)){
#   normalized_row <- as.numeric(log_clinicals.a[i,2:ncol(log_clinicals.a)])-as.numeric(baseline_means[baseline_mean_row[i],3:ncol(baseline_means)])
#   normalized_clinicals_df <- rbind(normalized_clinicals_df, normalized_row)
# }
# 
# rownames(normalized_clinicals_df) <- clinical_rownames
# colnames(normalized_clinicals_df) <- colnames(log_clinicals.a[2:ncol(log_clinicals.a)])
# 
# save(normalized_clinicals_df, clinical_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_",fiber_subset,"_clinicals_df.RData",sep=""))
# 

#logaggregnorm#

clinical_df3<- as.data.frame(clinical_metadata$week)
clinical_df3<- cbind(clinical_df3, normalized_clinicals_df[,-1])
clinical_df3<- aggregate(clinical_df3,list(Visit=clinical_df3[,1]), mean)
clinical_df3<- clinical_df3[ , -2]
clinical_df3<- clinical_df3[match(sorted_order, clinical_df3$Visit),]
rownames(clinical_df3) <- clinical_df3$Visit
clinical_df3 <- clinical_df3[,-2]
clinical_df3<- t(clinical_df3)
clinical_df3<- na.omit(clinical_df3)
clinical_df3<- data.frame(t(clinical_df3))
clinical_df3 <- clinical_df3[,-1]
aggregnorm_clinicals <- clinical_df3

save(aggregnorm_clinicals, clinical_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_clinicals_df.RData",sep=""))
}


}

} #End fiber subset


####Utils
# 
# aggregating_dataframes <- function(dataframe,week_metadata){
# 
# df<- as.data.frame(week_metadata)
# df<- cbind(df, dataframe)
# df<- aggregate(df,list(Visit=df[,1]), mean)
# df<- df[ , -2]
# df<- df[match(sortedorder, df$Visit),]
# rownames(df) <- df$Visit
# df <- df
# df <- df[,-1]
# 

