#!/usr/bin/R
#This code will be to split the metaphlan output into the different taxanomic level
#In the metaphlan documentation you can only do this when initially running metaphlan and I don't think there are tools to run it afterwards
#So this will  bring it down to whatever level we like
#Mike likes it at the tax levels that he known what to work with

source("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Software/Utils.r")


bugs_list_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat.pcl", sep="\t", header=TRUE, row.names = 1)
metaphlan_metadata_rows = 5

metaphlan_metadata = head(bugs_list_df,metaphlan_metadata_rows)
bugs_list_df2 <- tail(bugs_list_df, -metaphlan_metadata_rows) #Now getting rid of the metadata
bugs_list_df2 <- as.matrix(bugs_list_df2)
#bugs_list_df2 <- bugs_list_df2[!rownames(bugs_list_df2) %in% "k__Archaea.p__Euryarchaeota.c__Methanobacteria", ]
class(bugs_list_df2) <- "numeric"
#bugs_list_df2 <- bugs_list_df2/100 #This is old code. I used to do this before taking the log2
bugs_list_df2 <- data.frame(bugs_list_df2)

taxa <- rownames(bugs_list_df2)

taxa2 <- strsplit(taxa,"\\|")

kingdom <- bugs_list_df2[FALSE,]
phylum <- bugs_list_df2[FALSE,]
class_ <- bugs_list_df2[FALSE,]
order_ <- bugs_list_df2[FALSE,]
family_ <- bugs_list_df2[FALSE,]
genus <- bugs_list_df2[FALSE,]
species <- bugs_list_df2[FALSE,]
strain <- bugs_list_df2[FALSE,]

for (i in 1:length(taxa2)){
  if (length(taxa2[[i]]) == 1)
  {
    kingdom <- rbind(kingdom,bugs_list_df2[i,])
    rownames(kingdom)[length(rownames(kingdom))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 2)
  {
    phylum <- rbind(phylum,bugs_list_df2[i,])
    rownames(phylum)[length(rownames(phylum))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 3)
  {
    class_ <- rbind(class_,bugs_list_df2[i,])
    rownames(class_)[length(rownames(class_))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 4)
  {
    order_ <- rbind(order_,bugs_list_df2[i,])
    rownames(order_)[length(rownames(order_))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 5)
  {
    family_ <- rbind(family_,bugs_list_df2[i,])
    rownames(family_)[length(rownames(family_))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 6)
  {
    genus <- rbind(genus,bugs_list_df2[i,])
    rownames(genus)[length(rownames(genus))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 7)
  {
    species <- rbind(species,bugs_list_df2[i,])
    rownames(species)[length(rownames(species))] <- taxa2[[i]][length(taxa2[[i]])]
  }
  if (length(taxa2[[i]]) == 8)
  {
    strain <- rbind(strain,bugs_list_df2[i,])
    rownames(strain)[length(rownames(strain))] <- taxa2[[i]][length(taxa2[[i]])]
  }
}



kingdom <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(kingdom)))
phylum <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(phylum)))
order_ <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(order_)))
class_ <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(class_)))
family_ <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(family_)))
genus <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(genus)))
species <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(species)))
strain <- data.frame(rbind(as.matrix(metaphlan_metadata), as.matrix(strain)))

write.table(kingdom, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-kingdom.pcl",sep="\t")
write.table(phylum, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-phylum.pcl",sep="\t")
write.table(order_, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-order.pcl",sep="\t")
write.table(class_, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-class.pcl",sep="\t")
write.table(family_, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-family.pcl",sep="\t")
write.table(genus, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-genus.pcl",sep="\t")
write.table(species, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-species.pcl",sep="\t")
write.table(strain, "/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-strain.pcl",sep="\t")

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  print(fiber_subset)
  
  if (fiber_subset == "Mix") {
    sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
  } else {
    sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
  }

for (i in c("kingdom","phylum","order","class","family","genus","species","strain")){
  print(i)
  
  bugs_list_df = read.csv(paste("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/bugs-list-final-combat-",i,".pcl",sep=""), sep="\t", header=TRUE, row.names = 1)
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
  
  save(metaphlan_metadata, log_metaphlan, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/Full_Log_",fiber_subset,i,"_metaphlan_df.RData",sep=""))
  
  #Aggregate#
  metaphlan_df2 <- as.matrix(metaphlan_df2)
  class(metaphlan_df2) <- "numeric"
  metaphlan_df3 <- as.data.frame(metaphlan_metadata["Visit",])
  metaphlan_df3 <- as.data.frame(t(metaphlan_df3))
  metaphlan_df3 <- cbind(metaphlan_df3, metaphlan_df2)
  metaphlan_df3 <- aggregate(metaphlan_df3,list(Visit=metaphlan_df3$Visit), mean)
  metaphlan_df3 <- metaphlan_df3[ , -2]
  metaphlan_df3 <- metaphlan_df3[match(sorted_order, metaphlan_df3$Visit),]
  metaphlan_df3 <- na.omit(metaphlan_df3)
  rownames(metaphlan_df3) <- metaphlan_df3$Visit
  metaphlan_df3 <- metaphlan_df3[,-1]
  
  save(metaphlan_metadata, metaphlan_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/Aggregate_",fiber_subset,i,"_metaphlan_df.RData",sep=""))
  
  
} }

finding_baseline_means <- function(taxonomic_data) {
metaphlan_metadata_rows = 5
metaphlan_metadata = data.frame(head(taxonomic_data,metaphlan_metadata_rows))
bugs_list_df <- tail(taxonomic_data, -metaphlan_metadata_rows) #Now getting rid of the metadata
bugs_list_df <- as.matrix(bugs_list_df)
class(bugs_list_df) <- "numeric"
#bugs_list_df <- bugs_list_df/100

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
}

baseline_means <- finding_baseline_means(kingdom)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/kingdom_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(phylum)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/phylum_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(order_)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/order_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(class_)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/class_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(family_)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/family_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(genus)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/genus_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(species)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/species_baselines.RData",sep=""))
baseline_means <- finding_baseline_means(strain)
save(baseline_means, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/strain_baselines.RData",sep=""))


for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  print(fiber_subset)
  
  if (fiber_subset == "Mix") {
    sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
  } else {
    sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
  }
  
  for (i in c("kingdom","phylum","order","class","family","genus","species","strain")){
    print(i)
    
    #the baseline ones are with just the genus, and the metaphlan_df has the entire bugs list.
    load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/",i,"_baselines.RData",sep=""))
    load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/Full_Log_",fiber_subset,i,"_metaphlan_df.RData",sep=""))
    
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
    for (j in 1:length(baseline_mean_row)){
      normalized_row <- as.numeric(log_metaphlan[j,2:ncol(log_metaphlan)])-as.numeric(baseline_means[baseline_mean_row[j],3:ncol(baseline_means)])
      normalized_metaphlan_df <- rbind(normalized_metaphlan_df, normalized_row)
    }
    
    rownames(normalized_metaphlan_df) <- rownames(log_metaphlan)
    colnames(normalized_metaphlan_df) <- colnames(log_metaphlan[2:ncol(log_metaphlan)])
    
    save(normalized_metaphlan_df, metaphlan_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/Normalized_Log_",fiber_subset,i,"_metaphlan_df.RData",sep=""))
    
    #logaggregnorm#
    
    metaphlan_df3<- as.data.frame(metaphlan_metadata$Visit)
    metaphlan_df3<- cbind(metaphlan_df3, normalized_metaphlan_df)
    metaphlan_df3<- aggregate(metaphlan_df3,list(Visit=metaphlan_df3[,1]), mean)
    metaphlan_df3<- metaphlan_df3[ , -2]
    metaphlan_df3<- metaphlan_df3[match(sorted_order, metaphlan_df3$Visit),]
    rownames(metaphlan_df3) <- metaphlan_df3$Visit
    metaphlan_df3<- t(metaphlan_df3)
    metaphlan_df3<- na.omit(metaphlan_df3)
    metaphlan_df3<- data.frame(t(metaphlan_df3))
    metaphlan_df3 <- metaphlan_df3[,-1]
    logaggregnorm_metaphlan <- metaphlan_df3
    
    save(logaggregnorm_metaphlan, metaphlan_metadata, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/NormAggreg_Log_",fiber_subset,i,"_metaphlan_df.RData",sep=""))
    
  }
}
