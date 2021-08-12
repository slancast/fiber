#!/urs/bin/R
#This function will be to break up a particular analyte into responders and non-responders
#It will add up the values at 10g, 20g, and 30g.
#Then will break up participants into thirds. Top third, second third, and final third.
#Top third will be responders, and final thrid will be non-responders. Middle third will be considered
#Some moderate responders, or an in between category.
#
#The data frame needs to have the structure where the columns are the analytes
#Initially I will toss the time points that are not 10, 20 and 30 and then sum by participant.

data_frame <- data.frame(normalized_clinicals_df)
analyte <- "LDL..Calculated."
metadata <- clinical_metadata

responders_nonresponders <-  function(data_frame, analyte, metadata){
  #data frame is the data frame
  #analyte is the string name of the analyte
  #metadata is the metadata matrix
  analyte_data <- data_frame[,analyte]
  analyte_data2 <- cbind(as.character(metadata[,"week"]),as.character(metadata[,"participant"]),analyte_data )
  analyte_data3 <- analyte_data2[which(analyte_data2[,1] %in% c("10","20","30")),]
  class(analyte_data3) <- "numeric"
  analyte_data4 <- aggregate(analyte_data3,list(Participant = analyte_data3[,2]), mean )
  analyte_data5 <- analyte_data4[order(analyte_data4[,"analyte_data"]),]
  responders <- analyte_data5[,c(1,4)]
  colnames(responders)[2] <- analyte
  
}

#This funciton will take the responders and the non-responders
#it will the iterate over a data frame with all the baselines.
#I will do it with the baseline means because I already have that dataset. It might also
#be the correct way. I worry that three for each participant will confound the data.
#
#It will then do a ttest between responders and non-repsonders
#by first making two dataframes, one with responders and one with nonresponders,
#and then performing a ttest between them. This will be done with baseline_mean data


ttest_topvsbottom <- function(responders,baseline_means,metadata){
  top_third <- responders$Participant[1:6]
  bottom_third <- responders$Participant[13:18]
  baseline_means[,1] <- as.character(baseline_means[,1] )
  class(baseline_means[,1]) <- "numeric"
  top_third_data <- baseline_means[which(baseline_means[,1] %in% top_third),]
  bottom_third_data <- baseline_means[which(baseline_means[,1] %in% bottom_third),]
  top_third_data <- as.matrix(top_third_data[,-2])
  bottom_third_data <- as.matrix(bottom_third_data[,-2])
  class(top_third_data) <- "numeric"
  class(bottom_third_data) <- "numeric"
  
  pvalues <- data.frame(analyte=character(),pvalue=double(),top_average=double(),bottom_average=double())
  pvalues <- as.matrix(pvalues)
  for (i in 1:ncol(bottom_third_data)){
    #print(colnames(bottom_third_data)[i])
    #tested <- t.test(top_third_data[,i],bottom_third_data[,i]) #Neither pvalue histogram looks great
    tested <- wilcox.test(top_third_data[,i],bottom_third_data[,i])
    #print(tested)
    pvalues <- rbind(pvalues, c(colnames(bottom_third_data)[i], tested$p.value, mean(top_third_data[,i]), mean(bottom_third_data[,i])))
    
  }
  
  hist(as.numeric(pvalues[,2]),breaks = 20)
  
  pvalues <- pvalues[order(pvalues[,2]),]
  padjusted <- p.adjust(pvalues[,2], method="fdr")
  pvalues <- cbind(pvalues, padjusted)
  
}
  
# This funciton will do the same as above but it won't involve baseline means
# Except this time it will then see the correlates between the responders and 
# non responders from the first function but then see the differences
# between the fiber treatmetns rather than the baselines

data_frame2 <- normalized_metaphlan_df
metadata2 <- metabolomics_metadata
colnames(metadata2)
  
ttest_topvsbottom_dose <- function(responders,data_frame2,metadata2){
  top_third <- responders$Participant[1:6]
  bottom_third <- responders$Participant[13:18]

# Append dosage info to data_frame
# Number one in here needs to be changed to the appropriate column names
data_with_dose <- cbind(as.character(metadata2[,"timepoint"]), as.character(metadata2[,"subject_id"]), data_frame2)
#data <- data_with_dose[which(data_with_dose[,1] %in% c("10","20","30")), ]
data <- data_with_dose[which(data_with_dose[,1] %in% c("Baseline")), ]

#tossing dose, making participant numeric
data <- data[,-1]
data <- as.matrix(data)
class(data) <- "numeric"

#sub setting top third data and bottom third data
top_third_data <- data[which(data[,1] %in% top_third),]
bottom_third_data <- data[which(data[,1] %in% bottom_third),]

#ttest across bottom to top third
pvalues <- data.frame(analyte=character(),pvalue=double(),top_average=double(),bottom_average=double())
pvalues <- as.matrix(pvalues)
for (i in 1:ncol(bottom_third_data)){
  #print(colnames(bottom_third_data)[i])
  tested <- t.test(top_third_data[,i],bottom_third_data[,i]) #Wilcoxen looks a little better here
  #tested <- wilcox.test(top_third_data[,i],bottom_third_data[,i])
  #print(mean(top_third_data[,i]))
  #print(tested$p.value)
  pvalues <- rbind(pvalues, c(colnames(bottom_third_data)[i], tested$p.value, mean(top_third_data[,i]), mean(bottom_third_data[,i])))
}

hist(as.numeric(pvalues[,2]),breaks = 20)

pvalues <- pvalues[order(pvalues[,2]),]
padjusted <- p.adjust(pvalues[,2], method="fdr")
pvalues <- cbind(pvalues, padjusted)

pvalues <- pvalues[-1,] #Getting rid of row 1 which is the participant averages

}

metaphlan_pvalues <- pvalues
padjusted <- p.adjust(metaphlan_pvalues[,1], method = "BH", n = nrow(metaphlan_pvalues))
metaphlan_pvalues <- cbind(metaphlan_pvalues, padjusted)
metaphlan_pvalues <- metaphlan_pvalues[order(metaphlan_pvalues[,4]),]
write.table(metaphlan_pvalues, file = "/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/responders_vs_nonresponders/LDL_responders_metaphlan.csv",sep=",")
  
pathabundance_pvalues <- pvalues
padjusted <- p.adjust(pathabundance_pvalues[,1], method = "BH", n = nrow(pathabundance_pvalues))
pathabundance_pvalues <- cbind(pathabundance_pvalues, padjusted)
pathabundance_pvalues <- pathabundance_pvalues[order(pathabundance_pvalues[,4]),]
write.table(pathabundance_pvalues, file = "/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/responders_vs_nonresponders/LDL_responders_pathabundance.csv",sep=",")

clinical_pvalues <- pvalues
padjusted <- p.adjust(clinical_pvalues[,1], method = "BH", n = nrow(clinical_pvalues))
clinical_pvalues <- cbind(clinical_pvalues, padjusted)
clinical_pvalues <- clinical_pvalues[order(clinical_pvalues[,4]),]
write.table(clinical_pvalues, file = "/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/responders_vs_nonresponders/LDL_responders_clinical.csv",sep=",")

metabolomics_pvalues <- pvalues
padjusted <- p.adjust(metabolomics_pvalues[,1], method = "BH", n = nrow(metabolomics_pvalues))
metabolomics_pvalues <- cbind(metabolomics_pvalues, padjusted)
metabolomics_pvalues <- metabolomics_pvalues[order(metabolomics_pvalues[,4]),]
write.table(metabolomics_pvalues, file = "/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/responders_vs_nonresponders/LDL_responders_metabolomics.csv",sep=",")

lipids_pvalues <- pvalues
padjusted <- p.adjust(lipids_pvalues[,1], method = "BH", n = nrow(lipids_pvalues))
lipids_pvalues <- cbind(lipids_pvalues, padjusted)
lipids_pvalues <- lipids_pvalues[order(lipids_pvalues[,4]),]
write.table(lipids_pvalues, file = "/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/responders_vs_nonresponders/LDL_responders_lipids.csv",sep=",")

genef_pvalues <- pvalues
padjusted <- p.adjust(genef_pvalues[,2], method = "BH", n = nrow(genef_pvalues))
genef_pvalues <- cbind(genef_pvalues, padjusted)
genef_pvalues <- genef_pvalues[order(genef_pvalues[,5]),]
write.table(genef_pvalues, file = "/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/responders_vs_nonresponders/LDL_responders_genef.csv",sep=",", row.names = FALSE)


