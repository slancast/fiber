
for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {

if (fiber_subset == "Mix") {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}

library(openxlsx)
clinical_df <- read.xlsx("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Data/total_metadata_cleaned.xlsx", colNames = TRUE)

clinical_df$fiber <- as.character(clinical_df$fiber)
clinical_df <-  clinical_df[which( clinical_df$fiber==fiber_subset),] #subsetting by fiber
clinical_df <- data.frame(t( clinical_df))

#Now I will clean up the heterogeneous clinical data. Get rid of non-numeric characters, turn categorical variables into factors, etc.
clinical_df <- clinical_df[-which(rownames(clinical_df) %in% c("date_time_blood_draw_clinicals","Alb/Creat.Interp")),]

factors <- c("Ethnicity","IR_IS_from_SSPG")
clinical_df <- clinical_df[-which(rownames(clinical_df) %in% factors),]

#There is a problem with missing data where if there is a single value missing during the aggregation step
#the entire row becomes NA. I see two ways around this: 1) impute the values or 2) throw out rows/columns
#that contain the mising values. My intuition is that imputing the values will be a large challenge, 
#so I will begin by throwing out rows and columns that seem to be the least informative, while 
#contributing the most to the NA values.

clinical_df2 <- clinical_df[-which(rownames(clinical_df) %in% c("CREATININE,.URINE","Albumin,.Urine","Alb/Creat.Ratio")),] #Most NAs in here
clinical_df2 <- clinical_df2[-which(rownames(clinical_df2) %in% c("Insulin,.Fasting","HOMA.IR","BMI","SSPG","IR_IS_from_SSPG")),] #Second category of NAs; after this set it looks like most values are preserved

#Beginning aggregation

clinical_df2 <- data.frame(t(clinical_df2))
clinical_df2 <- na.omit(clinical_df2)
sample_names <- clinical_df2[,1]
clinical_metadata <- cbind(clinical_df2[,2:4],clinical_df2[,"Sex"])
colnames(clinical_metadata)[4] <- "Sex"

columns <- ncol(clinical_df2)
clinical_df2[,6:columns] <- gsub("[^0-9\\.]", "", as.matrix(clinical_df2[,6:columns]))
clinical_df2[,6:columns] <- as.numeric(as.matrix(clinical_df2[,6:columns])) #absolutely critical. For some fucking reason the code here turns the float numbers into strings. No idea why, but I have to change the to flaots before performing aggregate.
clinical_df3 <- aggregate(data.frame(clinical_df2),list(Visit=clinical_df2$week), FUN=mean, na.action = na.omit)
clinical_df3 <- clinical_df3[ ,-which(colnames(clinical_df3) %in% c("run","fiber","week","DOB.month.year","Sex","Ethnicity"))]
clinical_df3 <- clinical_df3[match(sorted_order,clinical_df3$Visit),]
clinical_df3 <- na.omit(t(clinical_df3))
clinical_df3 <- data.frame(t(clinical_df3))
rownames(clinical_df3) <- clinical_df3$Visit
clinical_df3 <-  clinical_df3[,-1]

Aggreg_clinicals <- clinical_df3
save(clinical_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Aggreg_",fiber_subset,"_clinical_df.RData",sep=""))

#Since I'm not taking the log of the clinicals I will just run the baseline normalization
#on them here. For the other omics measurements, for the sake of clarity, I only did normalization
#on the log datasets. Otherwise I was worried about proliferation of datasets becoming confusing,
#potentially leading to mistakes.

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/clinical_baselines.RData")

#There are some values missing between clinical df2 and baseline_means
#I think this has to do with NAs that are in the dataset where we then
#throw out the entire analyte. This is necessary for the pipeline to work.
clinical_df4 <- clinical_df2[ ,which(colnames(clinical_df2) %in% colnames(baseline_means))]
clinical_df4 <- cbind(as.character(clinical_df2$participant),clinical_df4)

if (length(colnames(clinical_df4)) != length(colnames(baseline_means))) {
  print("Warning column names are not the same size")
}

baseline_mean_row <- apply(clinical_df4, 1, function(x) match(x[1],baseline_means[,1]))
#baseline_means[,3:ncol(baseline_means)] <- log(baseline_means[,3:ncol(baseline_means)],2)

clinical_df4 <- data.frame(lapply(clinical_df4, as.character), stringsAsFactors=FALSE, row.names = rownames(clinical_df4)) #To change to a double, this needs to go through character first
clinical_df4 <- data.frame(lapply(clinical_df4, as.numeric), stringsAsFactors=FALSE, row.names = rownames(clinical_df4)) #then numeric

baseline_means <- data.frame(lapply(baseline_means, as.character), stringsAsFactors=FALSE, row.names = rownames(baseline_means)) #To change to a double, this needs to go through character first
baseline_means <- data.frame(lapply(baseline_means, as.numeric), stringsAsFactors=FALSE, row.names = rownames(baseline_means)) #then numeric

normalized_clinical_df <- c()
for (i in 1:length(baseline_mean_row)){
  normalized_row <- as.numeric(clinical_df4[i,2:ncol(clinical_df4)])-as.numeric(baseline_means[baseline_mean_row[i],2:ncol(baseline_means)])
  normalized_clinical_df <- rbind(normalized_clinical_df, normalized_row)
}

normalized_clinical_df <- data.frame(normalized_clinical_df)
normalized_clinical_df <- cbind(clinical_df2$participant,normalized_clinical_df)
rownames(normalized_clinical_df) <- rownames(clinical_df4)
colnames(normalized_clinical_df) <- colnames(clinical_df4[1:ncol(clinical_df4)])

normalized_clinical_df <- cbind(sample_names, normalized_clinical_df)
save(clinical_metadata, normalized_clinical_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_",fiber_subset,"_clinical_df.RData",sep=""))

clinical_df5 <- cbind(clinical_df2$week, normalized_clinical_df)
colnames(clinical_df5)[1] <- "week"
clinical_df5 <- aggregate(data.frame(clinical_df5),list(Visit=clinical_df5$week), FUN=mean, na.action = na.omit)
clinical_df5 <- clinical_df5[ ,-2]
clinical_df5 <- clinical_df5[ ,-2]
clinical_df5 <- clinical_df5[match(sorted_order,clinical_df5$Visit),]
clinical_df5 <- na.omit(t(clinical_df5))
clinical_df5 <- data.frame(t(clinical_df5))
rownames(clinical_df5) <- clinical_df5$Visit
clinical_df5 <-  clinical_df5[,-1]

NormAggreg_clinicals <- clinical_df5
save(clinical_df3, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_",fiber_subset,"_clinical_df.RData",sep=""))

} #ending fiber loop

