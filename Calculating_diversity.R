for (fiber_subset in c("Arabinoxylan", "LCInulin", "Mix")) {

load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/Normalized_Log_",fiber_subset,"species_metaphlan_df.RData",sep=""))

library(vegan)

normalized_metaphlan_df <- data.frame(t(normalized_metaphlan_df))
normalized_metaphlan_df <- 0.001 + normalized_metaphlan_df

shannon_diversity <- c()
simpson_diversity <- c()

#Generating the diversity calculations for each time point
for (i in 1:ncol(normalized_metaphlan_df)) {
  sample_diversity <- diversity(as.numeric(as.matrix(normalized_metaphlan_df[,i])))
  print(sample_diversity)
  shannon_diversity <- c(shannon_diversity, sample_diversity)
  sample_diversity <- diversity(as.numeric(as.matrix(normalized_metaphlan_df[,i])), index = "simpson")
  simpson_diversity <- c(simpson_diversity, sample_diversity)
}


if (fiber_subset == "Mix") {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")
} else {
  sorted_order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","WashoutFinal")
}

#Testing the differences in diversity statistically
dim(metaphlan_metadata)
for_manova <- cbind(shannon_diversity, simpson_diversity, metaphlan_metadata)

save(for_manova,file = paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/taxa/species_diversity_",fiber_subset,".RData",sep=""))

manova_dose <- as.factor(for_manova[which(for_manova$Dose %in% c("Baseline", "30")),]$Dose)
manova_diversity <- as.numeric(as.matrix(for_manova[which(for_manova$Dose %in% c("Baseline", "30")),]$shannon_diversity))

#manova_dose <- as.factor(for_manova$Dose)
#manova_diversity <- as.numeric(as.matrix(for_manova$shannon_diversity))

fit <- anosim(manova_diversity, manova_dose)

plot(manova_dose, manova_diversity)

}



