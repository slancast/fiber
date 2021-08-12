fiber_subset = "Arabinoxylan"

load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))

for_manova <- data.frame(t(metaphlan_metadata))

#manova_dose <- as.factor(for_manova[which(for_manova$Dose %in% c("Baseline", "30")),]$Dose)
#manova_diversity <- as.numeric(as.matrix(for_manova[which(for_manova$Dose %in% c("Baseline", "30")),]$simpson_diversity))

manova_dose <- as.factor(for_manova$Dose)
manova_diversity <- as.numeric(as.matrix(for_manova$shannon_diversity))

#fit <- anosim(manova_diversity, manova_dose)

plot(manova_dose, manova_diversity)
