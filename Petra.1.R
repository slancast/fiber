#Code snippet for Petra
#Iterating over columns to perform tests, adjusting for FDR, and saving RData

clinical_pvalues <- c()
for (i in colnames(normalized_clinical_df)[3:ncol(normalized_clinical_df)]) {
  print(i)
  female <- normalized_clinical_df[which(normalized_clinical_df2$interaction %in% c("F.30")), i]
  male <- normalized_clinical_df[which(normalized_clinical_df2$interaction %in% c("M.30")), i]
  test <- t.test(female, male)
  print(test$p.value)
  row <- c(mean(female),mean(male),sd(female),sd(male),test$p.value)
  clinical_pvalues <- rbind(clinical_pvalues,row)
}

rownames(clinical_pvalues) <- colnames(normalized_clinical_df)[3:ncol(normalized_clinical_df)]
adjusted <- p.adjust(as.numeric(clinical_pvalues[,3]), method = "BH", n = nrow(clinical_pvalues))
clinical_pvalues <- cbind(clinical_pvalues, adjusted)
colnames(clinical_pvalues) <- c("female_mean","male_mean","female_sd","male_sd","pvalue","BH_fdr")
save(clinical_pvalues, file="/Users/SLancaster/Desktop/R01_Admin_Stuff/clinical_sex_differences.RData")

