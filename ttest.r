# I will follow the outline of the Effect_Size.r program
# However, I will instead use the unnormalized values to determine differences between
# baseline and any timepoint like 30g.

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_LCInulin_cytokine_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_LCInulin_cytokine_df.RData")

cytokine_metadata2 <- data.frame(t(cytokine_metadata))
normalized_cytokine_df2 <- cbind(cytokine_metadata2$Week, data.frame(normalized_cytokine_df))
colnames(normalized_cytokine_df2)[1] <- "Week"

cytokine_pvalues <- c()
for (i in colnames(normalized_cytokine_df2)[2:ncol(normalized_cytokine_df2)]) {
  print(i)
  baselines <- normalized_cytokine_df2[which(normalized_cytokine_df2$Week %in% c("20")), i]
  fiber <- normalized_cytokine_df2[which(normalized_cytokine_df2$Week %in% c("30")), i]
  test <- t.test(baselines, fiber)
  print(test$p.value)
  row <- c(mean(baselines),mean(fiber),test$p.value)
  cytokine_pvalues <- rbind(cytokine_pvalues,row)
}

rownames(cytokine_pvalues) <- colnames(normalized_cytokine_df2)[2:ncol(normalized_cytokine_df2)]
adjusted <- p.adjust(as.numeric(cytokine_pvalues[,3]), method = "BH", n = nrow(cytokine_pvalues))
cytokine_pvalues <- cbind(cytokine_pvalues, adjusted)
colnames(cytokine_pvalues) <- c("baselines_mean","fiber_mean","pvalue","BH_fdr")

# Now that I've tested the individual cytokines, I want to make a test of all the cytokines to each other.

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/NormAggreg_Log_Arabinoxylan_cytokine_df.RData")

aggregated <- data.frame(t(logaggregnorm_cytokine))
test <- t.test(aggregated[,2], aggregated[,3])

