library(Rtsne)
time_point <- "Baseline"

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Arabinoxylan_metabolome_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_Arabinoxylan_metabolomics_df.RData")
merged_metabolomics_df <- normalized_metabolomics_df
merged_metabolomics_metadata <- t(metabolomics_metadata)
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_LCInulin_metabolome_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_LCInulin_metabolomics_df.RData")
merged_metabolomics_df <- rbind(merged_metabolomics_df, normalized_metabolomics_df)
merged_metabolomics_metadata <- rbind(merged_metabolomics_metadata, t(metabolomics_metadata))
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Mix_metabolome_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_Mix_metabolomics_df.RData")
merged_metabolomics_df <- rbind(merged_metabolomics_df, normalized_metabolomics_df)
merged_metabolomics_metadata <- rbind(merged_metabolomics_metadata, t(metabolomics_metadata))
merged_metabolomics_metadata <- data.frame(merged_metabolomics_metadata)

merged_metabolomics_df <- t(merged_metabolomics_df)
merged_metabolomics_df <- na.omit(merged_metabolomics_df)
merged_metabolomics_df <- data.frame(t(merged_metabolomics_df))
rownames(merged_metabolomics_df) <- make.names(paste(rownames(merged_metabolomics_df),merged_metabolomics_metadata$timepoint,merged_metabolomics_metadata$fiber,sep=""), unique=TRUE)
merged_metabolomics_df <- cbind(merged_metabolomics_df, merged_metabolomics_metadata$subject_id,merged_metabolomics_metadata$timepoint)
merged_metabolomics_df <- data.frame(merged_metabolomics_df)
merged_metabolomics_df  <- merged_metabolomics_df[which(merged_metabolomics_df$merged_metabolomics_metadata.timepoint==time_point),]


tsne_metabolomics <- Rtsne(merged_metabolomics_df, perplexity=5)
tsne_plots <- data.frame(tsne_metabolomics$Y)

library(randomcoloR)
subject_ids <- unique(merged_metabolomics_df$merged_metabolomics_metadata.subject_id)
palette <- distinctColorPalette(length(subject_ids), altCol=TRUE) 
subject_colors <- cbind(data.frame(subject_ids), palette)
rownames(subject_colors) <- subject_colors$subject_ids

merged_metabolomics_metadata  <- merged_metabolomics_metadata[which(merged_metabolomics_metadata$timepoint==time_point),]
meatdata_colors <- c()
for (i in 1:nrow(merged_metabolomics_metadata)) {
  print(i)
  subject_for_appending <- merged_metabolomics_metadata[i,"subject_id"]
  print(subject_for_appending)
  meatdata_colors <- c(meatdata_colors, as.character(subject_colors[subject_for_appending,2]))
}

tsne_plots <- cbind(tsne_plots, meatdata_colors)
library(ggplot2)
ggplot(data=tsne_plots, aes(x=as.numeric(tsne_plots[,1]), y=as.numeric(tsne_plots[,2]))) +
  geom_point(size=2, aes(col=tsne_plots$meatdata_colors)) +
  scale_colour_manual(values = palette) +
  theme(legend.position = "none") +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20), legend.position = "none") 

