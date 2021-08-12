rna_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/DeseqNormRNAwithMetadata.txt", sep="\t", row.names = 1, stringsAsFactors=FALSE)

rna_metadata_rows = 64
rna_df <- data.frame(lapply(rna_df, as.character), stringsAsFactors=FALSE, row.names = rownames(rna_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
rna_df <- data.frame(t(rna_df))
fiber_subset = "Arabinoxylan"
rna_df <- rna_df[which(rna_df$fiber==fiber_subset),] #subsetting by fiber
rna_df <- data.frame(lapply(rna_df, as.factor), stringsAsFactors=FALSE, row.names = rownames(rna_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
rna_df <- data.frame(t(rna_df))
rna_metadata = head(rna_df,rna_metadata_rows)
rna_df <- data.frame(t(rna_df))
save(rna_metadata, rna_df, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Data/Ternary/",fiber_subset,"_rna_df.RData",sep=""))

for (i in colnames(rna_metadata)) {
df = c(NA,NA,NA,NA)
library(lme4)
#iterate over the columns of the data that we're interested in using
for (i in (rna_metadata_rows+1):ncol(rna_df)) { 
  #print(shapiro.test(as.numeric(rna_df[,i])))
  #qqnorm(as.numeric(rna_df[,i])
 to_bind <- tryCatch({ 
  model <- glmer(as.numeric(rna_df[,i]) ~ (1|week) + (1|participant) + rna_df[[i]], data=rna_df)
  variance <- as.data.frame(VarCorr(model))
  to_bind <- variance$vcov/sum(variance$vcov)*100
 },
 error = function(err){c(NA,NA,NA)})
 to_bind <- c(to_bind,colnames(rna_df[i]))
 print(to_bind)
 df <- rbind(df,to_bind)
} 


#
write.table(df, file = paste("/Users/SLancaster/Desktop/Ternary/RNAAxWeekPart",i,".txt",sep=), sep="\t")
#df <- read.table("/Users/SLancaster/Desktop/varianceAxHomaPartWeek.txt", row.names=NULL)
#df <- df[,-1] #Only to be used when reading a df from a table
#df <- as.matrix(df) #Only to be used when reading a df from a table
df <- na.omit(df)
rownames(df) <- df[,4]
df <- df[,-4]
class(df) <- "numeric"
colnames(df) <- c("Participant", "Week", "Other")
df <- data.frame(df)
library(ggtern)

pdf(paste("/Users/SLancaster/Desktop/Ternary/RNAAxWeekPart",i,".txt",sep=), width=14)
ggtern(data = df,aes(x = Week ,y = Participant, z = Other)) +
  #geom_density_tern(aes(color=..level..), bins=4500) +
  #geom_density_tern(aes(fill=..level.., alpha=abs(..level..)),bins=500, binwidth=100) +
  #geom_point(size=0.5, alpha=0.5) +
  stat_density_tern(geom="polygon", n=4000, bins=150, aes(fill=..level..)) +
  scale_fill_gradient(low="yellow",high="red")   + 
  scale_color_gradient(low="yellow",high="red")   +
  theme_bw()  +
  theme(legend.justification=c(0,1), legend.position=c(0,1), panel.background = element_rect(fill = "light gray", colour = "light gray")) + 
  theme_nogrid() +
  guides(fill = guide_colorbar(order=1),
         alpha= guide_legend(order=2),
         color="none") + 
  labs(  title= "RNA",
         fill = "Value, V",alpha="|V - 0|")
dev.off()

}
