library(Mfuzz)
pcl_df = read.csv("/Users/SLancaster/Desktop/relab-pathabundance.pcl", sep="\t", header=TRUE, row.names = 1)
bugs_list_df = read.csv("/Users/SLancaster/Desktop/bugs_list.pcl", sep="\t", header=TRUE, row.names = 1)

pcl_metadata_rows = 64
pcl_df <- data.frame(t(pcl_df))
fiber_subset = "Mix"
pcl_df <- pcl_df[which(pcl_df$Fiber==fiber_subset),] #subsetting by fiber
pcl_df <- data.frame(t(pcl_df))
pcl_metadata = data.frame(head(pcl_df,pcl_metadata_rows))
bugs_list_df <- data.frame(t(bugs_list_df))
bugs_list_df<- bugs_list_df[which(bugs_list_df$Fiber==fiber_subset),] #subsetting by fiber
bugs_list_df <- data.frame(t(bugs_list_df))
bugs_list_df2 <- tail(bugs_list_df, -pcl_metadata_rows) #Now getting rid of the metadata
bugs_list_df2 <- as.matrix(bugs_list_df2)
class(bugs_list_df2) <- "numeric"
bugs_list_df2 <- bugs_list_df2/100

pcl_df2 <- tail(pcl_df, -pcl_metadata_rows) #Now getting rid of the metadata
colnames(bugs_list_df2) <- colnames(pcl_df2)
pcl_df2 <- rbind(bugs_list_df2,pcl_df2)
pcl_df2 <- as.matrix(pcl_df2)
class(pcl_df2) <- "numeric"
pcl_df2 <- asin(pcl_df2)
pcl_df2 <- data.frame(pcl_df2)
pcl_df2 <- data.frame(lapply(pcl_df2, as.character), stringsAsFactors=FALSE, row.names= rownames(pcl_df2)) #This needs to be done for the transcript data, but causes problems with the pcl data.
pcl_df2 <- rbind(pcl_metadata,pcl_df2)
pcl_df2 <- data.frame(t(pcl_df2))


save(pcl_metadata, pcl_df2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/data/Ternary/",fiber_subset,"_pcl_df.RData",sep=""))


library(lme4)
df = c(NA,NA,NA,NA)
for (i in (pcl_metadata_rows+1):ncol(pcl_df2)) { 
 # print(shapiro.test(as.numeric(pcl_df2[,i])))
 to_bind <- tryCatch({ 
    model <- lmer(as.numeric(pcl_df2[,i]) ~ 1 + (1|Dose) + (1|Participant) + Cholesterol.HDL.Ratio, data=pcl_df2, REML = FALSE)
  variance <- as.data.frame(VarCorr(model))
  print(variance)
  print(variance$vcov)
  to_bind <- variance$vcov/sum(variance$vcov)*100
 },
 error = function(err){print(err)
   return(c(NA,NA,NA))})
 to_bind <- c(to_bind,colnames(pcl_df2[i]))
 df <- rbind(df,to_bind)
} 
#
write.table(df, file = paste("/Users/SLancaster/Desktop/pclVariance",fiber_subset,"WeekPartCHDLR.txt",sep=""), sep="\t")
#df <- read.table("/Users/SLancaster/Desktop/varianceAxHomaPartWeek.txt", row.names=NULL)
#df <- df[,-1] #Only to be used when reading a df from a table
#df <- as.matrix(df) #Only to be used when reading a df from a table
df <- na.omit(df)
rownames(df) <- df[,4]
df <- df[,-4]
class(df) <- "numeric"
colnames(df) <- c("Participant", "Week", "Other")
#install.packages("ggtern")
df <- data.frame(df)
library(ggtern)
pdf(paste("/Users/SLancaster/Desktop/pclTernary",fiber_subset,"WeekPartCHDLR.pdf",sep=""), width=14)
ggtern(data = df,aes(x = Week ,y = Participant, z = Other)) +
  #geom_density_tern(aes(color=..level..), bins=4500) +
  #geom_density_tern(aes(fill=..level.., alpha=abs(..level..)),bins=500, binwidth=100) +
  #geom_point(size=0.5, alpha=0.5) +
  stat_density_tern(geom="polygon", n=4000, bins=500, aes(fill=..level..)) +
  scale_fill_gradient(low="yellow",high="red")   + 
  scale_color_gradient(low="yellow",high="red")   +
  theme_bw()  +
  theme(legend.justification=c(0,1), legend.position=c(0,1), panel.background = element_rect(fill = "light gray", colour = "light gray")) + 
  theme_nogrid() +
  guides(fill = guide_colorbar(order=1),
         alpha= guide_legend(order=2),
         color="none") + 
  labs(  title= "pcl",
         fill = "Value, V",alpha="|V - 0|")
dev.off()


