fiber_subset = "Arabinoxylan"
load(paste("/home/slancast/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""))
log_pcl <- cbind(t(pcl_metadata),log_pcl)

for (j in rownames(pcl_metadata)) {
  print(j)
df = c(NA,NA,NA,NA)
library(lme4)
for (i in nrow(pcl_metadata):ncol(log_pcl)) { 
 #print(shapiro.test(as.numeric(log_pcl[,i])))
  
 to_bind <- tryCatch({ 
    model <- glmer(as.numeric(log_pcl[,i]) ~ 1 + (1|Dose) + (1|Participant) + as.numeric(log_pcl[,i]), data=log_pcl)
  variance <- as.data.frame(VarCorr(model))
  to_bind <- variance$vcov/sum(variance$vcov)*100
 },
 error = function(err){print(err)
   return(c(NA,NA,NA))})
 to_bind <- c(to_bind,colnames(log_pcl[i]))
 df <- rbind(df,to_bind)
} 
#
write.table(df, file = paste("~/pclTernary/pcl",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
#df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/pclVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
#df <- df[,-1] #Only to be used when reading a df from a table
#df <- as.matrix(df) #Only to be used when reading a df from a table
df <- na.omit(df)
rownames(df) <- df[,4]
df <- df[,-4]
class(df) <- "numeric"
if  ( (nrow(df) == 0) | (ncol(df) <= 2) ) {
  next
} 
colnames(df) <- c("Participant", "Week", "Other")
#install.packages("ggtern")
df <- data.frame(df)
library(ggtern)
plot1 <- ggtern(data = df,aes(x = Week ,y = Participant, z = Other)) +
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
ggsave(plot1, file=paste("~/pclTernary/",fiber_subset,j,"WeekPartHematoTri.pdf",sep=""))
}
