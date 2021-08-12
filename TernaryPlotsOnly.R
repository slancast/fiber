
##Only plotting ternary plot
#the following need to be pasted into the command line
#if (FALSE) {
#for i in *.txt
#do
#./TernaryPlotsOnly.R $i
#done
#}

files <- list.files(pattern = "RNA.*.txt")

for (i in files) { 
print(i)
j = substr(i, 4, nchar(i)-4)
df <- read.table(i, row.names=NULL)
df <- df[,-1] #Only to be used when reading a df from a table
df <- as.matrix(df) #Only to be used when reading a df from a table
df <- na.omit(df)
rownames(df) <- df[,4]
df <- df[,-4]
class(df) <- "numeric"
print(nrow(df))

if  ( (nrow(df) == 0) | (ncol(df) <= 2) ) {
  next
} 

colnames(df) <- c("Participant", "Week", "Other")
df <- data.frame(df)
print(df[1:50,])
library(ggtern)

print("plotting")
plot1 <- ggtern(data = df,aes(x = Week ,y = Participant, z = Other)) +
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

ggsave(plot1, file=paste("./",j,".pdf",sep=""))

}

