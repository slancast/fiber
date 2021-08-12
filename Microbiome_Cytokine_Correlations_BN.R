load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Mix_cytokine_df.RData")
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Mix_metaphlan_df.RData")

#Generating row names where the the formatting is the same between the two datasets
#This involves only getting rid of some characters from from the pcl row names
a <- rownames(normalized_metaphlan_df)
tmp <- gsub("_FM.*kneaddata_metaphlan_bugs_list","",a)
tmp <- as.character(tmp)
rownames(normalized_metaphlan_df) <- tmp
a <- colnames(normalized_metaphlan_df)
tmp <- gsub("k__.*g__","",a)
tmp <- as.character(tmp)
colnames(normalized_metaphlan_df) <- tmp

a <- rownames(normalized_metaphlan_df)
tmp <- gsub("LC_Inulin","LCInulin",a)
tmp <- as.character(tmp)
rownames(normalized_metaphlan_df) <- tmp

a <- rownames(normalized_cytokine_df)
tmp <- gsub("LC_Inulin","LCInulin",a)
tmp <- as.character(tmp)
rownames(normalized_cytokine_df) <- tmp

##################################################################

#This will find anyting that matches the strings in x and throw out those entries in y
#A way to get rid of all the long term values 
#grepl is supposed to be quicker because it uses boolean functions
y <- rownames(normalized_cytokine_df)
x <- c("Final", "Week", "Month")
tmp_cyto <- y[!grepl(paste0(x, collapse = "|"), y)]
##################################################################

# Finding the samples that match
matching_samples <- !is.na( match(rownames(normalized_cytokine_df), rownames(normalized_metaphlan_df)))
cytokine_matching_df <- normalized_cytokine_df[matching_samples,]

normalized_metaphlan_df <- cbind(names=rownames(normalized_metaphlan_df),normalized_metaphlan_df)
matching_samples2 <- !is.na( match(rownames(normalized_metaphlan_df),rownames(cytokine_matching_df)))
pcl_matching_df <- normalized_metaphlan_df[matching_samples2,]
##################################################################

# The pcl dataframe has extra information in it
# First getting rid of the pathways
# The reducing the dubpicated rows by their average
pcl_matching_df <- pcl_matching_df[,1:ncol(pcl_matching_df)] #Keeping only the genera, and not pathways
names <- pcl_matching_df[,1] #Just for formatting concerns I had to split the matrix into the part that was character and the part that was numeric
pcl_matching_df <- pcl_matching_df[,2:ncol(pcl_matching_df)]
class(pcl_matching_df) <- "numeric" #For some reason the code requires this formatting
pcl_matching_df <- cbind(data.frame(names),data.frame(pcl_matching_df)) #Creating a data frame with the correct poritons as numeric or character

deduped.pcl <- aggregate(pcl_matching_df,list(names=pcl_matching_df$names),FUN=mean) #Aggregating to get rid of redundant rownames
rownames(deduped.pcl) <- deduped.pcl$names 
deduped.pcl <- deduped.pcl[,-c(1,2)] #Getting rid of extra columns
##################################################################

match(rownames(cytokine_matching_df),rownames(deduped.pcl)) #Double checking all the rows match
combined_df <- cbind(cytokine_matching_df,deduped.pcl) #finally creating the combined df

# Now I will run the rcorr on these data
# I believe the correct thing to do is run standardise before correlation
# If the difference in magnitude of the two cytokines is large than the line might be close to 0/1 with little correlation
# even though they change together

# This simple function will run standardise
# Takes in a data.frame object and will return a matrix
standardise_data <- function(data_frame) {
require(Mfuzz)
n <- as.matrix(data_frame)
class(n) <- "numeric"
eset1 <- ExpressionSet(n)
eset1 <- standardise(eset1) #Running standarise 
o <- exprs(eset1)
o <- na.omit(o)
return(o)
}
combined_df <- data.frame(t(combined_df))
a <- rownames(combined_df)
tmp <- gsub("k__.*g__","",a)
tmp <- as.character(tmp)
rownames(combined_df) <- tmp
o <- standardise_data(combined_df)
##################################################################

#Now that standardise has run I want to correlate the various data with one another
#This function will take in a matrix 
correlations <- function(to_correlate_df) {
require(Hmisc)
cor <- rcorr(format(t(to_correlate_df),digits=20), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
FDR.data <- apply(pval.data,2,p.adjust,method="bonferroni", n = length(pval.data))
pdf("./pval_bonferonni_hist.pdf")
hist(pval.data, breaks = 100, col="darkblue")
dev.off()
pdf("./FDR_bonferonni_hist.pdf")
hist(FDR.data, breaks = 100, col="darkblue")
dev.off()
pdf("./cor_bonferonni_hist.pdf")
hist(cor.data, breaks = 10, col="red")
dev.off()
cor.data[FDR.data > 0.75]=0
return(cor.data)
}

cor.data <- correlations(o)
##################################################################


library("igraph")
network=graph.adjacency(cor.data, weighted=T, mode="undirected", diag=F)
#V(network)$color <- lab[,2]
print(network)
ly <- layout.fruchterman.reingold(network,dim=2,grid="nogrid")

pdf("/Users/SLancaster/Desktop/Correlations.pdf")
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network,
     vertex.size=10,
     #vertex.color=V(subgraph)$color,
     vertex.label.cex=0.5,
     vertex.label.color="black",
     vertex.frame.color="black",
     #vertex.label = NA,
     layout = ly
)
dev.off()

