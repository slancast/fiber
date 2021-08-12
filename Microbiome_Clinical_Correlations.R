# Mike is very interested in the clinical values so I will correlate them with the microbiome
# 
# 

for (fiber_subset in c("Arabinoxylan", "LCInulin", "Mix")) {
  
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_",fiber_subset,"_clinical_df.RData",sep=""))
load(paste("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Full_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))

#Generating row names where the the formatting is the same between the two datasets
#This involves only getting rid of some characters from from the pcl row names
a <- rownames(log_metaphlan) #For the pcl files
tmp <- gsub("_FM.*kneaddata_Abundance","",a)
tmp <- as.character(tmp)
rownames(log_metaphlan) <- tmp
a <- rownames(log_metaphlan) #For the bugs list
tmp <- gsub("_FM.*kneaddata_metaphlan_bugs_list","",a)
tmp <- as.character(tmp)
rownames(log_metaphlan) <- make.names(tmp, unique = TRUE)
a <- rownames(log_metaphlan) #For the bugs list
tmp <- gsub("X","",a)
tmp <- as.character(tmp)
rownames(log_metaphlan) <- tmp
a <- colnames(log_metaphlan)
tmp <- gsub("k__.*g__","",a)
tmp <- as.character(tmp)
colnames(log_metaphlan) <- tmp
a <- rownames(log_metaphlan)
tmp <- gsub("LC_Inulin","LCInulin",a)
tmp <- as.character(tmp)
rownames(log_metaphlan) <- tmp

#Clinicals
rownames(clinical_df2) <- make.names(as.character(unlist(clinical_metadata[1,])), unique=TRUE)
a <- rownames(clinical_df2) #For the bugs list
tmp <- gsub("X","",a)
tmp <- as.character(tmp)
rownames(clinical_df2) <- tmp
a <- rownames(clinical_df2)
tmp <- gsub("WashoutD","Washout_D",a)
tmp <- as.character(tmp)
rownames(clinical_df2) <- tmp
to_toss <- which(colnames(clinical_df2) %in% c("date_time_blood_draw_clinicals"))
clinical_df2 <- clinical_df2[,-to_toss]

##################################################################

# Finding the samples that match
matching_samples <- !is.na( match(rownames(clinical_df2), rownames(log_metaphlan)))
clinical_matching_df <- clinical_df2[matching_samples,]

names <- rownames(log_metaphlan)
matching_samples2 <- !is.na( match(rownames(log_metaphlan),rownames(clinical_matching_df)))
pcl_matching_df <- log_metaphlan[matching_samples2,]

##################################################################

match(rownames(clinical_matching_df),rownames(pcl_matching_df)) #Double checking all the rows match
combined_df <- cbind(clinical_matching_df,pcl_matching_df) #finally creating the combined df

# Now I will run the rcorr on these data
# I believe the correct thing to do is run standardise before correlation
# If the difference in magnitude of the two cytokines is large than the line might be close to 0/1 with little correlation
# even though they change together

# This simple function will run standardise
# Takes in a data.frame object and will return a matrix
source("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Software/Utils.r")
combined_df <- as.matrix(t(combined_df))
#The clinical data is rife with missing values so I will try to impute them here
#Need to do this before standardise because often NAs are generated there
library(impute)
combined_df.i <- impute.knn(combined_df)
o <- standardise_data(combined_df.i$data)
##################################################################

#Now that standardise has run I want to correlate the various data with one another
#This function will take in a matrix 
cor.data <- correlations(o)
##################################################################


library("igraph")
#Need to make abs when I do this because otherwise it can ignore the negative.
network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)
#V(network)$color <- lab[,2]
print(network)
ly <- layout.fruchterman.reingold(network,dim=2,grid="nogrid")

#The following lines will find if  the correlations are negative or positive
#and then add the properly colored edge 
correlation_values <- c()
edges <- data.frame(get.edgelist(network))
for (i in 1:nrow(edges)){
  j <- as.matrix(edges[i,])
  correlation_value <- cor.data[as.character(j[2]),as.character(j[1])]
  correlation_values <- c(correlation_values,correlation_value)
}
E(network)$color <- ifelse(correlation_values < 0, "red","blue")

colors <- c(rep("yellow", length=ncol(clinical_matching_df)), rep("lightblue",length=(ncol(cor.data)-ncol(clinical_df2))))
V(network)$color <- colors

pdf(paste("/Users/SLancaster/Desktop/",fiber_subset,"Correlations.pdf",sep=""))
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network,
     vertex.size=5,
     #vertex.color=V(subgraph)$color,
     vertex.label.cex=0.5,
     vertex.label.color="black",
     vertex.frame.color="black",
     #vertex.label = NA,
     layout = ly,
     edge.width=E(network)$weight
)
dev.off()

}

