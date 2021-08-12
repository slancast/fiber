
#The data should be in the standard data format. 
load("/home/slancast/NormAggreg_Arabinoxylan_aggregnorm_cytokine.RData")
load("/home/slancast/NormAggreg_Arabinoxylan_metabolomics_df.RData")
load("/home/slancast/NormAggreg_Arabinoxylan_pcl.RData")
load("/home/slancast/NormAggreg_Arabinoxylan_rna.RData")

aggregnorm_rna <- data.frame(aggregnorm_rna)
aggregnorm_rna <- cbind(aggregnorm_rna, Gene = rownames(aggregnorm_rna))

aggregnorm_rna <- entrezid(aggregnorm_rna)
rownames(aggregnorm_rna) <- make.names(aggregnorm_rna$SYMBOL, unique=TRUE)
aggregnorm_rna <- subset(aggregnorm_rna, select=-c(SYMBOL,ENTREZID,Gene))
aggregnorm_rna <- as.matrix(aggregnorm_rna)
colnames(aggregnorm_rna) <- colnames(aggregnorm_pcl)

m <- rbind(aggregnorm_pcl, aggregnorm_rna, aggregnorm_metabolomics, aggregnorm_cytokine) 
lab <- c(rep("Microbiome",nrow(aggregnorm_pcl)),rep("RNA",nrow(aggregnorm_rna)),rep("Metabolites",nrow(aggregnorm_metabolomics)),rep("Cytokine",nrow(aggregnorm_cytokine)))
color <- c(rep("blue",nrow(aggregnorm_pcl)),rep("red",nrow(aggregnorm_rna)),rep("yellow",nrow(aggregnorm_metabolomics)),rep("orange",nrow(aggregnorm_cytokine)))
lab <- cbind(lab,color)
lab <- data.frame(lab)
rownames(lab) <- rownames(m)

library(Mfuzz)
m <- as.matrix(m)
class(m) <- "numeric"
m <- m[rowMeans(m)>0,] #If the rows sum to more than 0
rownames(m) <- make.names(rownames(m), unique=TRUE) 
rna_eset1 <- ExpressionSet(na.omit(m)) #Creating the type expression set with the metadata rows as a different argument
rna_eset1 <- standardise(rna_eset1) #Running standarise 
m <- exprs(rna_eset1)

library("Hmisc")
# Plot correlation graph
cor <- rcorr(format(t(m),digits=20), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
#FDR.data <- apply(pval.data,2,p.adjust,method="BY", n = length(pval.data))
pdf("./pval_bonferonni_hist.pdf")
hist(pval.data, breaks = 100, col="darkblue")
dev.off()
pdf("./FDR_bonferonni_hist.pdf")
hist(FDR.data, breaks = 100, col="darkblue")
dev.off()
pdf("./cor_bonferonni_hist.pdf")
hist(cor.data, breaks = 10, col="red")
dev.off()
cor.data[pval.data > 0.00000001]=0
save(cor.data, file="bonderroni_corrected_cor.data.RData")
#load("spear_bonferroni__corrected_cor.data.RData")
#cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
library("igraph")
network=graph.adjacency(cor.data, weighted=T, mode="undirected", diag=F)
V(network)$color <- lab[,2]
print(network)

if (FALSE) {
for (i in rownames(aggregnorm_cytokine)) {
#finding subcluster by vertex name
print(i)
vertex <- match(i, V(network)$name)
cluster_number <- clusters(network)$membership[vertex]}}
dg <- decompose.graph(network)
for (j in 1:length(dg)) {
	subgraph <- dg[[j]]
	counter <- 0 
for (k in 0:(length(V(subgraph))-1)) {
	counter = counter + k
	}
if (length(E(subgraph)) == counter) {
print(E(subgraph))
print(length(V(subgraph))-1)
print(counter)
ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")

pdf(paste("./p-8/",j,"_subgraph.pdf",sep=""))
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(subgraph,
     vertex.size=10,
     vertex.color=V(subgraph)$color,
     vertex.label.cex=0.5,
     vertex.label.color="black",
     vertex.frame.color="black",
     #vertex.label = NA,
     layout = ly
)
dev.off()
}}
