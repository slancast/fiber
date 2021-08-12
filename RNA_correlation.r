library(Mfuzz)
#The data should be in the standard data format. 
rna_datamatrix = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/DeseqNormRNAwithMetadata.txt", sep="\t", header=TRUE, row.names = 1, stringsAsFactors=FALSE)

rna_metadata_rows = 64
rna_datamatrix <- data.frame(t(rna_datamatrix))
rna_datamatrix <- data.frame(lapply(rna_datamatrix, as.character), stringsAsFactors=FALSE) #This needs to be done for the transcript data, but causes problems with the pcl data.
fiber_subset = "Arabinoxylan"
rna_datamatrix <- rna_datamatrix[which(rna_datamatrix$fiber==fiber_subset),] #subsetting by fiber
rna_datamatrix <- data.frame(t(rna_datamatrix))
rna_metadata = head(rna_datamatrix,rna_metadata_rows)



rna_datamatrix2 = tail(rna_datamatrix, -rna_metadata_rows) #Now getting rid of the metadata
rna_datamatrix2 <- as.matrix(t(rna_datamatrix2))
class(rna_datamatrix2) <- "numeric"

rna_datamatrix3 <- as.data.frame(rna_metadata["week",])
rna_datamatrix3 <- as.data.frame(t(rna_datamatrix3))
rna_datamatrix3 <- cbind(rna_datamatrix3, rna_datamatrix2)
rna_datamatrix3 <- aggregate(rna_datamatrix3,list(Visit=rna_datamatrix3$week), mean)
rna_datamatrix3 <- rna_datamatrix3[ , -2]
order = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10","Washoutweek3","WashoutFinal")
rna_datamatrix3 <- rna_datamatrix3[match(order, rna_datamatrix3$Visit),]
rna_datamatrix3 <- na.omit(rna_datamatrix3)
rownames(rna_datamatrix3) <- rna_datamatrix3$Visit


rna_datamatrix5 <- rna_datamatrix3
rna_datamatrix5 <- data.frame(t(rna_datamatrix5[,-1]))#Needs to be a matrix to be an expressionset, and cytokines need to be as rows for standardization
rna_datamatrix5 <- as.matrix(rna_datamatrix5)
class(rna_datamatrix5) <- "numeric"
rna_datamatrix5 <- rna_datamatrix5[rowMeans(rna_datamatrix5)>0,] #If all the rows are 0, 
rna_eset1 <- ExpressionSet(na.omit(rna_datamatrix5)) #Creating the type expression set with the metadata rows as a different argument
rna_eset1 <- standardise(rna_eset1) #Running standarise 
rna_m <- exprs(rna_eset1)


library("Hmisc")
# Plot correlation graph
m <- rna_m
cor <- rcorr(t(m), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
FDR.data <- apply(pval.data,2,p.adjust,method="fdr", n = length(pval.data))
pdf("./pval_hist.pdf")
hist(pval.data, breaks = 100, col="darkblue")
dev.off()
pdf("./FDR_hist.pdf")
hist(FDR.data, breaks = 100, col="darkblue")
dev.off()
pdf("./cor_hist.pdf")
hist(cor.data, breaks = 10, col="red")
dev.off()
cor.data[FDR.data > 0.05]=0
#cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
library("igraph")
network=graph.adjacency(cor.data, weighted=T, mode="undirected", diag=F)
network <- delete.vertices(network,which(degree(network)<1))
print(network)
cl <- clusters(network)
gs <- induced.subgraph(network, cl$membership %in% which(cl$csize>40))

attr <- cbind(id=1:vcount(gs), val=(vcount(gs)+1))
g <- gs + vertices(unique(attr[,2])) 
g <- g + igraph::edges(unlist(t(attr)), weight=0.1)
ly <- layout.fruchterman.reingold(g,dim=2,grid="nogrid",weights=E(g)$weight)
ly <- ly[1:(nrow(ly)-1),]
# plot
pdf("./network_large_networks.pdf")
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(g,
     vertex.size=0.3,
     #vertex.color=V(gs)$color,
     vertex.label.cex=0.5,
     vertex.label.color="black",
     vertex.frame.color="black",
     vertex.label = NA,
     layout = ly
)
dev.off()


lab <- data.frame(c(rep("Metabolites", length(Assay[which(Assay == "Metabolomics")])),
                    rep("Cytokines", length(Assay[which(Assay == "Cytokine")])),
                    rep("Clinical", length(Assay[which(Assay == "Clinical")])),
                    rep("Proteins", length(Assay[which(Assay == "Proteomics")]))),
                  rownames(data))
colnames(lab) <- c("Category","Analyte")
V(gs)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Metabolites")] <- "blue"
V(gs)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Cytokines")] <- "orange"
V(gs)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Clinical")] <- "green"
V(gs)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Proteins")] <- "purple"
