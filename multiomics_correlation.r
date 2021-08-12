library(Mfuzz)
#The data should be in the standard data format. 
load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/data/Clustering/Arabinoxylan_pcl_df.RData")

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/data/Clustering/Arabinoxylan_rna_df.RData")

m <- rbind(pcl_df5, rna_df5) 



library("Hmisc")
# Plot correlation graph
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
