
# Plot correlation graph
cor <- rcorr(t(data), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
FDR.data <- apply(pval.data,2,p.adjust,method="bonferroni", n = length(pval.data))
hist(pval.data, breaks = 100, col="darkblue")
hist(FDR.data, breaks = 100, col="darkblue")
hist(cor.data, breaks = 10, col="red")
cor.data[FDR.data > 0.05]=0
network=graph.adjacency(cor.data, weighted=T, mode="undirected", diag=F)
network <- delete.vertices(network,which(degree(network)<1))
network
lab <- data.frame(c(rep("Metabolites", length(Assay[which(Assay == "Metabolomics")])), 
                    rep("Cytokines", length(Assay[which(Assay == "Cytokine")])),
                    rep("Clinical", length(Assay[which(Assay == "Clinical")])),
                    rep("Proteins", length(Assay[which(Assay == "Proteomics")]))),
                  rownames(data))
colnames(lab) <- c("Category","Analyte")
V(network)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Metabolites")] <- "blue"
V(network)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Cytokines")] <- "orange"
V(network)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Clinical")] <- "green"
V(network)$color[which(lab$Category[lab$Analyte %in% V(network)$name] == "Proteins")] <- "purple"

# plot
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network, 
     vertex.size=8,
     vertex.color=V(network)$color, 
     vertex.label.cex=0.5,
     vertex.label.color="black",
     vertex.frame.color="black",
     #vertex.label = NA,
     layout = layout.fruchterman.reingold(network,dim=2)
)