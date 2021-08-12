fiber_subsets = c("Arabinoxylan","LCInulin","Mix")
time_points = c("Baseline", "10", "20", "30", "WashoutD3", "WashoutD10")

for (fiber_subset in fiber_subsets) {
  
load(paste("~/Normalized_Log_",fiber_subset,"_cytokine_df.RData",sep=""))
load(paste("~/Normalized_Log_",fiber_subset,"_pcl_df.RData",sep=""))
load(paste("~/Normalized_Log_",fiber_subset,"_rna_df.RData",sep=""))
load(paste("~/Normalized_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))
load(paste("~/Normalized_Log_",fiber_subset,"_lipids_df.RData",sep=""))
load(paste("~/Normalized_Log_",fiber_subset,"_proteomics_df.RData",sep=""))

load(paste("~/Full_Log_",fiber_subset,"_cytokine_df.RData",sep=""))
load(paste("~/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""))
load(paste("~/Full_Log_",fiber_subset,"_rna_df.RData",sep=""))
load(paste("~/Full_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))
load(paste("~/Full_Log_",fiber_subset,"_lipids_df.RData",sep=""))
load(paste("~/Full_Log_",fiber_subset,"_proteomics_df.RData",sep=""))

normalized_rna_df <- data.frame(t(normalized_rna_df))
normalized_rna_df <- cbind(normalized_rna_df, Gene = rownames(normalized_rna_df))
entrezid <- function( resdata ) {
  require(EnsDb.Hsapiens.v79)
  require(AnnotationDbi)  #the column to iterate over will be different if I'm using res vs resdata
  a = resdata$Gene #the column to iterate over will be different if I'm using res vs resdat
  tmp=gsub("\\..*","",a)
  tmp <- as.character(tmp)
  txdb <- EnsDb.Hsapiens.v79
  df <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "ENTREZID")
  df2 <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "SYMBOL")
  
  ENTREZID <- c()
  SYMBOL <- c()
  counter1 <- 0
  for (i in tmp) {
    counter1 <- counter1 + 1
    j <- match(i,df$GENEID)
    ENTREZID <- c(ENTREZID, toString(df[j,][2]))
    SYMBOL <- c(SYMBOL, toString(df2[j,][2]))}
  resdata$ENTREZID <- ENTREZID
  resdata$SYMBOL <- SYMBOL
  resdata$Gene_woversion <- tmp
  
  resdata
}

normalized_rna_df <- entrezid(normalized_rna_df)
rownames(normalized_rna_df) <- make.names(normalized_rna_df$Gene_woversion, unique=TRUE)
normalized_rna_df <- subset(normalized_rna_df, select=-c(SYMBOL,ENTREZID,Gene,Gene_woversion))
normalized_rna_df <- as.matrix(normalized_rna_df)
colnames(normalized_rna_df) <- colnames(normalized_rna_df)

#For this section I am concatenating the the time value to the log normalized values
rna_week <- data.frame(rna_metadata["week",]) #The syntax was easier to do it in two steps like this
pcl_week <- data.frame(pcl_metadata["Dose",])
cytokine_week <- data.frame(cytokine_metadata["Week",])
metabolomics_week <- data.frame(metabolomics_metadata["timepoint",])
lipids_week <- data.frame(lipids_metadata["timepoint",])
proteomics_week <- data.frame(proteomics_metadata["week",])
normalized_rna_df <- cbind(Dose = t(rna_week),t(normalized_rna_df))
normalized_pcl_df <- cbind(Dose = t(pcl_week),normalized_pcl_df)
normalized_cytokine_df <- cbind(Dose = t(cytokine_week),normalized_cytokine_df)
normalized_metabolomics_df <- cbind(Dose = t(metabolomics_week),normalized_metabolomics_df)
normalized_lipids_df <- cbind(Dose = t(lipids_week),normalized_lipids_df)
normalized_proteomics_df <- cbind(Dose = t(proteomics_week),normalized_proteomics_df)
normalized_rna_df <- data.frame(normalized_rna_df)
normalized_pcl_df <- data.frame(normalized_pcl_df)
normalized_cytokine_df <- data.frame(normalized_cytokine_df)
normalized_metabolomics_df <- data.frame(normalized_metabolomics_df)
normalized_lipids_df <- data.frame(normalized_lipids_df)
normalized_proteomics_df <- data.frame(normalized_proteomics_df)

for (time_point in time_points) {

rna_df  <- normalized_rna_df[which(normalized_rna_df$week==time_point),]
pcl_df  <- normalized_pcl_df[which(normalized_pcl_df$Dose==time_point),]
cytokine_df  <- normalized_cytokine_df[which(normalized_cytokine_df$Week==time_point),]
metabolomics_df  <- normalized_metabolomics_df[which(normalized_metabolomics_df$timepoint==time_point),]
lipids_df  <- normalized_lipids_df[which(normalized_lipids_df$timepoint==time_point),]
proteomics_df <- normalized_proteomics_df[which(normalized_proteomics_df$timepoint==time_point),]

##########################################################################
############################## lipids ####################################
##########################################################################
#if (FALSE){
library(network)
library(Mfuzz)
library(matrixStats)
n <- as.matrix(t(lipids_df))
class(n) <- "numeric"
eset1 <- ExpressionSet(n)
eset1 <- standardise(eset1) #Running standarise
o <- exprs(eset1)
o <- na.omit(o)
lipid_class_tmp <- lipid_class[which(names(lipid_class) %in% rownames(o))]

fdr_method = "BH"
library(icesTAF)
mkdir(paste("./lognormlipids",fiber_subset,time_point,fdr_method, sep=""))

library("Hmisc")
# Plot correlation graph
print(time_point)
cor <- rcorr(format(t(o),digits=20), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
FDR.data <- apply(pval.data,2,p.adjust,method=fdr_method, n = length(pval.data))
pdf(paste("./lognormlipids",fiber_subset,time_point,fdr_method,"/pval_hist.pdf",sep=""))
hist(pval.data, breaks = 100, col="darkblue")
dev.off()
pdf(paste("./lognormlipids",fiber_subset,time_point,fdr_method,"/FDR_hist.pdf",sep=""))
hist(FDR.data, breaks = 100, col="darkblue")
dev.off()
pdf(paste("./lognormlipids",fiber_subset,time_point,fdr_method,"/cor_hist.pdf",sep=""))
hist(cor.data, breaks = 10, col="red")
dev.off()
cor.data[FDR.data > 0.05]=0
write.table(cor.data,file=paste("./lognormlipids",fiber_subset,time_point,fdr_method,"/cor_data.txt",sep=""), sep="\t")
#load("spear_bonferroni__corrected_cor.data.RData")
#cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
library("igraph")
network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)


#The following lines will find if  the correlations are negative or positive
#and then add the properly colored edge 
correlation_values <- c()
edges <- data.frame(get.edgelist(network))
for (i in 1:nrow(edges)){
  j <- as.matrix(edges[i,])
  correlation_value <- cor.data[as.character(j[2]),as.character(j[1])]
  correlation_values <- c(correlation_values,correlation_value)
}
E(network)$color <- ifelse(correlation_values < 0, "red","green")


colors = rep("red", length=length(lipid_class_tmp))
colors[which(lipid_class_tmp=="CE")] = "red"
colors[which(lipid_class_tmp=="CER")] = "gray"
colors[which(lipid_class_tmp=="DAG")] = "blue"
colors[which(lipid_class_tmp=="FFA")] = "yellow"
colors[which(lipid_class_tmp=="HCER")] = "darkorange"
colors[which(lipid_class_tmp=="LCER")] = "chartreuse"
colors[which(lipid_class_tmp=="LPC")] = "cyan"
colors[which(lipid_class_tmp=="LPE")] = "brown"
colors[which(lipid_class_tmp=="PC")] = "purple"
colors[which(lipid_class_tmp=="PE")] = "limegreen"
colors[which(lipid_class_tmp=="SM")] = "aquamarine"
colors[which(lipid_class_tmp=="TAG")] = "deeppink1"
colors <- as.character(colors)

V(network)$color <- colors
print(network)

ly <- layout.fruchterman.reingold(network,dim=2,grid="nogrid")

pdf(paste("./lognormlipids",fiber_subset,time_point,fdr_method,"/Full_network.pdf",sep=""))
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network,
     vertex.size=3,
     vertex.color=V(network)$color,
     #vertex.label.cex=0.5,
     #vertex.label.color="black",
     #vertex.frame.color="black",
     vertex.label = NA,
     layout = ly
)
dev.off()

##########################################################################
#Printing All subgraphs###################################################
##########################################################################

dg <- decompose.graph(network)

for (j in 1:length(dg)) {
  subgraph <- dg[[j]]
  if (length(V(subgraph)) == 1) { next }
  else {
  counter <- 0
  ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")

  pdf(paste("./lognormlipids",fiber_subset,time_point,fdr_method,"/",j,fiber_subset,time_point,"_subgraph.pdf",sep=""))
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
} }

##########################################################################
########################### END ##########################################
##########################################################################
#}

##########################################################################
########################### Metabolomics #################################
##########################################################################
#Same length as metabolomics data
# 
library(network)
library(Mfuzz)
library(matrixStats)
n <- as.matrix(t(metabolomics_df))
class(n) <- "numeric"
eset1 <- ExpressionSet(n)
eset1 <- standardise(eset1) #Running standarise
o <- exprs(eset1)
o <- na.omit(o)
metabolome_superpathway_temp <- metabolome_superpathway[,which(colnames(metabolome_superpathway) %in% rownames(o))]


fdr_method = "BH"
library(icesTAF)
mkdir(paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method, sep=""))

library("Hmisc")
# Plot correlation graph
print(time_point)
cor <- rcorr(format(t(o),digits=20), type="spearman")
cor.data <- cor$r
cor.data[upper.tri(cor.data, diag = T)]<- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)]<- NA
FDR.data <- apply(pval.data,2,p.adjust,method=fdr_method, n = length(pval.data))
pdf(paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method,"/pval_hist.pdf",sep=""))
hist(pval.data, breaks = 100, col="darkblue")
dev.off()
pdf(paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method,"/FDR_hist.pdf",sep=""))
hist(FDR.data, breaks = 100, col="darkblue")
dev.off()
pdf(paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method,"/cor_hist.pdf",sep=""))
hist(cor.data, breaks = 10, col="red")
dev.off()
cor.data[FDR.data > 0.05]=0
write.table(cor.data,file=paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method,"/cor_data.txt",sep=""), sep="\t")
#load("spear_bonferroni__corrected_cor.data.RData")
#cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
library("igraph")
network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)
colors = rep("red", length=length(metabolome_superpathway_temp))
colors[which(metabolome_superpathway_temp[1,]=="Amino Acid")] = "red"
colors[which(metabolome_superpathway_temp[1,]=="Peptide")] = "red"
colors[which(metabolome_superpathway_temp[1,]=="Other")] = "blue"
colors[which(metabolome_superpathway_temp[1,]=="Lipid")] = "yellow"
colors[which(metabolome_superpathway_temp[1,]=="Cofactors and Vitamins")] = "darkorange"
colors[which(metabolome_superpathway_temp[1,]=="Nucleotide")] = "chartreuse"
colors[which(metabolome_superpathway_temp[1,]=="Xenobiotics")] = "cyan"
colors[which(metabolome_superpathway_temp[1,]=="Carbohydrate")] = "brown"
colors <- as.character(colors)

V(network)$color <- colors
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
E(network)$color <- ifelse(correlation_values < 0, "red","green")


pdf(paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method,"/Full_network.pdf",sep=""))
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network,
     vertex.size=3,
     vertex.color=V(network)$color,
     #vertex.label.cex=0.5,
     #vertex.label.color="black",
     #vertex.frame.color="black",
     vertex.label = NA,
     layout = ly
)
dev.off()

##########################################################################
#Printing All subgraphs###################################################
##########################################################################

dg <- decompose.graph(network)
#
for (j in 1:length(dg)) {
  subgraph <- dg[[j]]
  counter <- 0
  ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
#
  pdf(paste("./lognormmetabolomics",fiber_subset,time_point,fdr_method,"/",j,fiber_subset,time_point,"_subgraph.pdf",sep=""))
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
}

##########################################################################
########################### END ##########################################
##########################################################################


##########################################################################
########################### RNA ##########################################
##########################################################################
# 
# library(Mfuzz)
# library(matrixStats)
# n <- as.matrix(t(rna_df))
# class(n) <- "numeric"
# eset1 <- ExpressionSet(n)
# eset1 <- standardise(eset1) #Running standarise 
# o <- exprs(eset1)
# o <- na.omit(o)
# 
# library("Hmisc")
# # Plot correlation graph
# cor <- rcorr(format(t(o),digits=20), type="spearman")
# cor.data <- cor$r
# cor.data[upper.tri(cor.data, diag = T)]<- 0
# pval.data <- cor$P
# pval.data[upper.tri(pval.data, diag = T)]<- NA
# FDR.data <- apply(pval.data,2,p.adjust,method="bonferroni", n = length(pval.data))
# pdf("./pval_bonferonni_hist.pdf")
# hist(pval.data, breaks = 100, col="darkblue")
# dev.off()
# pdf("./FDR_bonferonni_hist.pdf")
# hist(FDR.data, breaks = 100, col="darkblue")
# dev.off()
# pdf("./cor_bonferonni_hist.pdf")
# hist(cor.data, breaks = 10, col="red")
# dev.off()
# cor.data[FDR.data > 0.05]=0
# save(cor.data, file="bonderroni_corrected_cor.data.RData")
# #load("spear_bonferroni__corrected_cor.data.RData")
# #cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
# library("igraph")
# network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)
# #V(network)$color <- lab[,2]
# print(network)
# 
# library(icesTAF)
# mkdir(paste("./lognormrna",fiber_subset,time_point,sep=""))
# 
# ##########################################################################
# #Printing All subgraphs###################################################
# ##########################################################################
# 
# dg <- decompose.graph(network)
# total_genes <- as.list(V(network))
# write.table(total_genes, file = paste("./lognormrna",fiber_subset,time_point,"/",fiber_subset,time_point,"_total_genes.txt",sep=""),sep="\t")
# 
# for (j in 1:length(dg)) {
#   subgraph <- dg[[j]]
#   vertex_list <- as.list(V(subgraph))
#   write.table(vertex_list, file = paste("./lognormrna",fiber_subset,time_point,"/",j,fiber_subset,time_point,"_subgraph.txt",sep=""),sep="\t")
# }
# 
# for (j in 1:length(dg)) {
#   subgraph <- dg[[j]]
#   counter <- 0
#   ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
#   
#   pdf(paste("./lognormrna",fiber_subset,time_point,"/",j,fiber_subset,time_point,"_subgraph.pdf",sep=""))
#   par(bg="white", mar=c(0,0,0,0))
#   set.seed(4)
#   plot(subgraph,
#        vertex.size=10,
#        vertex.color=V(subgraph)$color,
#        vertex.label.cex=0.5,
#        vertex.label.color="black",
#        vertex.frame.color="black",
#        vertex.label = NA,
#        layout = ly
#   )
#   dev.off()
# }
##########################################################################
########################### END ##########################################
##########################################################################

##########################################################################
########################### PCL ##########################################
##########################################################################
# 
# library(Mfuzz)
# library(matrixStats)
# n <- as.matrix(t(pcl_df))
# class(n) <- "numeric"
# eset1 <- ExpressionSet(n)
# eset1 <- standardise(eset1) #Running standarise 
# o <- exprs(eset1)
# o <- na.omit(o)
# 
# library("Hmisc")
# # Plot correlation graph
# cor <- rcorr(format(t(o),digits=20), type="spearman")
# cor.data <- cor$r
# cor.data[upper.tri(cor.data, diag = T)]<- 0
# pval.data <- cor$P
# pval.data[upper.tri(pval.data, diag = T)]<- NA
# FDR.data <- apply(pval.data,2,p.adjust,method="bonferroni", n = length(pval.data))
# pdf("./pval_bonferonni_hist.pdf")
# hist(pval.data, breaks = 100, col="darkblue")
# dev.off()
# pdf("./FDR_bonferonni_hist.pdf")
# hist(FDR.data, breaks = 100, col="darkblue")
# dev.off()
# pdf("./cor_bonferonni_hist.pdf")
# hist(cor.data, breaks = 10, col="red")
# dev.off()
# cor.data[FDR.data > 0.05]=0
# save(cor.data, file="bonderroni_corrected_cor.data.RData")
# #load("spear_bonferroni__corrected_cor.data.RData")
# #cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
# library("igraph")
# network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)
# #V(network)$color <- lab[,2]
# print(network)
# 
# library(icesTAF)
# mkdir(paste("./lognormpcl",fiber_subset,time_point,sep=""))
# 
##########################################################################
#Printing All subgraphs###################################################
##########################################################################
# 
# dg <- decompose.graph(network)
# 
# for (j in 1:length(dg)) {
#   subgraph <- dg[[j]]
#   counter <- 0
#   ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
#   
#   pdf(paste("./lognormpcl",fiber_subset,time_point,"/",j,fiber_subset,time_point,"_subgraph.pdf",sep=""))
#   par(bg="white", mar=c(0,0,0,0))
#   set.seed(4)
#   plot(subgraph,
#        vertex.size=10,
#        vertex.color=V(subgraph)$color,
#        vertex.label.cex=0.5,
#        vertex.label.color="black",
#        vertex.frame.color="black",
#        #vertex.label = NA,
#        layout = ly
#   )
#   dev.off()
# }
##########################################################################
########################### END ##########################################
##########################################################################


##########################################################################
########################### Cytokine #####################################
##########################################################################
# 
# library(Mfuzz)
# library(matrixStats)
# n <- as.matrix(t(cytokine_df))
# class(n) <- "numeric"
# eset1 <- ExpressionSet(n)
# eset1 <- standardise(eset1) #Running standarise 
# o <- exprs(eset1)
# o <- na.omit(o)
# 
# library("Hmisc")
# # Plot correlation graph
# cor <- rcorr(format(t(o),digits=20), type="spearman")
# cor.data <- cor$r
# cor.data[upper.tri(cor.data, diag = T)]<- 0
# pval.data <- cor$P
# pval.data[upper.tri(pval.data, diag = T)]<- NA
# FDR.data <- apply(pval.data,2,p.adjust,method="bonferroni", n = length(pval.data))
# pdf("./pval_bonferonni_hist.pdf")
# hist(pval.data, breaks = 100, col="darkblue")
# dev.off()
# pdf("./FDR_bonferonni_hist.pdf")
# hist(FDR.data, breaks = 100, col="darkblue")
# dev.off()
# pdf("./cor_bonferonni_hist.pdf")
# hist(cor.data, breaks = 10, col="red")
# dev.off()
# cor.data[FDR.data > 0.05]=0
# #load("spear_bonferroni__corrected_cor.data.RData")
# #cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
# library("igraph")
# network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)
# #V(network)$color <- lab[,2]
# print(network)
# 
# library(icesTAF)
# mkdir(paste("./lognormcytokine",fiber_subset,time_point,sep=""))
# 
##########################################################################
#Printing All subgraphs###################################################
##########################################################################
# 
# dg <- decompose.graph(network)
# 
# for (j in 1:length(dg)) {
#   subgraph <- dg[[j]]
#   counter <- 0
#   ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
#   
#   pdf(paste("./lognormcytokine",fiber_subset,time_point,"/",j,fiber_subset,time_point,"_subgraph.pdf",sep=""))
#   par(bg="white", mar=c(0,0,0,0))
#   set.seed(4)
#   plot(subgraph,
#        vertex.size=10,
#        vertex.color=V(subgraph)$color,
#        vertex.label.cex=0.5,
#        vertex.label.color="black",
#        vertex.frame.color="black",
#        #vertex.label = NA,
#        layout = ly
#   )
#   dev.off()
# }
##########################################################################
########################### END ##########################################
##########################################################################

##########################################################################
########################### Proteomics ###################################
##########################################################################



##########################################################################
########################### END ##########################################
##########################################################################

} }
