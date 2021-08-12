fiber_subsets = c("Arabinoxylan","LCInulin","Mix")
fdr_method = "BH"

for (fiber_subset in fiber_subsets) {
  library(icesTAF)
  mkdir(paste("~/lognormmultiomics",fiber_subset,fdr_method, sep=""))
  
  load(paste("~/NormAggreg_Log_",fiber_subset,"_cytokine_df.RData",sep=""))
  load(paste("~/NormAggreg_Log_",fiber_subset,"_pathabundance_df.RData",sep=""))
  #load(paste("~/NormAggreg_Log_",fiber_subset,"_RNA_df.RData",sep=""))
  load(paste("~/NormAggreg_Log_",fiber_subset,"_metabolomics_df.RData",sep=""))
  load(paste("~/NormAggreg_Log_",fiber_subset,"_lipids_df.RData",sep=""))
  load(paste("~/NormAggreg_Log_",fiber_subset,"_proteomics_df.RData",sep=""))
  load(paste("~/NormAggreg_Log_",fiber_subset,"species_metaphlan_df.RData",sep=""))
  load(paste("~/Aggregate_",fiber_subset,"_clinical_df.RData",sep=""))
  
  m <- rbind(t(logaggregnorm_pathabundance), t(logaggregnorm_lipids), t(logaggregnorm_metabolomics), t(logaggregnorm_cytokine), t(logaggregnorm_proteomics), t(logaggregnorm_metaphlan), t(aggregated_clinical_df))
  lab <- c(rep("Microbiome",nrow(t(logaggregnorm_pathabundance))),rep("lipids",nrow(t(logaggregnorm_lipids))),rep("Metabolites",nrow(t(logaggregnorm_metabolomics))),rep("Cytokine",nrow(t(logaggregnorm_cytokine))),rep("proteomics",nrow(t(logaggregnorm_proteomics))),rep("metaphlan",nrow(t(logaggregnorm_metaphlan))),rep("clinicals",nrow(t(aggregated_clinical_df))))
  color <- c(rep("blue",nrow(t(logaggregnorm_pathabundance))),rep("red",nrow(t(logaggregnorm_lipids))),rep("yellow",nrow(t(logaggregnorm_metabolomics))),rep("orange",nrow(t(logaggregnorm_cytokine))),rep("green",nrow(t(logaggregnorm_proteomics))),rep("blue",nrow(t(logaggregnorm_metaphlan))),rep("orange",nrow(t(aggregated_clinical_df))))
  lab <- cbind(lab,color)
  rownames(lab) <- rownames(m)
  
  
  
  library(network)
  library(Mfuzz)
  library(matrixStats)
  n <- as.matrix(m)
  class(n) <- "numeric"
  eset1 <- ExpressionSet(n)
  eset1 <- standardise(eset1) #Running standarise 
  o <- exprs(eset1)
  o <- na.omit(o)
  lab <- lab[which(rownames(lab) %in% rownames(o)),]

  library("Hmisc")
  # Plot correlation graph
  
  cor <- rcorr(format(t(o),digits=20), type="spearman")
  cor.data <- cor$r
  cor.data[upper.tri(cor.data, diag = T)]<- 0
  pval.data <- cor$P
  pval.data[upper.tri(pval.data, diag = T)]<- NA
  FDR.data <- apply(pval.data,2,p.adjust,method=fdr_method, n = length(pval.data))
  pdf(paste("~/lognormmultiomics",fiber_subset,fdr_method,"/pval_hist.pdf",sep=""))
  hist(pval.data, breaks = 100, col="darkblue")
  dev.off()
  pdf(paste("~/lognormmultiomics",fiber_subset,fdr_method,"/FDR_hist.pdf",sep=""))
  hist(FDR.data, breaks = 100, col="darkblue")
  dev.off()
  pdf(paste("~/lognormmultiomics",fiber_subset,fdr_method,"/cor_hist.pdf",sep=""))
  hist(cor.data, breaks = 10, col="red")
  dev.off()
  cor.data[FDR.data > 0.05]=0
  write.table(cor.data,file=paste("~/lognormmultiomics",fiber_subset,fdr_method,"/cor_data.txt",sep=""), sep="\t")
  #load("spear_bonferroni__corrected_cor.data.RData")
  #cor.data <- cor.data[which(rowSums(cor.data)>0),which(colSums(cor.data)>0)]
  library("igraph")
  network=graph.adjacency(cor.data, weighted=T, mode="undirected", diag=F)
  
  
  V(network)$color <- lab[,2]
  print(network)
  
  ly <- layout.fruchterman.reingold(network,dim=2,grid="nogrid")

  pdf(paste("~/lognormmultiomics",fiber_subset,fdr_method,"/Full_network.pdf",sep=""))
  par(bg="white", mar=c(0,0,0,0))
  set.seed(4)
  plot(network,
       vertex.size=0.3,
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
  subgraph_composition <- data.frame(Name=character(),
                   Subgraph=integer(),
                   stringsAsFactors=FALSE)
  
  for (j in 1:length(dg)) {
    subgraph <- dg[[j]]
    counter <- 0
    ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")
    
    pdf(paste("~/lognormmultiomics",fiber_subset,fdr_method,"/",j,fiber_subset,"_subgraph.pdf",sep=""))
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
    
    #I want to write out all the subgraphs by the vertexes too.
    
    vertex_list <- data.frame(V(subgraph)$name)
    vertex_list <- cbind(vertex_list,rep(j, length(vertex_list)))
    
    subgraph_composition<- rbind(subgraph_composition,vertex_list)
  }
  
  write.table(subgraph_composition, file=paste("~/lognormmultiomics",fiber_subset,fdr_method,"/",fiber_subset,"_subgraph_composition.txt",sep=""),sep="\t")
  
}


  