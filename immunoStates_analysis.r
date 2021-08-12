#---
#title: "iPOP exercise immunoStates"
#output: html_notebook
#---

#This analysis is to perform cell type deconvolution using immunoStates
#```{r}
library("MetaIntegrator")
library("DESeq2")
library("tidyverse")
#```

#Load dds (DESeqDataSet object) from previously saved
#There is a function for saveRDS that we can apparently use to save the dds file.
#To run this program it looks like I do need the dds file
#```{r}
dds <- readRDS('/home/slancast/RNA/dds.rds')
#```

#Filter out rRNA genes and select timepoints E11, E12, E13, E14, E15, and E19
#```{r}
# remove rRNA gene_ids
# rRNA <- read.table('files/rRNA.v28.txt')
# rRNAkeys <- as.vector(rRNA$V1)
# dds <- dds[!(row.names(dds) %in% rRNAkeys),]

# filter for E11, E12, E13, E14, E15, and E19 Timepoints
# dds <- dds[,row.names(colData(dds)[colData(dds)$Timepoint %in% c("E11", "E12", "E13", "E14", "E15", "E19"),])]
# dds$Timepoint <- relevel(dds$Timepoint, ref = "E11")
# dds$Timepoint <- droplevels(dds$Timepoint)

# remove first trial of 69-001
# dds <- dds [,row.names(colData(dds)[colData(dds)$Patient != "69-001_1",])]
# dds$Patient <- droplevels(dds$Patient)
#```

#Variance stabilizing transformation
#```{r}
vsd <- vst(dds, blind=TRUE)
#```

#Create datasetObject
#```{r}
dataObj <- list()

# $expr
dataObj$expr <- assay(vsd)
# remove version from ensemble id in rownames
rownames(dataObj$expr) <- gsub("\\..*","", rownames(dataObj$expr))

# $keys
dataObj$keys <- ens_ensgID_table[match(rownames(dataObj$expr), ens_ensgID_table$stable_id),]$display_label
names(dataObj$keys) <- rownames(dataObj$expr)

# $formattedName
dataObj$formattedName <- "iPOP fiber"

# $pheno
dataObj$pheno <- as.data.frame(colData(dds))
#```

#```{r}
# checks if it is a valid datasetObject
checkDataObject(dataObj, "Dataset")
#```

#Create metaObject
#```{r}
metaObj <- list()
metaObj$originalData <- list(dataObj)
#```

#Run immunoStates deconvolution
#```{r}
immunoStates <- immunoStatesDecov(metaObj) 
#```

#Results from immunoStatesDecov are stored in $immunoStates
#```{r}
print(immunoStates$immunoStates[[1]])
#```

#Generate plot of immune cell proportions at each timepoint
#```{r}
immunoStates$immunoStates[[1]]$timepoint <- immunoStates$originalData[[1]]$pheno[,"fiber_week"]
pdf("/home/slancast/RNA/immunostates.pdf")
immunoStates$immunoStates[[1]] %>% 
  group_by(timepoint) %>% 
  summarise(
    natural_killer_cell = mean(natural_killer_cell),
    monocyte = mean(monocyte),
    B_cell = mean(B_cell),
    T_cell = mean(T_cell),
    granulocyte = mean(granulocyte)
  ) %>% 
  mutate(other = 1 - (natural_killer_cell + monocyte + B_cell + T_cell + granulocyte)) %>%
  gather(natural_killer_cell, monocyte, B_cell, T_cell, granulocyte, other, key="cell_type", value="proportion") %>%
  mutate(cell_type = fct_relevel(cell_type, "other", "granulocyte", "monocyte", "B_cell", "T_cell", "natural_killer_cell"))%>%
  filter(timepoint %in% c("Baseline","Arabinoxylan.10","Arabinoxylan.20","Arabinoxylan.30", "Arabinoxylan.WashoutD3", "Arabinoxylan.WashoutD10", "LCInulin.10", "LCInulin.20", "LCInulin.30","LCInulin.WashoutD3", "LCInulin.WashoutD10", "Mix.10", "Mix.20", "Mix.30","Mix.WashoutD3","Mix.WashoutD10")) %>%
  mutate(timepoint = fct_relevel(timepoint, "Baseline","Arabinoxylan.10","Arabinoxylan.20","Arabinoxylan.30", "Arabinoxylan.WashoutD3", "Arabinoxylan.WashoutD10", "LCInulin.10", "LCInulin.20", "LCInulin.30","LCInulin.WashoutD3", "LCInulin.WashoutD10", "Mix.10", "Mix.20", "Mix.30","Mix.WashoutD3","Mix.WashoutD10")) %>%
  ggplot() + geom_bar(aes(x=timepoint, y=proportion, fill=cell_type), stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

  dev.off()
  
#```geom_bar(aes(x=timepoint, y=proportion, fill=cell_type), stat="identity") 
#  mutate(cell_type = fct_relevel(cell_type, "other", "granulocyte", "monocyte", "B_cell", "T_cell", "natural_killer_cell"))
# This portion will make one figure with all the cell values in it
# the ,legend.position="none" is how to determine whether or not to use the legend.
#   
#   
#   
#   
  tmp <- immunoStates$immunoStates[[1]] %>% 
    group_by(timepoint) %>% 
    summarise(
      natural_killer_cell = mean(natural_killer_cell),
      monocyte = mean(monocyte),
      B_cell = mean(B_cell),
      T_cell = mean(T_cell),
      granulocyte = mean(granulocyte),
      CD14_positive_monocyte = mean(CD14_positive_monocyte),
      CD16_positive_monocyte = mean(CD16_positive_monocyte),
      CD4_positive_alpha_beta_T_cell = mean(CD4_positive_alpha_beta_T_cell),
      CD56bright_natural_killer_cell = mean(CD56bright_natural_killer_cell),
      CD56dim_natural_killer_cell = mean(CD56dim_natural_killer_cell),
      CD8_positive_alpha_beta_T_cell = mean(CD8_positive_alpha_beta_T_cell),
      MAST_cell = mean(MAST_cell),
      basophil = mean(basophil),
      eosinophil = mean(eosinophil),
      gamma_delta_T_cell = mean(gamma_delta_T_cell),
      hematopoietic_progenitor = mean(hematopoietic_progenitor),
      macrophage_m0 = mean(macrophage_m0),
      macrophage_m1 = mean(macrophage_m1),
      macrophage_m2 = mean(macrophage_m2),
      memory_B_cell = mean(memory_B_cell),
      myeloid_dendritic_cell = mean(myeloid_dendritic_cell),
      naive_B_cell = mean(naive_B_cell),
      neutrophil = mean(neutrophil),
      plasma_cell = mean(plasma_cell),
      plasmacytoid_dendritic_cell = mean(plasmacytoid_dendritic_cell)
    ) %>% 
    
   gather(natural_killer_cell, monocyte, B_cell, T_cell, granulocyte, CD14_positive_monocyte, CD16_positive_monocyte, CD4_positive_alpha_beta_T_cell, CD56bright_natural_killer_cell, CD56dim_natural_killer_cell, CD8_positive_alpha_beta_T_cell, MAST_cell, basophil, eosinophil, gamma_delta_T_cell, hematopoietic_progenitor, macrophage_m0, macrophage_m1, macrophage_m2, memory_B_cell, myeloid_dendritic_cell, naive_B_cell, neutrophil, plasma_cell, plasmacytoid_dendritic_cell, key="cell_type", value="proportion") %>%
    
    print(immunoStates$immunoStates[[1]]$cell_type) %>%
    
    mutate(cell_type = fct_relevel(cell_type, "granulocyte", "monocyte", "B_cell", "T_cell", "natural_killer_cell", "CD14_positive_monocyte", "CD16_positive_monocyte", "CD4_positive_alpha_beta_T_cell", "CD56bright_natural_killer_cell", "CD56dim_natural_killer_cell", "CD8_positive_alpha_beta_T_cell", "MAST_cell", "basophil", "eosinophil", "gamma_delta_T_cell", "hematopoietic_progenitor", "macrophage_m0", "macrophage_m1", "macrophage_m2", "memory_B_cell", "myeloid_dendritic_cell", "naive_B_cell", "neutrophil", "plasma_cell", "plasmacytoid_dendritic_cell"))%>%
    
    filter(timepoint %in% c("Baseline","Arabinoxylan.10","Arabinoxylan.20","Arabinoxylan.30", "Arabinoxylan.WashoutD3", "Arabinoxylan.WashoutD10", "LCInulin.10", "LCInulin.20", "LCInulin.30","LCInulin.WashoutD3", "LCInulin.WashoutD10", "Mix.10", "Mix.20", "Mix.30","Mix.WashoutD3","Mix.WashoutD10")) %>%
    
    mutate(timepoint = fct_relevel(timepoint, "Baseline","Arabinoxylan.10","Arabinoxylan.20","Arabinoxylan.30", "Arabinoxylan.WashoutD3", "Arabinoxylan.WashoutD10", "LCInulin.10", "LCInulin.20", "LCInulin.30","LCInulin.WashoutD3", "LCInulin.WashoutD10", "Mix.10", "Mix.20", "Mix.30","Mix.WashoutD3","Mix.WashoutD10")) %>%
    
    ggplot() + geom_bar(aes(x=timepoint, y=proportion, fill=cell_type), stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #,legend.position="none"
  
  ggsave(file=paste("/home/slancast/RNA/immunostates.pdf",sep=""), plot=tmp)
    
    
  
 # 
 #  
 #  
 #  
 #  
 #  
 #  
 #   This will create the table using tidyverse that has the summary statistics for the barplot. This table will feed into the loop below
 #   (which doesn't use tidyverse) to create all the table for the indivudal figures. To this I will add the standard deviations as well
 #   From Eric's code
  
  celltype_mean <- immunoStates$immunoStates[[1]] %>% 
    group_by(timepoint) %>% 
    summarise(
      natural_killer_cell = mean(natural_killer_cell),
      monocyte = mean(monocyte),
      B_cell = mean(B_cell),
      T_cell = mean(T_cell),
      granulocyte = mean(granulocyte),
      CD14_positive_monocyte = mean(CD14_positive_monocyte),
      CD16_positive_monocyte = mean(CD16_positive_monocyte),
      CD4_positive_alpha_beta_T_cell = mean(CD4_positive_alpha_beta_T_cell),
      CD56bright_natural_killer_cell = mean(CD56bright_natural_killer_cell),
      CD56dim_natural_killer_cell = mean(CD56dim_natural_killer_cell),
      CD8_positive_alpha_beta_T_cell = mean(CD8_positive_alpha_beta_T_cell),
      MAST_cell = mean(MAST_cell),
      basophil = mean(basophil),
      eosinophil = mean(eosinophil),
      gamma_delta_T_cell = mean(gamma_delta_T_cell),
      hematopoietic_progenitor = mean(hematopoietic_progenitor),
      macrophage_m0 = mean(macrophage_m0),
      macrophage_m1 = mean(macrophage_m1),
      macrophage_m2 = mean(macrophage_m2),
      memory_B_cell = mean(memory_B_cell),
      myeloid_dendritic_cell = mean(myeloid_dendritic_cell),
      naive_B_cell = mean(naive_B_cell),
      neutrophil = mean(neutrophil),
      plasma_cell = mean(plasma_cell),
      plasmacytoid_dendritic_cell = mean(plasmacytoid_dendritic_cell)
    ) %>% 
    
    gather(natural_killer_cell, monocyte, B_cell, T_cell, granulocyte, CD14_positive_monocyte, CD16_positive_monocyte, CD4_positive_alpha_beta_T_cell, CD56bright_natural_killer_cell, CD56dim_natural_killer_cell, CD8_positive_alpha_beta_T_cell, MAST_cell, basophil, eosinophil, gamma_delta_T_cell, hematopoietic_progenitor, macrophage_m0, macrophage_m1, macrophage_m2, memory_B_cell, myeloid_dendritic_cell, naive_B_cell, neutrophil, plasma_cell, plasmacytoid_dendritic_cell, RMSE, key="cell_type", value="proportion")
    
  # Second portion adding the standard deviations to see how the error bars look.
  
  celltype_se <- immunoStates$immunoStates[[1]] %>% 
    group_by(timepoint) %>% 
    summarise(
      natural_killer_cell = sd(natural_killer_cell) / sqrt(n()) * 1.96,
      monocyte = sd(monocyte) / sqrt(n()) * 1.96,
      B_cell = sd(B_cell) / sqrt(n()) * 1.96,
      T_cell = sd(T_cell) / sqrt(n()) * 1.96,
      granulocyte = sd(granulocyte) / sqrt(n()) * 1.96,
      CD14_positive_monocyte = sd(CD14_positive_monocyte) / sqrt(n()) * 1.96,
      CD16_positive_monocyte = sd(CD16_positive_monocyte) / sqrt(n()) * 1.96,
      CD4_positive_alpha_beta_T_cell = sd(CD4_positive_alpha_beta_T_cell) / sqrt(n()) * 1.96,
      CD56bright_natural_killer_cell = sd(CD56bright_natural_killer_cell) / sqrt(n()) * 1.96,
      CD56dim_natural_killer_cell = sd(CD56dim_natural_killer_cell) / sqrt(n()) * 1.96,
      CD8_positive_alpha_beta_T_cell = sd(CD8_positive_alpha_beta_T_cell) / sqrt(n()) * 1.96,
      MAST_cell = sd(MAST_cell) / sqrt(n()) * 1.96,
      basophil = sd(basophil) / sqrt(n()) * 1.96,
      eosinophil = sd(eosinophil) / sqrt(n()) * 1.96,
      gamma_delta_T_cell = sd(gamma_delta_T_cell) / sqrt(n()) * 1.96,
      hematopoietic_progenitor = sd(hematopoietic_progenitor) / sqrt(n()) * 1.96,
      macrophage_m0 = sd(macrophage_m0) / sqrt(n()) * 1.96,
      macrophage_m1 = sd(macrophage_m1) / sqrt(n()) * 1.96,
      macrophage_m2 = sd(macrophage_m2) / sqrt(n()) * 1.96,
      memory_B_cell = sd(memory_B_cell) / sqrt(n()) * 1.96,
      myeloid_dendritic_cell = sd(myeloid_dendritic_cell) / sqrt(n()) * 1.96,
      naive_B_cell = sd(naive_B_cell) / sqrt(n()) * 1.96,
      neutrophil = sd(neutrophil) / sqrt(n()) * 1.96,
      plasma_cell = sd(plasma_cell) / sqrt(n()) * 1.96,
      plasmacytoid_dendritic_cell = sd(plasmacytoid_dendritic_cell) / sqrt(n()) * 1.96
    ) %>%   
    
    gather(natural_killer_cell, monocyte, B_cell, T_cell, granulocyte, CD14_positive_monocyte, CD16_positive_monocyte, CD4_positive_alpha_beta_T_cell, CD56bright_natural_killer_cell, CD56dim_natural_killer_cell, CD8_positive_alpha_beta_T_cell, MAST_cell, basophil, eosinophil, gamma_delta_T_cell, hematopoietic_progenitor, macrophage_m0, macrophage_m1, macrophage_m2, memory_B_cell, myeloid_dendritic_cell, naive_B_cell, neutrophil, plasma_cell, plasmacytoid_dendritic_cell, key="cell_type", value="sd")
    
  
  
    tmp <- full_join(celltype_mean, celltype_se) %>%
    mutate(cell_type = fct_relevel(cell_type, "granulocyte", "monocyte", "B_cell", "T_cell", "natural_killer_cell", "CD14_positive_monocyte", "CD16_positive_monocyte", "CD4_positive_alpha_beta_T_cell", "CD56bright_natural_killer_cell", "CD56dim_natural_killer_cell", "CD8_positive_alpha_beta_T_cell", "MAST_cell", "basophil", "eosinophil", "gamma_delta_T_cell", "hematopoietic_progenitor", "macrophage_m0", "macrophage_m1", "macrophage_m2", "memory_B_cell", "myeloid_dendritic_cell", "naive_B_cell", "neutrophil", "plasma_cell", "plasmacytoid_dendritic_cell"))%>%
    filter(timepoint %in% c("Baseline","Arabinoxylan.10","Arabinoxylan.20","Arabinoxylan.30", "Arabinoxylan.WashoutD3", "Arabinoxylan.WashoutD10", "LCInulin.10", "LCInulin.20", "LCInulin.30","LCInulin.WashoutD3", "LCInulin.WashoutD10", "Mix.10", "Mix.20", "Mix.30","Mix.WashoutD3","Mix.WashoutD10")) %>%
    mutate(timepoint = fct_relevel(timepoint, "Baseline","Arabinoxylan.10","Arabinoxylan.20","Arabinoxylan.30", "Arabinoxylan.WashoutD3", "Arabinoxylan.WashoutD10", "LCInulin.10", "LCInulin.20", "LCInulin.30","LCInulin.WashoutD3", "LCInulin.WashoutD10", "Mix.10", "Mix.20", "Mix.30","Mix.WashoutD3","Mix.WashoutD10"))
    
    
    
for (i in colnames(immunoStates$immunoStates[[1]])) {
  print(i)
  tmpa <- tmp[tmp$cell_type==i,]
  gpl <- ggplot(tmpa, aes(x=timepoint, y=proportion, color=cell_type, group=cell_type)) + 
    geom_line(stat="identity", position=position_dodge(1)) + 
    geom_point(stat="identity", position=position_dodge(1)) +  
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.5, position=position_dodge(1), alpha=0.5)

  ggsave(file=paste("/home/slancast/RNA/immunostates",i,".pdf",sep=""))
}


