dir <- "/Users/SLancaster/Desktop/Projects/Fiber/Data/rsem_renamed/"

#Example data
library(sva)
library(bladderbatch)
data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)

#The model below is for latent variables
mod = model.matrix(~as.factor(cancer),data=pheno)
mod0 = model.matrix(~1, data=pheno) 
#n.sv = num.sv(edata,mod,method="leek") #This is used to estimate the number of latent factors, if I'm sure of the number of factors I can determine that too.
sva1 = sva(edata,mod,mod0,n.sv=2)
expdata <- resdata[8:ncol(resdata)+9]
expdata <- as.matrix(expdata)
row.names(expdata) <- resdata[,1]
phendata <- samples
row.names(phendata) <- colnames(expdata)
mod = model.matrix(~batch,data=phendata)
mod0 = model.matrix(~1, data=phendata)
sva1 = sva(expdata,mod,mod0,n.sv=2)
summary(lm(sva1$sv ~ pheno$batch))

#This code is for known batches
#I think for fiber I should include fiber and week in the model matrix
#Unfortunately the example doesn't show this, but when googling it, it seems that this is the way to go.
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

##########################
###########PCL############
##########################
#,"../biom_files/genefamilies/genefamilies" #needs to be run on gcloud
for (i in c("pathabundance","pathcoverage","bugs-list","genefamilies-90ko-groupedonly")) {
  pcl_df = read.csv(paste("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/",i,"-final.pcl",sep=""), sep="\t", header=TRUE, row.names = 1)
  print(i)
  #Full#
  
  pcl_metadata_rows = 5
  pcl_metadata = head(pcl_df,pcl_metadata_rows)
  
  
  pcl_df2 <- tail(pcl_df, -pcl_metadata_rows) #Now getting rid of the metadata
  pcl_df2 <- data.frame(lapply(pcl_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(pcl_df2)) #To change to a double, this needs to go through character first
  pcl_df2 <- data.frame(lapply(pcl_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(pcl_df2)) #then numeric
  pcl_metadata = data.frame(t(pcl_metadata))
  pcl_df2 <- as.matrix(pcl_df2)
  pcl_df2 <- pcl_df2[apply(pcl_df2, 1, var) != 0,]
  pcl_df2 <- log(pcl_df2+1,2)
  pcl_df2 <- na.omit(pcl_df2)
  
  batch <- pcl_metadata$Metagenomic_seq_plate
  modcombat <- cbind(as.numeric(pcl_metadata$Fiber),as.numeric(pcl_metadata$Dose),as.numeric(pcl_metadata$Participant))
  rownames(modcombat) <- rownames(pcl_metadata)
  combat_edata = ComBat(dat=pcl_df2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  class(combat_edata) <- "character"
  
  pcl_metadata = data.frame(t(pcl_metadata))
 
  combat_edata2 <- rbind(pcl_metadata,combat_edata)
  
  write.table(combat_edata2, file=paste("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/",i,"-final-combat.pcl",sep=""), sep="\t")

}

