#Example data
library(sva)
i <- "/srv/gsfs0/projects/snyder/fiber-ipop/batch_effects/genefamilies-final.pcl"
pcl_df = read.csv(i, sep="\t", header=TRUE, row.names = 1)
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

write.table(combat_edata2, file="/srv/gsfs0/projects/snyder/fiber-ipop/batch_effects/genefamilies-final-combat.pcl", sep="\t")
