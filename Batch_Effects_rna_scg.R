#Example data
library(sva)
i <- "/labs/mpsnyder/slancast/fiber/rnaseq/fall2019/batch_correction/rna_raw_counts.csv"
rna_df = read.csv(i, sep=",", header=TRUE, row.names = 1)
print(i)
#Full#

rna_metadata_rows = 5
rna_metadata = head(rna_df,rna_metadata_rows)


rna_df2 <- tail(rna_df, -rna_metadata_rows) #Now getting rid of the metadata
rna_df2 <- data.frame(lapply(rna_df2, as.character), stringsAsFactors=FALSE, row.names = rownames(rna_df2)) #To change to a double, this needs to go through character first
rna_df2 <- data.frame(lapply(rna_df2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(rna_df2)) #then numeric
rna_metadata = data.frame(t(rna_metadata))
rna_df2 <- as.matrix(rna_df2)
rna_df2 <- rna_df2[apply(rna_df2, 1, var) != 0,]
rna_df2 <- log(rna_df2+1,2)
rna_df2 <- na.omit(rna_df2)

batch <- rna_metadata$batch
modcombat <- cbind(as.numeric(rna_metadata$condition),as.numeric(rna_metadata$week),as.numeric(rna_metadata$participant))
rownames(modcombat) <- rownames(rna_metadata)
combat_edata = ComBat(dat=rna_df2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
class(combat_edata) <- "character"

rna_metadata = data.frame(t(rna_metadata))

combat_edata2 <- rbind(rna_metadata,combat_edata)

write.table(combat_edata2, file="/labs/mpsnyder/slancast/fiber/rnaseq/fall2019/batch_correction/rna-final-combat.txt", sep="\t")
