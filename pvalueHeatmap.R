#In this program I'll create a script to create a pvalue histogram
#Kevin made one of these for the exercise data, and I think it would be good for this project too
#Especially without the RNA data complete, I don't think I can get a final dataset,
#but it might be with some of the DESeq data we alraedy have

##This data will match up the genes
resdata1 = read.table("/Volumes/My_Book_Duo/RNAseq_results_scg/RNAseq_results/fiber_week_Arabinoxylan.Baseline_vs_Arabinoxylan.20-diffexpr-resultsb.csv", sep=",", header=TRUE,row.names=1)
resdata2 = read.table("/Volumes/My_Book_Duo/RNAseq_results_scg/RNAseq_results/fiber_week_Arabinoxylan.Baseline_vs_Arabinoxylan.30-diffexpr-resultsb.csv", sep=",", header=TRUE,row.names=1)
Summary_Data <- cbind(Gene1 = as.character(resdata1$Gene), pvalue1 = resdata1$padj)
Summary_Data <- data.frame(Summary_Data)
Gene2 <- resdata2$Gene[match(resdata1$Gene, resdata2$Gene)]
pvalue2 <- resdata2$padj[match(resdata1$Gene, resdata2$Gene)]
Summary_Data2 <- data.frame(cbind(Gene2 = as.character(Gene2), pvalue2 = pvalue1))
Summary_Data <- data.frame(cbind(Summary_Data, pvalue2))
rownames(Summary_Data) <- Summary_Data$Gene1
Summary_Data <- data.frame(Summary_Data)

library(reshape2)
Summary_Data.melted <- melt(Summary_Data[1:20,])
p <- ggplot(Summary_Data.melted, aes(y = Gene, x = variable, fill = value)) + geom_tile()

###This portion is to use the heatmap function rather than ggplot
if (FALSE) {
Summary_Data <- Summary_Data[,-1]
Summary_Data <- as.matrix(Summary_Data)
class(Summary_Data) <- "numeric"
heatmap(Summary_Data)
Summary_Data <- data.frame(cbind(Gene = rownames(Summary_Data), Summary_Data))
Summary_Data <- data.frame(cbind(Gene = rownames(Summary_Data), Summary_Data))
}

