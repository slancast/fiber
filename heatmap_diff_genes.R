##In this program I will take the most differentially expressed genes by log fold 
##change that are in the arabinoxylan treatments, and then plot a heatmap of those

week1 <- read.csv("/Users/SLancaster/Desktop/deseq_results_combined_baseline/Arabinoxylan.10Baseline-diffexpr-resultsb.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
week2 <- read.csv("/Users/SLancaster/Desktop/deseq_results_combined_baseline/Arabinoxylan.20Baseline-diffexpr-resultsb.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
week3 <- read.csv("/Users/SLancaster/Desktop/deseq_results_combined_baseline/Arabinoxylan.30Baseline-diffexpr-resultsb.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
week4 <- read.csv("/Users/SLancaster/Desktop/deseq_results_combined_baseline/Arabinoxylan.WashoutD3Baseline-diffexpr-resultsb.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
week5 <- read.csv("/Users/SLancaster/Desktop/deseq_results_combined_baseline/Arabinoxylan.WashoutD10Baseline-diffexpr-resultsb.csv", sep=",", row.names = 1, stringsAsFactors=FALSE, header=TRUE)

#There are 6 caegories of expression we're looking for:
#1
week1_change <- week1[week1$log2FoldChange > 0.5,]$Gene
week2_change <- week2[week2$log2FoldChange > 0.5,]$Gene
week3_change <- week3[week3$log2FoldChange > 0.5,]$Gene

diff_genes_1 <- Reduce(intersect, list(week1_change, week2_change, week3_change))

#2
week1_change <- week1[week1$log2FoldChange > -0.25 & week1$log2FoldChange < 0.5,]$Gene
week2_change <- week2[week2$log2FoldChange > 0.5,]$Gene
week3_change <- week3[week3$log2FoldChange > 0.5,]$Gene

diff_genes_2 <- Reduce(intersect, list(week1_change, week2_change, week3_change))

if (FALSE) {
#3
week3_change <- week3[week3$log2FoldChange > 0.5,]$Gene
week2_change <- week2[abs(week2$log2FoldChange) < 0.25,]$Gene
week1_change <- week1[abs(week1$log2FoldChange) < 0.25,]$Gene

diff_genes_3 <- Reduce(intersect, list(week1_change, week2_change, week3_change))
}

#4
week1_change <- week1[week1$log2FoldChange < -0.5,]$Gene
week2_change <- week2[week2$log2FoldChange < -0.5,]$Gene
week3_change <- week3[week3$log2FoldChange < -0.5,]$Gene

diff_genes_4 <- Reduce(intersect, list(week1_change, week2_change, week3_change))

#5
week1_change <- week1[week1$log2FoldChange < 0.25 & week1$log2FoldChange > -0.5,]$Gene
week2_change <- week2[week2$log2FoldChange < -0.5,]$Gene
week3_change <- week3[week3$log2FoldChange < -0.5,]$Gene

diff_genes_5 <- Reduce(intersect, list(week1_change, week2_change, week3_change))

if (FALSE) {
#6
week3_change <- week3[week3$log2FoldChange < -0.5,]$Gene
week2_change <- week2[abs(week2$log2FoldChange) < 0.25,]$Gene
week1_change <- week1[abs(week1$log2FoldChange) < 0.25,]$Gene

diff_genes_6 <- Reduce(intersect, list(week1_change, week2_change, week3_change))
}

diff_genes <- c(diff_genes_1, diff_genes_2, diff_genes_4, diff_genes_5)

week1_values <- as.matrix(week1[week1$Gene %in% diff_genes_1,]$log2FoldChange)
rownames(week1_values) <- week1[week1$Gene %in% diff_genes_1,]$Gene
week2_values <- as.matrix(week2[week2$Gene %in% diff_genes_1,]$log2FoldChange)
rownames(week2_values) <- week2[week2$Gene %in% diff_genes_1,]$Gene
week3_values <- as.matrix(week3[week3$Gene %in% diff_genes_1,]$log2FoldChange)
rownames(week3_values) <- week3[week3$Gene %in% diff_genes_1,]$Gene
week4_values <- as.matrix(week4[week4$Gene %in% diff_genes_1,]$log2FoldChange)
rownames(week4_values) <- week4[week4$Gene %in% diff_genes_1,]$Gene
week5_values <- as.matrix(week5[week5$Gene %in% diff_genes_1,]$log2FoldChange)
rownames(week5_values) <- week5[week5$Gene %in% diff_genes_1,]$Gene

heatmap_matrix1 = cbind(week1_values,week2_values,week3_values,week4_values,week5_values)

week1_values <- as.matrix(week1[week1$Gene %in% diff_genes_2,]$log2FoldChange)
rownames(week1_values) <- week1[week1$Gene %in% diff_genes_2,]$Gene
week2_values <- as.matrix(week2[week2$Gene %in% diff_genes_2,]$log2FoldChange)
rownames(week2_values) <- week2[week2$Gene %in% diff_genes_2,]$Gene
week3_values <- as.matrix(week3[week3$Gene %in% diff_genes_2,]$log2FoldChange)
rownames(week3_values) <- week3[week3$Gene %in% diff_genes_2,]$Gene
week4_values <- as.matrix(week4[week4$Gene %in% diff_genes_2,]$log2FoldChange)
rownames(week4_values) <- week4[week4$Gene %in% diff_genes_2,]$Gene
week5_values <- as.matrix(week5[week5$Gene %in% diff_genes_2,]$log2FoldChange)
rownames(week5_values) <- week5[week5$Gene %in% diff_genes_2,]$Gene

heatmap_matrix2 = cbind(week1_values,week2_values,week3_values,week4_values,week5_values)

week1_values <- as.matrix(week1[week1$Gene %in% diff_genes_4,]$log2FoldChange)
rownames(week1_values) <- week1[week1$Gene %in% diff_genes_4,]$Gene
week2_values <- as.matrix(week2[week2$Gene %in% diff_genes_4,]$log2FoldChange)
rownames(week2_values) <- week2[week2$Gene %in% diff_genes_4,]$Gene
week3_values <- as.matrix(week3[week3$Gene %in% diff_genes_4,]$log2FoldChange)
rownames(week3_values) <- week3[week3$Gene %in% diff_genes_4,]$Gene
week4_values <- as.matrix(week4[week4$Gene %in% diff_genes_4,]$log2FoldChange)
rownames(week4_values) <- week4[week4$Gene %in% diff_genes_4,]$Gene
week5_values <- as.matrix(week5[week5$Gene %in% diff_genes_4,]$log2FoldChange)
rownames(week5_values) <- week5[week5$Gene %in% diff_genes_4,]$Gene

heatmap_matrix4 = cbind(week1_values,week2_values,week3_values,week4_values,week5_values)

week1_values <- as.matrix(week1[week1$Gene %in% diff_genes_5,]$log2FoldChange)
rownames(week1_values) <- week1[week1$Gene %in% diff_genes_5,]$Gene
week2_values <- as.matrix(week2[week2$Gene %in% diff_genes_5,]$log2FoldChange)
rownames(week2_values) <- week2[week2$Gene %in% diff_genes_5,]$Gene
week3_values <- as.matrix(week3[week3$Gene %in% diff_genes_5,]$log2FoldChange)
rownames(week3_values) <- week3[week3$Gene %in% diff_genes_5,]$Gene
week4_values <- as.matrix(week4[week4$Gene %in% diff_genes_5,]$log2FoldChange)
rownames(week4_values) <- week4[week4$Gene %in% diff_genes_5,]$Gene
week5_values <- as.matrix(week5[week5$Gene %in% diff_genes_5,]$log2FoldChange)
rownames(week5_values) <- week5[week5$Gene %in% diff_genes_5,]$Gene

heatmap_matrix5 = cbind(week1_values,week2_values,week3_values,week4_values,week5_values)

library(gplots)
heatmap_matrix1 <- heatmap_matrix1[-match("gSpikein_phiX174",rownames(heatmap_matrix1)),]
rownames(heatmap_matrix1) <- gsub("\\..*","",rownames(heatmap_matrix1))
pdf("/Users/SLancaster/Desktop/diff_expressed_gene_heatmaps/Ax_heatmap_up1.pdf")
heatmap.2(heatmap_matrix1, Colv = FALSE, Rowv = FALSE, col=bluered(75), dendrogram = 'none', margins = c(5,7), cexRow=.5, trace = 'none')
dev.off()

rownames(heatmap_matrix2) <- gsub("\\..*","",rownames(heatmap_matrix2))
pdf("/Users/SLancaster/Desktop/diff_expressed_gene_heatmaps/Ax_heatmap_up2.pdf")
heatmap.2(heatmap_matrix2, Colv = FALSE, Rowv = FALSE, col=bluered(75), dendrogram = 'none', margins = c(5,7), cexRow=.5, trace = 'none')
dev.off()

heatmap_matrix_down <- rbind(heatmap_matrix4, heatmap_matrix5)
rownames(heatmap_matrix_down) <- gsub("\\..*","",rownames(heatmap_matrix_down))
pdf("/Users/SLancaster/Desktop/diff_expressed_gene_heatmaps/Ax_heatmap_down.pdf")
heatmap.2(heatmap_matrix_down, Colv = FALSE, Rowv = FALSE, col=bluered(75), dendrogram = 'none', margins = c(5,7), cexRow=.5, trace = 'none')
dev.off()

library(RDAVIDWebService)
david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setTimeOut(david, 100000)
background_genes <- gsub("\\..*","",week1$Gene)
background <- addList(david, background_genes, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")
result <- addList(david, c(rownames(heatmap_matrix1), rownames(heatmap_matrix2)), idType="ENSEMBL_GENE_ID",listName="heatmap_up", listType="Gene")
getFunctionalAnnotationChartFile(david, fileName = paste("/Users/SLancaster/Desktop/diff_expressed_gene_heatmaps/diff_genes_heatmap_AXup_AnnotationChart.txt",sep=""))

david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setTimeOut(david, 100000)
background_genes <- gsub("\\..*","",week1$Gene)
background <- addList(david, background_genes, idType="ENSEMBL_GENE_ID",listName="Total_genes", listType="Background")
result <- addList(david, rownames(heatmap_matrix_down), idType="ENSEMBL_GENE_ID",listName="heatmap_down", listType="Gene")
getFunctionalAnnotationChartFile(david, fileName = paste("/Users/SLancaster/Desktop/diff_expressed_gene_heatmaps/diff_genes_heatmap_AXdown_AnnotationChart.txt",sep=""))

