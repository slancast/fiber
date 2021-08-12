#This will be to make a heapmap of the pvalues for 

week1 <- read.csv("/Users/SLancaster/Desktop/david_results_combined_baseline/Mix.10BaselineAnnotationChart.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
week2 <- read.csv("/Users/SLancaster/Desktop/david_results_combined_baseline/Mix.20BaselineAnnotationChart.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
week3 <- read.csv("/Users/SLancaster/Desktop/david_results_combined_baseline/BaselineMix.30AnnotationChart.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
week4 <- read.csv("/Users/SLancaster/Desktop/david_results_combined_baseline/Mix.WashoutD3BaselineAnnotationChart.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
week5 <- read.csv("/Users/SLancaster/Desktop/david_results_combined_baseline/Mix.WashoutD10BaselineAnnotationChart.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

go_subsets <- c("GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT")
week1 <- week1[week1$Category %in% go_subsets,]
week2 <- week2[week2$Category %in% go_subsets,]
week3 <- week3[week3$Category %in% go_subsets,]
week4 <- week4[week4$Category %in% go_subsets,]
week5 <- week5[week5$Category %in% go_subsets,]

week1_change <- week1[week1$Benjamini < 0.001,]$Term
week2_change <- week2[week2$Benjamini < 0.001,]$Term
week3_change <- week3[week3$Benjamini < 0.001,]$Term

diff_pathways_1 <- Reduce(intersect, list(week1_change, week2_change, week3_change))

week2_change <- week2[week2$Benjamini < 0.001,]$Term
week3_change <- week3[week3$Benjamini < 0.001,]$Term

diff_pathways_2 <- Reduce(intersect, list(week2_change, week3_change))

week3_change <- week3[week3$Benjamini < 0.001,]$Term

diff_pathways_3 <- week3_change

diff_pathways <- Reduce(union, list(diff_pathways_1, diff_pathways_2, diff_pathways_3))


week1_values <- as.matrix(week1[week1$Term %in% diff_pathways,]$Benjamini)
rownames(week1_values) <- week1[week1$Term %in% diff_pathways,]$Term
colnames(week1_values) <- "A10"
week2_values <- as.matrix(week2[week2$Term %in% diff_pathways,]$Benjamini)
rownames(week2_values) <- week2[week2$Term %in% diff_pathways,]$Term
colnames(week2_values) <- "A20"
week3_values <- as.matrix(week3[week3$Term %in% diff_pathways,]$Benjamini)
rownames(week3_values) <- week3[week3$Term %in% diff_pathways,]$Term
colnames(week3_values) <- "A30"
week4_values <- as.matrix(week4[week4$Term %in% diff_pathways,]$Benjamini)
rownames(week4_values) <- week4[week4$Term %in% diff_pathways,]$Term
colnames(week4_values) <- "AW3"
week5_values <- as.matrix(week5[week5$Term %in% diff_pathways,]$Benjamini)
rownames(week5_values) <- week5[week5$Term %in% diff_pathways,]$Term
colnames(week5_values) <- "AW10"


c <- match(diff_pathways, rownames(unique(week1_values)))
d <- diff_pathways[is.na(c)]
A10 <- c()
for (i in is.na(c)) {
  if (isTRUE(i)) {
    A10 <- c(A10, 1)
  }}
A10 <- data.frame(A10)
rownames(A10) <- d
week1_values <- rbind(unique(week1_values), A10)

c <- match(diff_pathways, rownames(unique(week2_values)))
d <- diff_pathways[is.na(c)]
A20 <- c()
for (i in is.na(c)) {
  if (isTRUE(i)) {
    A20 <- c(A20, 1)
  }}
A20 <- data.frame(A20)
rownames(A20) <- d
week2_values <- rbind(unique(week2_values), A20)

c <- match(diff_pathways, rownames(unique(week3_values)))
d <- diff_pathways[is.na(c)]
A30 <- c()
for (i in is.na(c)) {
  if (isTRUE(i)) {
    A30 <- c(A30, 1)
  }}
A30 <- data.frame(A30)
rownames(A30) <- d
week3_values <- rbind(unique(week3_values), A30)

c <- match(diff_pathways, rownames(unique(week4_values)))
d <- diff_pathways[is.na(c)]
AW3 <- c()
for (i in is.na(c)) {
  if (isTRUE(i)) {
    AW3 <- c(AW3, 1)
  }}
AW3 <- data.frame(AW3)
rownames(AW3) <- d
week4_values <- rbind(unique(week4_values), AW3)

c <- match(diff_pathways, rownames(unique(week5_values)))
d <- diff_pathways[is.na(c)]
AW10 <- c()
for (i in is.na(c)) {
  if (isTRUE(i)) {
    AW10 <- c(AW10, 1)
  }}
AW10 <- data.frame(AW10)
rownames(AW10) <- d
week5_values <- rbind(unique(week5_values), AW10)

heatmap_matrix <- cbind(week1_values,week2_values,week3_values,week4_values,week5_values)
heatmap_matrix <- as.matrix(heatmap_matrix)
class(heatmap_matrix) = "numeric"

library(gplots)

heatmap_matrix <- -log10(heatmap_matrix)
pdf('/Users/SLancaster/Desktop/pvalue-heatmap.pdf')
heatmap.2(heatmap_matrix, dendrogram="row", Colv = FALSE, margins = c(5,20), cexRow=.5)
graphics.off()



