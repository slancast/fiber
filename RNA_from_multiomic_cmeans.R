#!/usr/bin/r

#This program will be to analyze the transctipt data from the multiomic clusters
#First it will pull them out, and then it will add the gene symbol and perform
#enrichments etc.

fiber_subset <- "Arabinoxylan"
df <- read.csv(paste0("/Users/SLancaster/Desktop/multiomics_clusters/clusters16multiomics_",fiber_subset,"_cluster_membership.txt"), header = TRUE, sep="\t", row.names = NULL)

df2 <- df[which(startsWith(df[,1], "ENSG")),] #Pulling out only the ensemble genes

#running the program I wrote called entrezid
source('~/Desktop/Projects/Fiber/Multiomics/Software/Utils.r')
df2 <- entrezid(df2)

interesting_clusters <- c(3,6,13,15) #Put all the intersting clusters here and then it will iterate over them. Otherwise you can run them individually

for (cluster_number in interesting_clusters ) {
#pulling out the cluster of choice. This coule be updated to use all the clusters if necessary
#cluster_number <- 3 #Use this line instead of loop to run clusters one by one
df3 <- df2[which(df2[,2] %in% cluster_number),]

#setting up davidwebservice and the background
library(RDAVIDWebService)
david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setTimeOut(david, 1000000) #Large for large gene lists
#setAnnotationCategories(david, c("GOTERM_BP_DIRECT","GOTERM_MF_DIRECT","GOTERM_CC_DIRECT")) #This specficies which annotation category I want (i.e. GO, Kegg, UP, etc)
#The previous function doesn't seem to increase my power with the test set I was using, although I remember a proir test where it did increase the power.

#It seems that entrez works slightly better than ensembl but they are comparable. Symbol wokred the worst.
background <- addList(david, df2$ENTREZID, idType="ENTREZ_GENE_ID",listName="Total_genes", listType="Background")  #This needs to be reset before adding new result lists

#I will use this file name several times, so creating it as a object
AnnotChartFilename <- paste("/Users/SLancaster/Desktop/multiomics_clusters/cluster_",cluster_number,"_",fiber_subset,"AnnotationChart.txt",sep="")

result<-addList(david, df3$ENTREZID, idType="ENTREZ_GENE_ID",listName=paste("cluster",cluster_number,sep=""), listType="Gene")
getFunctionalAnnotationChartFile(david,  fileName = AnnotChartFilename)
getFunctionalAnnotationTableFile(david, fileName = paste("/Users/SLancaster/Desktop/multiomics_clusters/cluster_",cluster_number,"_",fiber_subset,"AnnotationTable.txt",sep=""))
getClusterReportFile(david, fileName = paste("/Users/SLancaster/Desktop/multiomics_clusters/cluster_",cluster_number,"_",fiber_subset,"ClusterReport.txt",sep=""))

#This section will readjust the pvalues. David does a terrible job of that.
chart <- read.csv(AnnotChartFilename, header = TRUE, sep="\t", row.names = NULL)
padjusted <- p.adjust(chart$PValue, method = "BH")
chart$Benjamini <- padjusted
write.table(chart, file=AnnotChartFilename, sep="\t", row.names = FALSE)

}

