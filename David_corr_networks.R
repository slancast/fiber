#This will be to perform enrichments in the correlation networks on the cluster. 
#It's difficult to install David on ubuntu, so they will all be performed in the cluster snapshot "multiomic-clustering"
#Even though it doesn't actually use any of the soft clustering.

for (sample in c("Arabinoxylan10", 
                 "Arabinoxylan20", 
                 "Arabinoxylan30",
                 "ArabinoxylanBaseline",
                 "ArabinoxylanWashoutD3",
                 "ArabinoxylanWashoutD10",
                 "LCInulin10", 
                 "LCInulin20", 
                 "LCInulin30",
                 "LCInulinBaseline",
                 "LCInulinWashoutD3",
                 "LCInulinWashoutD10",
                 "Mix10", 
                 "Mix20", 
                 "Mix30",
                 "MixBaseline",
                 "MixWashoutD3",
                 "MixWashoutD10")) {

results_files <- list.files(path = paste("/home/slancast/RNA_correlation_networks/lognormrna",sample,sep=""), pattern=".txt")
library(RDAVIDWebService)

pat <- paste("/home/slancast/RNA_correlation_networks/lognormrna",sample,"/",sample,"_total_genes.txt",sep="")
background_file = read.csv(pat, sep="\t", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
background_genes <- colnames(background_file)

david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setTimeOut(david, 1000000)
background <- addList(david, background_genes, idType="ENSEMBL_GENE_ID",listName="Total_background_genes", listType="Background")

for (i in results_files) { 
  if (i == pat) next
  name = substr(i,1,nchar(i)-4)
  result_file = read.csv(paste("/home/slancast/RNA_correlation_networks/lognormrna",sample,"/",i,sep=""), sep="\t", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
  sig_genes <- colnames(result_file)
  if (length(sig_genes) <= 2) next
  print(name)
  result<- tryCatch({addList(david, sig_genes, idType="ENSEMBL_GENE_ID",listName= name, listType="Gene")},
                    error = function(err){c(1)})
  getFunctionalAnnotationChartFile(david, fileName = paste("/home/slancast/RNA_correlation_networks/david_results/", name,"AnnotationChart.txt",sep=""))
  getFunctionalAnnotationTableFile(david, fileName = paste("/home/slancast/RNA_correlation_networks/david_results/", name ,"AnnotationTable.txt",sep=""))
  getClusterReportFile(david, fileName = paste("/home/slancast/RNA_correlation_networks/david_results/", name ,"ClusterReport.txt",sep=""))
}

}