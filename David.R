results_files <- list.files(path = "~/deseq_results_fall2019", pattern=".csv")
library(RDAVIDWebService)

for (i in results_files) { 
  print(i)
  name = substr(i,1,nchar(i)-22)
  print(name)
  result_file = read.csv(i, sep=",", row.names = 1, stringsAsFactors=FALSE, header=TRUE)
  background_genes <- result_file$Gene
  background_genes <- gsub("\\..*","",background_genes)
  sig_genes <- result_file[result_file$padj < 0.1,]$Gene
  sig_genes <- gsub("\\..*","",sig_genes)

  david<-DAVIDWebService(email="slancast@stanford.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setTimeOut(david, 1000000)
  setAnnotationCategories(david, c("GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"))
  background <- addList(david, background_genes, idType="ENSEMBL_GENE_ID",listName="Total_background_genes", listType="Background")

  result<- tryCatch({addList(david, sig_genes, idType="ENSEMBL_GENE_ID",listName= name, listType="Gene")},
                    error = function(err){c(1)})
  if (length(result) <= 1) next
  getFunctionalAnnotationChartFile(david, fileName = paste("/home/slancast/david_results_fall2019/", name,"AnnotationChart.txt",sep=""))
  getFunctionalAnnotationTableFile(david, fileName = paste("/home/slancast/david_results_fall2019/", name ,"AnnotationTable.txt",sep=""))
  getClusterReportFile(david, fileName = paste("/home/slancast/david_results_fall2019/", name ,"ClusterReport.txt",sep=""))
}

