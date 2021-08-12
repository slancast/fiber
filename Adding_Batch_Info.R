#I have to add some of the batch info into the RNA data
#It seems the most efficient way to do this will be using a code.

batch_info <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/RNABatchInfo.txt",sep="\t",header=TRUE)

rna_data <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/RNAseq_data_with_metadata.txt",sep="\t",header=TRUE,row.names=1)

colnames(rna_data) <- gsub("X","",colnames(rna_data))

batch_info <- t(batch_info)
colnames(batch_info) <- batch_info[1,]
batch_info <- data.frame(batch_info)
colnames(batch_info) <- gsub("X","",colnames(batch_info))
batch_info <- batch_info[-1,]

colnames(batch_info) <- gsub("Washout.D", "WashoutD",colnames(batch_info))
colnames(batch_info) <- gsub("lcInulin", "LCInulin",colnames(batch_info))


match(colnames(rna_data), colnames(batch_info))

batch <- c()
for (i in colnames(rna_data)) {
  print(i)
  if (is.na(match(i, colnames(batch_info)))){
    
    batch <- c(batch,NA)
  } else { 
    print(batch_info[1,match(i, colnames(batch_info))])
    batch <- c(batch,as.numeric(as.matrix(batch_info[1,match(i, colnames(batch_info))])))
  }
}

batch <- t(data.frame(batch))

colnames(batch) <- colnames(rna_data)
rna_data2 <- rbind(batch, rna_data)

write.table(data.frame("run"=rownames(rna_data2),rna_data2), file = "/Users/SLancaster/Desktop/data_with_metadata.txt", sep="\t", row.names=FALSE)


