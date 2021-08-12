#Now that I've created the master metadata file, I want to facilitate adding it to the data files
#This program won't be perfect, but is designed to facilitate
#The column names of the metadata and the data files must match.

data_file <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/MasterRNASeq_DeseqNormalized_w_Metadata_Ensembl.txt", sep="\t", row.names = 1, header=TRUE, stringsAsFactors=FALSE)

master_metadata <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Data/metadata_master.txt", sep="\t", row.names = 1, header=TRUE, stringsAsFactors=FALSE)

#Making data files conform to metadata names structure
colnames(data_file) <- gsub("Washout.D", "WashoutD",colnames(data_file))
colnames(data_file) <- gsub("lcInulin", "LCInulin",colnames(data_file))
colnames(data_file) <- gsub("X69", "X069",colnames(data_file))
colnames(data_file)[which(!colnames(data_file) %in% colnames(master_metadata))]

data_file_with_meta <- data.frame(matrix(nrow = 0, ncol = length(colnames(data_file))))

metadata_to_add <- data.frame(matrix(nrow = nrow(master_metadata), ncol = 0))

empty <- rep(NA, nrow(master_metadata))
for (i in colnames(data_file)) {
  print(i)
  if (is.na(match(i, colnames(master_metadata)))){
    metadata_to_add <- cbind(metadata_to_add,empty)
  } else { 
  print(match(i, colnames(master_metadata)))
  metadata_to_add <- cbind(metadata_to_add,master_metadata[,match(i, colnames(master_metadata))])
    }
}

rownames(metadata_to_add) <- rownames(master_metadata)
data_file_with_meta <- rbind(data_file_with_meta, metadata_to_add)
colnames(data_file_with_meta) <- colnames(data_file)

output <- rbind(data_file_with_meta, data_file)

write.table(data.frame("run"=rownames(output),output), file = "/Users/SLancaster/Desktop/data_with_metadata.txt", sep="\t", row.names=FALSE)

