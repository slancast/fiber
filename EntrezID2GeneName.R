
entrezid <- function( resdata ) {
  require(EnsDb.Hsapiens.v79)
  a = resdata$Gene #the column to iterate over will be different if I'm using res vs resdata
  tmp=gsub("\\..*","",a)
  tmp <- as.character(tmp)
  txdb <- EnsDb.Hsapiens.v79 
  df <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "ENTREZID")
  df2 <- AnnotationDbi::select(txdb, keys = tmp, keytype = "GENEID", columns = "SYMBOL")
  
  ENTREZID <- c() 
  SYMBOL <- c() 
  counter1 <- 0
  for (i in tmp) { 
    counter1 <- counter1 + 1
    j <- match(i,df$GENEID)
    ENTREZID <- c(ENTREZID, toString(df[j,][2]))
    SYMBOL <- c(SYMBOL, toString(df2[j,][2]))}
  resdata$ENTREZID <- ENTREZID
  resdata$SYMBOL <- SYMBOL
  
  resdata
}
