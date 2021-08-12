#This project will perform supervised clustering

library(supclust)

fiber_subset = "Arabinoxylan"

metabolomics <- read.table(paste("/Users/SLancaster/Desktop/supervisied_clustering/metabolomics_",fiber_subset,".csv",sep=""),sep=",",header=TRUE)
metabolomics <- as.matrix(metabolomics)
class(metabolomics) <- "numeric"
x <- metabolomics[,2:ncol(metabolomics)]
y <- metabolomics[,1]

fit <- wilma(x, y, noc = 3, trace = 1)

summary(fit)
plot(fit)
fitted(fit)

fit <- pelora(x, y, noc = 3)

summary(fit)
plot(fit)

