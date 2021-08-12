library(DESeq2)

resdata <- read.csv("/Users/SLancaster/Desktop/Projects/Fiber/RNAseq/Data/RNAseq_results_scg/RNAseq_results/fiber_week_Arabinoxylan.10_vs_Arabinoxylan.Baseline-diffexpr-resultsb.csv")
resdata <- data.frame(resdata)

png(paste("/Users/SLancaster/Desktop/diffexpr-volcanoplot.png",sep=""), 1200, 1000, pointsize=20)
volcanoplot(resdata,  sigthresh=0.1, lfcthresh=0.4, labelsig=FALSE, textcx=.8, xlim=c(-5, 5))
dev.off()

