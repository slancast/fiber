#!/usr/bin/R
#A tsne is something standard that we do for our data
#So this program will to try to attempt a standard version of doing that
#I think it should be pretty straigh forward, but I did if he had code but I haven't heard back
#I will plot them by participant 

load("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/RData/Normalized_Log_Arabinoxylan_metaphlan_df.RData")

require(Rtsne)

tsne <- Rtsne(normalized_metaphlan_df)
tsne_df <- data.frame(tsne$Y)

require(ggplot2)

p <- ggplot(tsne_df, aes(tsne_df[,1], tsne_df[,2]))
p + geom_point(aes(colour = factor(metaphlan_metadata$Participant)))


