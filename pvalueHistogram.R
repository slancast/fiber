#In this program I'll create a script to create a pvalue histogram
#Kevin made one of these for the exercise data, and I think it would be good for this project too
#Here is some code for an example pvalue historgram that I got from this webpage:
#https://github.com/andrewheiss/Attitudes-in-the-Arab-World/blob/master/figure12.R

if (FALSE) {
  #-----------------------------------------------------------
  # Figure 12: Grid plot for coefficient t-values (expanded)
  #-----------------------------------------------------------
  # Create dataframe of t-values for combined model and all individual models
  values <- data.frame(rbind(summary(all.countries.expanded)$coefficients[, "t value"],
                             t(sapply(seq(country.models.expanded), 
                                      FUN=function (x) summary(country.models.expanded[[x]])$coefficients[, "t value"]))
  ))
  
  # Clean up dataframe
  num.coefs <- length(attr(terms(form.expanded), "term.labels"))  # number of lhs coefficients
  values <- values[, seq(num.coefs)]  # Remove taus; only include actual coefficients
  values$country <- c("All countries", levels(barometer$country.name))  # Add column with country names
  
  # Build dataframe for plotting
  plot.data <- melt(values, id="country")  # Convert to long
  levels(plot.data$variable) <- nice.names.expanded  # Make coefficient names pretty
  plot.data$p.value <- pnorm(abs(plot.data$value), lower.tail=FALSE) * 2  # Calculate p-values for t-values
  plot.data$stars <- cut(plot.data$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
  plot.data$variable <- factor(plot.data$variable, levels=levels(with(plot.data[plot.data$country=="All countries",], reorder(variable, -value))))  # Sort coefficients by value
  plot.data$country <- factor(plot.data$country, levels=levels(with(plot.data[plot.data$variable==levels(plot.data$variable)[1],], reorder(country, -value))))  # Sort country by highest by most significant coefficient
  
  # Plot everything
  p <- ggplot(aes(x=country, y=variable, fill=value), data=plot.data)
  fig12 <- p + geom_tile() + scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") + 
    #   geom_text(aes(label=stars, color=value), size=8) + scale_colour_gradient(low="grey30", high="white", guide="none") +
    geom_text(aes(label=stars), color="black", size=5) + 
    labs(y=NULL, x=NULL, fill="t-value") + geom_vline(xintercept=1.5, size=1.5, color="grey50") + 
    theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
  fig12
}

#Especially without the RNA data complete, I don't think I can get a final dataset,
#but it might be with some of the DESeq data we alraedy have

##This data will match up the genes
resdata1 = read.table("/Volumes/My_Book_Duo/RNAseq_results_scg/RNAseq_results/fiber_week_Arabinoxylan.Baseline_vs_Arabinoxylan.20-diffexpr-resultsb.csv", sep=",", header=TRUE,row.names=1)
resdata2 = read.table("/Volumes/My_Book_Duo/RNAseq_results_scg/RNAseq_results/fiber_week_Arabinoxylan.Baseline_vs_Arabinoxylan.30-diffexpr-resultsb.csv", sep=",", header=TRUE,row.names=1)
Summary_Data <- cbind(Gene1 = as.character(resdata1$Gene), pvalue1 = resdata1$padj)
Summary_Data <- data.frame(Summary_Data)
Gene2 <- resdata2$Gene[match(resdata1$Gene, resdata2$Gene)]
pvalue2 <- resdata2$padj[match(resdata1$Gene, resdata2$Gene)]
Summary_Data2 <- data.frame(cbind(Gene2 = as.character(Gene2), pvalue2 = pvalue1))
Summary_Data <- data.frame(cbind(Summary_Data, pvalue2))
rownames(Summary_Data) <- Summary_Data$Gene1
Summary_Data <- data.frame(Summary_Data)

library(reshape2)
Summary_Data.melted <- melt(Summary_Data)
p <- ggplot(Summary_Data.melted, aes(y = Gene, x = variable, fill = value)) + geom_tile()

###This portion is to use the heatmap function rather than ggplot
if (FALSE) {
Summary_Data <- Summary_Data[,-1]
Summary_Data <- as.matrix(Summary_Data)
class(Summary_Data) <- "numeric"
heatmap(Summary_Data)
Summary_Data <- data.frame(cbind(Gene = rownames(Summary_Data), Summary_Data))
Summary_Data <- data.frame(cbind(Gene = rownames(Summary_Data), Summary_Data))
}

