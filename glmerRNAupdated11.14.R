
rna_df = read.csv("~/DeseqNormRNAwithMetadata.txt", sep="\t", row.names = 1, stringsAsFactors=FALSE)

rna_metadata_rows = 64
rna_df <- data.frame(lapply(rna_df, as.character), stringsAsFactors=FALSE, row.names = rownames(rna_df)) #This needs to be done for the transcript data, but causes problems with the pcl data.
rna_df <- data.frame(t(rna_df))
fiber_subset = "Mix"
rna_df <- rna_df[which(rna_df$fiber==fiber_subset),] #subsetting by fiber
rna_df <- data.frame(lapply(rna_df, as.factor), stringsAsFactors=FALSE, row.names = rownames(rna_df))
rna_df <- data.frame(t(rna_df))
rna_metadata = head(rna_df,rna_metadata_rows)
rna_df <- data.frame(t(rna_df))

log_rna <- cbind(t(rna_metadata),rna_df)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  print(fiber_subset)
  
  #load(paste("/home/slancast/Full_Log_",fiber_subset,"_rna_df.RData",sep=""))
  
  library(icesTAF)
  mkdir(paste("~/rnaTernary",fiber_subset,sep=""))
  
  for (j in colnames(rna_df)) {
    print(j)
    df = c(NA,NA,NA,NA)
    total_stats = c(NA,NA,NA,NA)
    total_variance <- c()
    library(lme4)
    for (i in nrow(rna_metadata):ncol(log_rna)) { 
      #print(shapiro.test(as.numeric(log_rna[,i])))
      
      to_bind <- tryCatch({ 
        model <- glmer(as.numeric(log_rna[,i]) ~ 1 + (1|week) + (1|participant) + as.numeric(log_rna[,j]), data=log_rna)
        variance <- as.data.frame(VarCorr(model)) ##Determining the variance composition of the random effect variables
        to_bind <- variance$vcov/sum(variance$vcov)*100 #Making the variance as a proportion of total.
        Vcov <- vcov(model, useScale = FALSE) #Covariance of the fixed effect variables
        betas <- fixef(model) #coefficients of the fixed effect variables
        se <- sqrt(diag(Vcov)) ####Beginning finding pvalue
        zval <- betas / se ###Continuing finding p-value
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) #Finding p-value
        stats <- cbind(betas, se, zval, pval,sum(variance$vcov))
        colnames(stats)[5] <- "tot_var" 
        rownames(stats) <- c("Intercept", j)
        return_list <- list("to_bind" = to_bind, "stats" = stats[2,])
      },
      error = function(err){print(err)
        return(c(NA,NA,NA))})
      to_bind <- return_list$to_bind
      stats <- return_list$stats
      to_bind <- c(to_bind,colnames(log_rna[i]))
      df <- rbind(df,to_bind)
      total_stats <- rbind(total_stats,stats)
    } 
    #
    rownames(total_stats) <- df[,4]
    total_stats <- total_stats[-1,]
    rownames(df) <- df[,4]
    df <- df[-1,]
    df <- df[,-4]
    colnames(df) <- c("Participant", "Week", "Other")
    
    write.table(df, file = paste("~/rnaTernary",fiber_subset,"/rnaCovar",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
    write.table(total_stats, file = paste("~/rnaTernary",fiber_subset,"/rnaStats",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
    #df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/rnaVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
    #df <- df[,-1] #Only to be used when reading a df from a table
    #df <- as.matrix(df) #Only to be used when reading a df from a table
    df2 <- na.omit(df)
    class(df2) <- "numeric"
    if  ( (nrow(df2) == 0) | (ncol(df2) <= 2) ) {
      next
    } 
    #install.packages("ggtern")
    df2 <- data.frame(df2)
    library(ggtern)
    plot1 <- ggtern(data = df2,aes(x = Week ,y = Participant, z = Other)) +
      #geom_density_tern(aes(color=..level..), bins=4500) +
      #geom_density_tern(aes(fill=..level.., alpha=abs(..level..)),bins=500, binwidth=100) +
      geom_point(size=0.5, alpha=0.1, color=c("red")) +
      #stat_density_tern(geom="polygon", n=4000, bins=500, aes(fill=..level..)) +
      scale_fill_gradient(low="yellow",high="red")   + 
      scale_color_gradient(low="yellow",high="red")   +
      theme_bw()  +
      theme(legend.justification=c(0,1), legend.position=c(0,1)) + #, panel.background = element_rect(fill = "light gray", colour = "light gray")
      theme_nogrid() +
      guides(fill = guide_colorbar(order=1),
             alpha= guide_legend(order=2),
             color="none") + 
      labs(  title= "rna",
             fill = "Value, V",alpha="|V - 0|")
    ggsave(plot1, file=paste("~/rnaTernary",fiber_subset,"/rna",fiber_subset,j,"GeomPointWeekPart.pdf",sep=""))
    
    plot2 <- ggtern(data = df2,aes(x = Week ,y = Participant, z = Other)) +
      stat_density_tern(geom="polygon", n=4000, bins=500, aes(fill=..level..)) +
      scale_fill_gradient(low="yellow",high="red")   + 
      scale_color_gradient(low="yellow",high="red")   +
      theme_bw()  +
      theme(legend.justification=c(0,1), legend.position=c(0,1)) + #, panel.background = element_rect(fill = "light gray", colour = "light gray")
      theme_nogrid() +
      guides(fill = guide_colorbar(order=1),
             alpha= guide_legend(order=2),
             color="none") + 
      labs(  title= "rna",
             fill = "Value, V",alpha="|V - 0|")
    ggsave(plot2, file=paste("~/rnaTernary",fiber_subset,"/rna",fiber_subset,j,"DensityWeekPart.pdf",sep=""))
    
  } #Ending the metadata loop
  
} #Ending the fiber_subset loop



