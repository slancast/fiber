
library("metaMA")
library(fdrtool)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  print(fiber_subset)
  
  load(paste("/home/slancast/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""))
  log_pcl <- as.matrix(t(log_pcl))
  class(log_pcl) <- "numeric"
  log_pcl <- log_pcl[rowVars(log_pcl) > 0,]
  
  log_pcl <- data.frame(cbind(t(pcl_metadata),t(log_pcl)))
  
  log_pcl <- log_pcl[!log_pcl$Dose %in% "WashoutFinal",]
  
  library(icesTAF)
  mkdir(paste("~/pclTernary",fiber_subset,sep=""))
  
  for (j in rownames(pcl_metadata)) {
    print(j)
    df = c(NA,NA,NA,NA)
    total_stats = c(NA,NA,NA,NA)
    total_variance <- c()
    library(lme4)
    for (i in nrow(pcl_metadata):ncol(log_pcl)) { 
      #print(shapiro.test(as.numeric(log_pcl[,i])))
      
      return_list <- tryCatch({ 
        model <- lmer(as.numeric(log_pcl[,i]) ~ 1 + (1|Participant) + (1|Dose) + (1|Metagenomic_seq_plate) +  (1|log_pcl$Alb.Creat.Ratio) + (1|log_pcl$eGFR) + as.numeric(log_pcl[,j]), data=log_pcl, singular.ok = TRUE) #
        #model <- lmer(as.numeric(log_pcl[,i]) ~ 1 + (1|Dose) + (1|Participant) + (1|as.numeric(log_pcl[,j])), data=log_pcl)
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
        stats <- data.frame(stats)
        
        return_list <- list("to_bind" = to_bind, "stats" = stats[2,])
      },
      error = function(err){print(err)
        return_list <- list("to_bind" = c(NA,NA,NA,NA), "stats" = c(NA,NA,NA,NA,NA))
        return(return_list)})
      to_bind <- return_list$to_bind
      stats <- return_list$stats
      to_bind <- c(to_bind,colnames(log_pcl[i]))
      df <- rbind(df,to_bind)
      total_stats <- rbind(total_stats,stats)
    } #end data loop
    
    stats <- data.frame(stats)
    padjusted <- p.adjust(total_stats$pval, method = "BH", n = nrow(total_stats))
    total_stats$padj <- padjusted
    total_stats <- total_stats[-1,]
    df <- df[-1,]
    rownames(total_stats) <- df[,4] #with the random effects interaction term
    rownames(df) <- df[,4]
    df <- df[-1,]
    df <- df[,-4]
    df <- data.frame(df)
    colnames(df) <- c("Participant", "Week", "Metagenomic_seq_plate", "Other")
    
    write.table(df, file = paste("~/pcl",fiber_subset,"/pclREMPCovar",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
    write.table(total_stats, file = paste("~/pcl",fiber_subset,"/pclREMPStats",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
    #df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/pclVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
    #df <- df[,-1] #Only to be used when reading a df from a table
    #df <- as.matrix(df) #Only to be used when reading a df from a table
    df2 <- na.omit(df)
    df2 <- as.matrix(df2)
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
      labs(  title= "pcl",
             fill = "Value, V",alpha="|V - 0|")
    ggsave(plot1, file=paste("~/pclTernary",fiber_subset,"/pclREMP",fiber_subset,j,"GeomPointWeekPart.pdf",sep=""))
    
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
      labs(  title= "pcl",
             fill = "Value, V",alpha="|V - 0|")
    ggsave(plot2, file=paste("~/pclTernary",fiber_subset,"/pclREMP",fiber_subset,j,"DensityWeekPart.pdf",sep=""))
    
  } #Ending the metadata loop
  
} #Ending the fiber_subset loop



