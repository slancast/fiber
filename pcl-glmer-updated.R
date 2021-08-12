library(icesTAF)
library("metaMA")
library(lme4)
library(ggplot2)
library(pbnm)

for (fiber_subset in c("Arabinoxylan","LCInulin","Mix")) {
  print(fiber_subset)
  
  load(paste("/home/slancast/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""))
  log_pcl <- as.matrix(t(log_pcl))
  class(log_pcl) <- "numeric"
  log_pcl <- log_pcl[rowVars(log_pcl) > 0,]
  
  metadata_rows <- nrow(pcl_metadata)
  log_pcl <- data.frame(cbind(t(pcl_metadata),t(log_pcl)))
  
  #This is only needed when combining terms in the mixed model
  #Because of the complication of our project throwing out dosages
  #Allows for fewer combinations of patient x dosage
  #log_pcl <- log_pcl[!log_pcl$Dose %in% "WashoutFinal",]
  
  #Only needed when directory is not present
  #library(icesTAF)
  #mkdir(paste("~/pclTernary",fiber_subset,sep=""))
  
  #From the metadata glmer program. This throws out the metadata variables that cannot
  #possibly be a part of the model, and those that are factors and can only be incorporated
  #as random variables.
  drops <- c("Fiber", "DOB.month.year", "date_time_blood_draw_clinicals", "Dose", "Participant","Alb.Creat.Interp")
  metadata_factors <- c("Metagenomic_seq_plate", "Sex", "Ethnicity", "IR_IS_from_SSPG","Alb.Creat.Interp")
  pcl_metadata <- pcl_metadata[!rownames(pcl_metadata) %in% drops,]
  pcl_metadata <- pcl_metadata[!rownames(pcl_metadata) %in% metadata_factors,]
  
  for (j in rownames(pcl_metadata)) {
    #There are two ways of determining which metadata to use
    #1 is to do it based on a priori medical knowledge of the event
    #My guess is this would include things like sequencing plate, kidney damage, etc.
    #The second would be to run over the metadata, and see which ones are
    #Associated with the most analytes, and then combine those
    #Perhaps a combination of these two approaches 
    mkdir(paste("~/pcl",fiber_subset,"/analyte-plots",j,sep=""))
    
    print(j)
    df = c(NA,NA,NA,NA)
    total_stats <- data.frame(intercept=double(),
                              metadata=double(),
                              total_var=double(),
                              analyte=double(),
                              dose_p=double(),
                              dose_p_unused=double(),
                              stringsAsFactors=FALSE)
    total_variance <- c()
    stats <- c()
    to_bind <- c()
    
    metadatum <- gsub("[^0-9\\.]", "", log_pcl[,j]) #removing any non-numeric characters from the data. 
    metadatum <- as.numeric(as.matrix(metadatum))
    
    for (i in ((metadata_rows+2):ncol(log_pcl))) { 
      #print(shapiro.test(as.numeric(log_pcl[,i]))) #testing for normality. Can be combined with lmer for increased statistical confidence
      print(colnames(log_pcl)[i])
      response <- as.numeric(as.matrix(log_pcl[,i]))
      response <- response * 1000000
      response <- round(response, digits=0)
      if (sum(response) <= 30) next 
      return_list <- tryCatch({ 
        model_nofiber <- glmer(response ~ 1 + (1|Participant) + (1|Metagenomic_seq_plate) + metadatum, data=log_pcl, family=poisson) #Poisson seems to work best at least for metagenome as the response variable; however only when there is a limited number of predictor variables. For the targeted program with several predictor variables lmer with the gaussian distribution is the best.
        model <- glmer(response ~ 1 + (1|Dose) + (1|Participant) + (1|Metagenomic_seq_plate) + metadatum, data=log_pcl, family=poisson)
        pbmmG2 <- pbnm(model,model_nofiber,nsim=100,tasks=10,cores=1,seed=1) #Using pbnm to assign a p-value to the ffect of dose 
        
        #Calculating varous statistics for the data. Statistics are labeled after the command
        variance <- as.data.frame(VarCorr(model)) ##Determining the variance composition of the random effect variables
        to_bind <- variance$vcov/sum(variance$vcov)*100 #Making the variance as a proportion of total.
        Vcov <- vcov(model, useScale = FALSE) #Covariance of the fixed effect variables
        betas <- fixef(model) #coefficients of the fixed effect variables
        se <- sqrt(diag(Vcov)) ####Beginning finding pvalue
        zval <- betas / se ###Continuing finding p-value
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) #Finding p-value
        
        #Adding pbnm pvalue to the stats, assuing it exists. If not, adding the confidence interval.
        if (is.na(summary(pbmmG2)$`P(>=obs)`)) {
          stats <- c(pval,sum(variance$vcov),colnames(log_pcl[i]),summary(pbmmG2)$`P(>=obs)`,summary(pbmmG2)$`#(>=obs)+noEst/runs`)
        } else {stats <- c(pval,sum(variance$vcov),colnames(log_pcl[i]),summary(pbmmG2)$`P(>=obs)`,"NA")
        plot <- ggplot(data=log_pcl, aes(x=log_pcl$Dose, y=response)) + #This seciton plots the data, so that we can better see the relationship, but only when a significant relationship is determined.
          geom_point() +
          scale_x_discrete(limits=c("Baseline","10","20","30","WashoutD3","WashoutD10","WashoutFinal"))
        ggsave(plot, file=paste("~/pcl",fiber_subset,"/analyte-plots",j,"/numeric",fiber_subset,j,colnames(log_pcl)[i],".pdf",sep=""))
        }
        
        #For matting thes "stats" object for the individual analyte
        stats <- data.frame(t(stats))
        colnames(stats) <- c("Intercept","metadata","total_var","analyte","dose_p","dose_p_unused")
        stats <- data.frame(stats)
        
        
        #Appending the stats to the data frames,
        to_bind <- c(to_bind)
        df <- rbind(df,to_bind)
        print("tot")
        total_stats <- rbind(total_stats,stats)
        
      },
      error = function(err){print(err)
        #Making sure that if there's an error the program will keep running
        to_bind = c(NA,NA,NA,NA)
        stats = c(NA,NA,NA,NA)
        })
    } #end data loop
    
    #Adding p-adjusted column, and polishing final data frame for output
    total_stats <- data.frame(total_stats)
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$metadata)), method = "BH", n = nrow(total_stats))
    total_stats$padj <- padjusted
    rownames(total_stats) <- make.names(total_stats[,4], unique=TRUE) #this is the analyte name, per design of the statistics matrix
    df <- df[-1,] #still don't like this, will replace with data frame asap
    df <- data.frame(df)
    rownames(df) <- make.names(total_stats[,4], unique=TRUE)
    total_stats <- subset(total_stats, select=-c(analyte))
    
    df <- df[,-3]
    colnames(df) <- c("Participant", "Week")
    
    write.table(df, file = paste("~/pcl",fiber_subset,"/numericCovar",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
    write.table(total_stats, file = paste("~/pcl",fiber_subset,"/numericStats",fiber_subset,j,"WeekPart.txt",sep=""), sep="\t")
    #df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/pclVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
    #df <- df[,-1] #Only to be used when reading a df from a table
    #df <- as.matrix(df) #Only to be used when reading a df from a table
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,j,"WeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$metadata)), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,j,"dose-p_WeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$dose_p_unused)), breaks=40)
    dev.off()

  } #Ending the metadata loop
  
} #Ending the fiber_subset loop



#This portion of the code is to make the ternary plots.
#However with the poisson distribution
#############################################################################
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#############################################################################


if (FALSE) {
#install.packages("ggtern")
  
  df2 <- na.omit(df)
  df2 <- as.matrix(df2)
  class(df2) <- "numeric"
  if  ( (nrow(df2) == 0) | (ncol(df2) <= 2) ) {
    next
  } 
  
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
ggsave(plot1, file=paste("~/pclTernary",fiber_subset,"/numeric",fiber_subset,j,"GeomPointWeekPart.pdf",sep=""))

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
ggsave(plot2, file=paste("~/pclTernary",fiber_subset,"/numeric",fiber_subset,j,"DensityWeekPart.pdf",sep=""))
}

