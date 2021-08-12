library(icesTAF)
library("metaMA")
library(lme4)
library(ggplot2)
library(pbnm)
library(RLRsim)

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
  
  #Dropping metadata that cannot possibly be a part of the model
  #Also since in this case the metadatum is the response variable it cannot be a factor
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
                              response_var=double(),
                              total_var=double(),
                              analyte=double(),
                              dose_p=double(),
                              dose_p_unused=double(),
                              stringsAsFactors=FALSE)
    
    total_variance <- c()
    
    #Here I will create a list of the metadata that should be factors vs those that should be numeric
    metadatum <- gsub("[^0-9\\.]", "", log_pcl[,j]) #removing any non-numeric characters from the data. 
    metadatum <- as.numeric(as.matrix(metadatum))
    metadatum <- metadatum * 10
    metadatum <- round(metadatum, digits=0)
    
    for (i in (metadata_rows + 2):ncol(log_pcl)) {
      #print(shapiro.test(as.numeric(log_pcl[,i])))
      #The response varaible (before the ~) with a poisson
      #distribution 
      response <- as.numeric(as.matrix(log_pcl[,i]))
      if (sum(response) <= 30/100000) next 
      plot = NA
      print(colnames(log_pcl)[i])
      return_list <- tryCatch({ 
        #these are the values that make the most Arabinoxylan values significant
        #although that was with them as factors, and not numeric. I will rerun everything with them as numeric.
        #the ones labeled with "numeric" are the ones run as numeric, the unlabeled ones are as factors
        model_nofiber <- lmer(metadatum ~ 1 + (1|Participant) + (1|Metagenomic_seq_plate) + response, data=log_pcl)
        model <- lmer(metadatum ~ 1 + (1|Dose) + (1|Participant) + (1|Metagenomic_seq_plate) + response, data=log_pcl)
        pbmmG2 <- pbnm(model,model_nofiber,nsim=100,tasks=10,cores=1,seed=1) #Using pbnm to assign a p-value to the ffect of dose 
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
          plot <- ggplot(data=log_pcl, aes(x=log_pcl$Dose, y=response)) + #This seciton plots the data, so that we can better see the relationship.
            geom_point() +
            scale_x_discrete(limits=c("Baseline","10","20","30","WashoutD3","WashoutD10","WashoutFinal"))
          ggsave(plot, file=paste("~/pcl",fiber_subset,"/analyte-plots",j,"/numeric",fiber_subset,j,colnames(log_pcl)[i],".pdf",sep=""))
          }
          
          
        stats <- data.frame(t(stats))
        colnames(stats) <- c("Intercept","response_var","total_var","analyte","dose_p","log_lik")
        stats <- data.frame(stats)
        
        to_bind <- c(to_bind)
        df <- rbind(df,to_bind)
        print("tot")
        total_stats <- rbind(total_stats,stats)
        
      },
      error = function(err){print(err)
        to_bind = c(NA,NA,NA,NA)
        stats = c(NA,NA,NA,NA)})
      
    } #end data loop
    
    total_stats <- data.frame(total_stats)
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$response_var)), method = "BH", n = nrow(total_stats))
    total_stats$padj <- padjusted
    rownames(total_stats) <- make.names(total_stats[,4], unique=TRUE) #with the random effects interaction term
    df <- df[-1,]
    df <- data.frame(df)
    rownames(df) <- make.names(total_stats[,4], unique=TRUE)
    total_stats <- subset(total_stats, select=-c(analyte))
    
    df <- df[,-3]
    colnames(df) <- c("Participant", "Week")
    
    write.table(df, file = paste("~/pcl",fiber_subset,"metadata-response/numericCovar",fiber_subset,"WeekPart.txt",sep=""), sep="\t")
    write.table(total_stats, file = paste("~/pcl",fiber_subset,"metadata-response/numericStats",fiber_subset,"WeekPart.txt",sep=""), sep="\t")
    #df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/pclVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
    #df <- df[,-1] #Only to be used when reading a df from a table
    #df <- as.matrix(df) #Only to be used when reading a df from a table
    pdf(paste("~/pcl",fiber_subset,"metadata-response/numericHist",fiber_subset,j,"WeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$response_var)), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"metadata-response/numericHist",fiber_subset,j,"dose-p_WeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$log_lik)), breaks=40)
    dev.off()
    
  } #Ending the metadata loop
  
} #Ending the fiber_subset loop


