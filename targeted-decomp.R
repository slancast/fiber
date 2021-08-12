library(icesTAF)
library("metaMA")
library(lme4)
library(ggplot2)
library(pbnm)

for (fiber_subset in c("Arabinoxylan", "LCInulin", "Mix")) { 

  print(fiber_subset)
  
  load(paste("/home/slancast/Full_Log_",fiber_subset,"_pcl_df.RData",sep=""))
  log_pcl <- as.matrix(t(log_pcl))
  class(log_pcl) <- "numeric"
  log_pcl <- log_pcl[rowVars(log_pcl) > 0,]
  
  log_pcl <- data.frame(cbind(t(pcl_metadata),t(log_pcl)))
  
  #This is only needed when combining terms in the mixed model
  #Because of the complication of our project throwing out dosages
  #Allows for fewer combinations of patient x dosage
  #log_pcl <- log_pcl[!log_pcl$Dose %in% "WashoutFinal",]
  
  #Only needed when directory is not present

    #There are two ways of determining which metadata to use
    #1 is to do it based on a priori medical knowledge of the event
    #My guess is this would include things like sequencing plate, kidney damage, etc.
    #The second would be to run over the metadata, and see which ones are
    #Associated with the most analytes, and then combine those
    #Perhaps a combination of these two approaches 
    mkdir(paste("~/pcl",fiber_subset,sep=""))
    
    total_variance <- c()
    eGFR <- gsub("[^0-9\\.]", "", log_pcl$eGFR) #removing any non-numeric characters from the data. 
    eGFR <- as.numeric(as.matrix(eGFR))
    eGFR <- eGFR
    eGFR_model <- round(eGFR, digits=0)
    
    Hematocrit <- gsub("[^0-9\\.]", "", log_pcl$Hematocrit) #removing any non-numeric characters from the data. 
    Hematocrit <- as.numeric(as.matrix(Hematocrit))
    Hematocrit <- Hematocrit * 10
    Hematocrit_model <- round(Hematocrit, digits=0)
    
    High.Sensitivity.CRP <- gsub("[^0-9\\.]", "", log_pcl$High.Sensitivity.CRP) #removing any non-numeric characters from the data. 
    High.Sensitivity.CRP <- as.numeric(as.matrix(High.Sensitivity.CRP))
    High.Sensitivity.CRP <- High.Sensitivity.CRP * 10
    High.Sensitivity.CRP_model <- round(High.Sensitivity.CRP, digits=0)
    
    LDL.HDL.Ratio <- gsub("[^0-9\\.]", "", log_pcl$LDL.HDL.Ratio) #removing any non-numeric characters from the data. 
    LDL.HDL.Ratio <- as.numeric(as.matrix(LDL.HDL.Ratio))
    LDL.HDL.Ratio <- LDL.HDL.Ratio * 10
    LDL.HDL.Ratio_model <- round(LDL.HDL.Ratio, digits=0)
    
    Hemoglobin.A1c <- gsub("[^0-9\\.]", "", log_pcl$Hemoglobin.A1c) #removing any non-numeric characters from the data. 
    Hemoglobin.A1c <- as.numeric(as.matrix(Hemoglobin.A1c))
    Hemoglobin.A1c <- Hemoglobin.A1c * 10
    Hemoglobin.A1c_model <- round(Hemoglobin.A1c, digits=0)
    
    RBC <- gsub("[^0-9\\.]", "", log_pcl$RBC) #removing any non-numeric characters from the data. 
    RBC <- as.numeric(as.matrix(RBC))
    RBC <- RBC * 10
    RBC_model <- round(RBC, digits=0)

    Triglyceride..Ser.Plas <- gsub("[^0-9\\.]", "", log_pcl$Triglyceride..Ser.Plas) #removing any non-numeric characters from the data. 
    Triglyceride..Ser.Plas <- as.numeric(as.matrix(Triglyceride..Ser.Plas))
    Triglyceride..Ser.Plas <- Triglyceride..Ser.Plas * 10
    Triglyceride..Ser.Plas_model <- round(Triglyceride..Ser.Plas, digits=0)
    
    df <- data.frame(Participant=double(),
                    Dose = double(),
                    Metagenomic_seq_plate = double(),
                    Other = double())
    
    total_stats <- data.frame(intercept=double(),
                              eGFR=double(),
                              Hematocrit=double(),
                              High.Sensitivity.CRP=double(),
                              LDL.HDL.Ratio=double(),
                              RBC=double(),
                              Hemoglobin.A1c=double(),
                              Triglyceride..Ser.Plas=double(),
                              total_var=double(),
                              analyte=double(),
                              dose_p=double(),
                              dose_p_unused=double(),
                              intercept_coef=double(),
                              eGFR_coef=double(),
                              Hematocrit_coef=double(),
                              High.Sensitivity.CRP_coef=double(),
                              LDL.HDL.Ratio_coef=double(),
                              RBC_coef=double(),
                              Hemoglobin.A1c_coef=double(),
                              Triglyceride..Ser.Plas_coef=double(),
                              stringsAsFactors=FALSE)
    column_names <- colnames(total_stats)
    
    stats <- c()
    to_bind <- c()
    counter <- 0
    
    for (i in (nrow(pcl_metadata)+2):ncol(log_pcl)) { 
      #print(shapiro.test(as.numeric(log_pcl[,i])))
      print(colnames(log_pcl)[i])
      response <- as.numeric(as.matrix(log_pcl[,i]))
      response <- response * 1000000
      response <- round(response, digits=0)
      if (sum(response) <= 30) next 
      return_list <- tryCatch({ 
        model_nofiber <- lmer(response ~ 1 + (1|Participant) + (1|Metagenomic_seq_plate) + eGFR_model + Hematocrit_model + High.Sensitivity.CRP_model + LDL.HDL.Ratio_model + RBC_model + Hemoglobin.A1c_model + Triglyceride..Ser.Plas_model, data=log_pcl) #, family=poisson(link=identity) Poisson seems to work best at least for metagenome as the response variable
        model <- lmer(response ~ 1 + (1|Dose) + (1|Participant) + (1|Metagenomic_seq_plate) + eGFR_model + Hematocrit_model + High.Sensitivity.CRP_model + LDL.HDL.Ratio_model + RBC_model + Hemoglobin.A1c_model + Triglyceride..Ser.Plas_model, data=log_pcl) #, family=poisson
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
        }
        
        stats <- c(stats, coef(summary(model))[,1])
        #For matting thes "stats" object for the individual analyte
        stats <- data.frame(t(stats))
        colnames(stats) <- column_names
        stats <- data.frame(stats)
        
      }, #End of the trycatch first clause
      error = function(err){print(err)
        #Making sure that if there's an error the program will keep running
        #Not sure why but this allows the program to keep running. If it ain't broke don't fix it.
        to_bind = c()
        stats = c()
        }) #End of the total trycatch
      
      #Appending the stats to the data frames,
      to_bind <- c(to_bind)
      df <- rbind(df,to_bind)
      counter = counter + 1
      print(counter)
      total_stats <- rbind(total_stats,stats)
    } #end data loop
    
    #Adding p-adjusted column, and polishing final data frame for output
    total_stats <- data.frame(total_stats)
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$eGFR)), method = "BH", n = nrow(total_stats))
    total_stats$eGFR_adj <- padjusted
    
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Hematocrit)), method = "BH", n = nrow(total_stats))
    total_stats$Hematocrit_adj <- padjusted
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$High.Sensitivity.CRP)), method = "BH", n = nrow(total_stats))
    total_stats$HHigh.Sensitivity.CRP_adj <- padjusted
    
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$LDL.HDL.Ratio)), method = "BH", n = nrow(total_stats))
    total_stats$LDL.HDL.Ratio_adj <- padjusted
    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Hemoglobin.A1c)), method = "BH", n = nrow(total_stats))
    total_stats$Hemoglobin.A1c_adj <- padjusted

    padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Triglyceride..Ser.Plas)), method = "BH", n = nrow(total_stats))
    total_stats$Triglyceride..Ser.Plas_adj <- padjusted
    
    rownames(total_stats) <- make.names(total_stats$analyte, unique=TRUE) #with the random effects interaction term
    df <- data.frame(df)
    rownames(df) <- make.names(total_stats$analyte, unique=TRUE)
    total_stats <- subset(total_stats, select=-c(analyte))
    
    df <- df[,-5] #This is used because I've appended the names to the dataframe
    
    write.table(df, file = paste("~/pcl",fiber_subset,"/numericCovar",fiber_subset,"WeekPart.txt",sep=""), sep="\t")
    write.table(data.frame("pathway"=rownames(total_stats),total_stats), file = paste("~/pcl",fiber_subset,"/numericStats",fiber_subset,"WeekPart.txt",sep=""), row.names=FALSE, sep="\t")
    #df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/pclVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
    #df <- df[,-1] #Only to be used when reading a df from a table
    #df <- as.matrix(df) #Only to be used when reading a df from a table
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"eGFRWeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$eGFR), main="eGFR"), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"HematocritWeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$Hematocrit), main="Hematocrit", xlab="p-value"), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"High.Sensitivity.CRP_modelWeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$High.Sensitivity.CRP), main="CRP", xlab="p-value"), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"LDL.HDL.Ratio_model.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$LDL.HDL.Ratio), main="LDL", xlab="p-value"), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"Hemoglobin.A1c_modelWeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$Hemoglobin.A1c), main="A1c", xlab="p-value"), breaks=40)
    dev.off()
 
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"Triglyceride..Ser.Plas_modelWeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$Triglyceride..Ser.Plas), main="Triglyceride", xlab="p-value"), breaks=40)
    dev.off()
    
    pdf(paste("~/pcl",fiber_subset,"/numericHist",fiber_subset,"dose-p_WeekPart.pdf",sep=""))
    hist(as.numeric(as.matrix(total_stats$dose_p_unused), main="dose_p", xlab="p-value"), breaks=40)
    dev.off()

    alarm()
    
} #Ending the fiber_subset loop

