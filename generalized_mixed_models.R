#This will be the same as the targeted decomp except using the fiber as a 
#quantitatve variable not a categorical variable with random effects

library(icesTAF)
library("metaMA")
library(lme4)
library(ggplot2)
library(pbnm)

for (fiber_subset in c("LCInulin","Mix")) { #"Arabinoxylan",
  print(fiber_subset) 
  
  load(paste("/home/slancast/Full_Log_",fiber_subset,"_metaphlan_df.RData",sep=""))
  log_metaphlan <- as.matrix(t(log_metaphlan))
  class(log_metaphlan) <- "numeric"
  log_metaphlan <- log_metaphlan[rowVars(log_metaphlan) > 0,]
  
  log_metaphlan <- data.frame(cbind(t(metaphlan_metadata),t(log_metaphlan)))
  
  #log_metaphlan <- log_metaphlan[!log_metaphlan$Dose %in% "WashoutD3",]
  log_metaphlan$Dose <- gsub("WashoutD3", 0, log_metaphlan$Dose)
  log_metaphlan$Dose <- gsub("Baseline", 0, log_metaphlan$Dose)
  log_metaphlan$Dose <- gsub("WashoutD10", 0, log_metaphlan$Dose)
  log_metaphlan$Dose <- gsub("WashoutFinal", 0, log_metaphlan$Dose)
  
  #Only needed when directory is not present
  #library(icesTAF)
  #mkdir(paste("~/metaphlanTernary",fiber_subset,sep=""))
  
  #There are two ways of determining which metadata to use
  #1 is to do it based on a priori medical knowledge of the event
  #My guess is this would include things like sequencing plate, kidney damage, etc.
  #The second would be to run over the metadata, and see which ones are
  #Associated with the most analytes, and then combine those
  #Perhaps a combination of these two approaches 
  mkdir(paste("~/metaphlan",fiber_subset,"-numericdose/",sep=""))
  
  total_variance <- c()
  
  Dose <- gsub("[^0-9\\.]", "", log_metaphlan$Dose) #removing any non-numeric characters from the data. 
  Dose <- as.numeric(as.matrix(Dose))
  Dose <- Dose
  Dose_model <- round(Dose, digits=0)
  
  eGFR <- gsub("[^0-9\\.]", "", log_metaphlan$eGFR) #removing any non-numeric characters from the data. 
  eGFR <- as.numeric(as.matrix(eGFR))
  eGFR <- eGFR
  eGFR_model <- round(eGFR, digits=0)
  
  Hematocrit <- gsub("[^0-9\\.]", "", log_metaphlan$Hematocrit) #removing any non-numeric characters from the data. 
  Hematocrit <- as.numeric(as.matrix(Hematocrit))
  Hematocrit <- Hematocrit * 10
  Hematocrit_model <- round(Hematocrit, digits=0)
  
  High.Sensitivity.CRP <- gsub("[^0-9\\.]", "", log_metaphlan$High.Sensitivity.CRP) #removing any non-numeric characters from the data. 
  High.Sensitivity.CRP <- as.numeric(as.matrix(High.Sensitivity.CRP))
  High.Sensitivity.CRP <- High.Sensitivity.CRP * 10
  High.Sensitivity.CRP_model <- round(High.Sensitivity.CRP, digits=0)
  
  LDL..Calculated. <- gsub("[^0-9\\.]", "", log_metaphlan$LDL..Calculated.) #removing any non-numeric characters from the data. 
  LDL..Calculated. <- as.numeric(as.matrix(LDL..Calculated.))
  LDL..Calculated. <- LDL..Calculated. * 10
  LDL..Calculated._model <- round(LDL..Calculated., digits=0)
  
  Hemoglobin.A1c <- gsub("[^0-9\\.]", "", log_metaphlan$Hemoglobin.A1c) #removing any non-numeric characters from the data. 
  Hemoglobin.A1c <- as.numeric(as.matrix(Hemoglobin.A1c))
  Hemoglobin.A1c <- Hemoglobin.A1c * 10
  Hemoglobin.A1c_model <- round(Hemoglobin.A1c, digits=0)
  
  RBC <- gsub("[^0-9\\.]", "", log_metaphlan$RBC) #removing any non-numeric characters from the data. 
  RBC <- as.numeric(as.matrix(RBC))
  RBC <- RBC * 10
  RBC_model <- round(RBC, digits=0)
  
  Glucose..Ser.Plas <- gsub("[^0-9\\.]", "", log_metaphlan$Glucose..Ser.Plas) #removing any non-numeric characters from the data. 
  Glucose..Ser.Plas <- as.numeric(as.matrix(Glucose..Ser.Plas))
  Glucose..Ser.Plas <- Glucose..Ser.Plas * 10
  Glucose..Ser.Plas_model <- round(Glucose..Ser.Plas, digits=0)
  
  Triglyceride..Ser.Plas <- gsub("[^0-9\\.]", "", log_metaphlan$Triglyceride..Ser.Plas) #removing any non-numeric characters from the data. 
  Triglyceride..Ser.Plas <- as.numeric(as.matrix(Triglyceride..Ser.Plas))
  Triglyceride..Ser.Plas <- Triglyceride..Ser.Plas * 10
  Triglyceride..Ser.Plas_model <- round(Triglyceride..Ser.Plas, digits=0)
  
  df = c(NA,NA,NA,NA)
  total_stats <- data.frame(intercept=double(),
                            eGFR=double(),
                            Dose=double(),
                            Hematocrit=double(),
                            High.Sensitivity.CRP=double(),
                            LDL.HDL.Ratio=double(),
                            RBC=double(),
                            Hemoglobin.A1c=double(),
                            Glucose..Ser.Plas=double(),
                            Triglyceride..Ser.Plas=double(),
                            #total_var=double(),
                            analyte=double(),
                            intercept_coef=double(),
                            eGFR_coef=double(),
                            Dose_coef=double(),
                            Hematocrit_coef=double(),
                            High.Sensitivity.CRP_coef=double(),
                            LDL.HDL.Ratio_coef=double(),
                            RBC_coef=double(),
                            Hemoglobin.A1c_coef=double(),
                            Glucose..Ser.Plas_coef=double(),
                            Triglyceride..Ser.Plas_coef=double(),
                            #prtic_p=double(),
                            #prtic_p_unused=double(),
                            stringsAsFactors=FALSE)
  column_names <- colnames(total_stats)
  
  stats <- c()
  to_bind <- c()
  counter <- 0
  
  for (i in (nrow(metaphlan_metadata)+1):ncol(log_metaphlan)) { 
    #print(shapiro.test(as.numeric(log_metaphlan[,i])))
    print(colnames(log_metaphlan)[i])
    response <- as.numeric(as.matrix(log_metaphlan[,i])) 
    # response <- response * 1000000 #This portion is only needed if the poisson distribution is used.
    # response <- round(response, digits=0)
    if (sum(response) <= 30/1000000) next 
    return_list <- tryCatch({ 
      model_nofiber <- lmer(response ~ 1 + (1|Metagenomic_seq_plate) + Dose_model + eGFR_model + Hematocrit_model + High.Sensitivity.CRP_model + LDL..Calculated._model + RBC_model + Hemoglobin.A1c_model + Glucose..Ser.Plas_model + Triglyceride..Ser.Plas_model, data=log_metaphlan) #, family=poisson(link=identity) Poisson seems to work best at least for metagenome as the response variable
      model <- lmer(response ~ 1 + (1|Participant) + (1|Metagenomic_seq_plate) + eGFR_model + Dose_model + Hematocrit_model + High.Sensitivity.CRP_model + LDL..Calculated._model + RBC_model + Hemoglobin.A1c_model + Glucose..Ser.Plas_model + Triglyceride..Ser.Plas_model, data=log_metaphlan) #, family=poisson
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
      #if (is.na(summary(pbmmG2)$`P(>=obs)`)) {
      #  stats <- c(pval,sum(variance$vcov),colnames(log_metaphlan[i]),summary(pbmmG2)$`P(>=obs)`,summary(pbmmG2)$`#(>=obs)+noEst/runs`)
      #} else {stats <- c(pval,sum(variance$vcov),colnames(log_metaphlan[i]),summary(pbmmG2)$`P(>=obs)`,"NA")}
      
      stats <- c(pval,colnames(log_metaphlan[i]), coef(summary(model))[,1])
      #For matting thes "stats" object for the individual analyte
      stats <- data.frame(t(stats))
      colnames(stats) <- column_names
      stats <- data.frame(stats)
      
    },
    error = function(err){print(err)
      #Makins sure that if there's an error the program will keep running
      to_bind = c()
      stats = c()
    })
    
    #Appending the stats to the data frames,
    counter <- counter + 1
    print(counter)
    total_stats <- rbind(total_stats,stats)
  } #end data loop
  
  #Adding p-adjusted column, and polishing final data frame for output
  total_stats <- data.frame(total_stats)
  
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$eGFR)), method = "BH", n = nrow(total_stats))
  total_stats$eGFR_adj <- padjusted
  
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Dose)), method = "BH", n = nrow(total_stats))
  total_stats$Dose_adj <- padjusted
  
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Hematocrit)), method = "BH", n = nrow(total_stats))
  total_stats$Hematocrit_adj <- padjusted
  
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$High.Sensitivity.CRP)), method = "BH", n = nrow(total_stats))
  total_stats$HHigh.Sensitivity.CRP_adj <- padjusted
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$LDL.HDL.Ratio)), method = "BH", n = nrow(total_stats))
  total_stats$LDL.HDL.Ratio_adj <- padjusted
  
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Hemoglobin.A1c)), method = "BH", n = nrow(total_stats))
  total_stats$Hemoglobin.A1c_adj <- padjusted
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Glucose..Ser.Plas)), method = "BH", n = nrow(total_stats))
  total_stats$Glucose..Ser.Plas_adj <- padjusted
  
  padjusted <- p.adjust(as.numeric(as.matrix(total_stats$Triglyceride..Ser.Plas)), method = "BH", n = nrow(total_stats))
  total_stats$Triglyceride..Ser.Plas_adj <- padjusted
  
  rownames(total_stats) <- make.names(total_stats$analyte, unique=TRUE) #with the random effects interaction term
  total_stats <- subset(total_stats, select=-c(analyte))

    write.table(data.frame("run"=rownames(total_stats),total_stats), file = paste("~/metaphlan",fiber_subset,"-numericdose/numericStats",fiber_subset,"-numericdoseWeekPart.txt",sep=""), sep="\t",row.names=FALSE)
  #df <- read.table("/Users/SLancaster/Desktop/Ternary/Data/metaphlanVarianceMixWeekPartHematoHOMAIRTRI.txt", row.names=NULL)
  #df <- df[,-1] #Only to be used when reading a df from a table
  #df <- as.matrix(df) #Only to be used when reading a df from a table
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseDoseWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$Dose)), breaks=40)
  dev.off()
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseeGFRWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$eGFR)), breaks=40)
  dev.off()
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseHematocritWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$Hematocrit)), breaks=40)
  dev.off()
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseHigh.Sensitivity.CRP_modelWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$High.Sensitivity.CRP)), breaks=40)
  dev.off()
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseLDL..Calculated._model.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$LDL.HDL.Ratio)), breaks=40)
  dev.off()
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseHemoglobin.A1c_modelWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$Hemoglobin.A1c)), breaks=40)
  dev.off()
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseGlucose..Ser.Plas_modelWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$Glucose..Ser.Plas)), breaks=40)
  dev.off()   
  
  pdf(paste("~/metaphlan",fiber_subset,"-numericdose/numericHist",fiber_subset,"-numericdoseTriglyceride..Ser.Plas_modelWeekPart.pdf",sep=""))
  hist(as.numeric(as.matrix(total_stats$Triglyceride..Ser.Plas)), breaks=40)
  dev.off()
 
  alarm()
  
} #Ending the fiber_subset loop


#Additional Code
#Plot each individual analyte. Unnecessary for now.
plot <- ggplot(data=log_metaphlan, aes(x=log_metaphlan$Dose, y=response)) + #This seciton plots the data, so that we can better see the relationship.
  geom_point() +
  scale_x_discrete(limits=c("Baseline","10","20","30","WashoutD3","WashoutD10","WashoutFinal"))
ggsave(plot, file=paste("~/metaphlan",fiber_subset,"-numericdose/analyte-plots/numeric",fiber_subset,j,colnames(log_metaphlan)[i],".pdf",sep=""))


