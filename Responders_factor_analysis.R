#/user/bin/R
#Here I will use the factor analysis on the responders vs non-responders to see how the factors are different between the two
#It looks like I can do a ttest between the factor scores.

source("/Users/SLancaster/Desktop/Projects/Fiber/Multiomics/Software/responders_vs_nonresponders.R")

ttest_topvsbottom <- function(responders,baseline_means,metadata){
  top_third <- responders$Participant[1:6]
  bottom_third <- responders$Participant[13:18]
  baseline_means[,1] <- as.character(baseline_means[,1] )
  class(baseline_means[,1]) <- "numeric"
  top_third_data <- baseline_means[which(baseline_means[,1] %in% top_third),]
  bottom_third_data <- baseline_means[which(baseline_means[,1] %in% bottom_third),]
  top_third_data <- as.matrix(top_third_data[,-2])
  bottom_third_data <- as.matrix(bottom_third_data[,-2])
  class(top_third_data) <- "numeric"
  class(bottom_third_data) <- "numeric"

  require('psych')
  
  combined <- rbind(top_third_data, bottom_third_data)
  combined_nopart <- combined[,-1]
  combined_nopart <- data.frame(t(combined_nopart ))
  
  #Ting the most abundant ones in case the lower ones are stochastically weird
  combined_nopart2 <- combined_nopart[order(rowSums(combined_nopart), decreasing = TRUE),]
  #combined_nopart2 <- combined_nopart2[1:100,]
  
  combined_parallel <- fa.parallel(combined_nopart2, fm = 'minres', fa = 'fa')
  
  combined_factor <- fa(combined_nopart2,nfactors = 5,rotate = "oblimin",fm="minres")
  
  groups <- data.frame(c(rep(1,6),rep(2,6)))
  
  plotting <- cbind(combined_factor$weights,Reponder=groups)
  # 
  # ggplot() +
  # geom_point(plotting, mapping = aes(x=MR1, y= MR2, colour=as.factor(plotting[,4]) ))
  # 
  fa_group1 <- plotting[1:6,]
  fa_group2 <- plotting[7:12,]
  
  #Why are the genes grouped together in latent variables?
  pvalues <- data.frame(analyte=character(),pvalue=double(),top_average=double(),bottom_average=double())
  pvalues <- as.matrix(pvalues)
  for (i in 1:(ncol(fa_group2)-1)){
    print(colnames(fa_group2)[i])
    tested <- t.test(fa_group1[,i],fa_group2[,i])
    print(mean(fa_group1[,i]))
    print(tested$p.value)
    pvalues <- rbind(pvalues, c(colnames(bottom_third_data)[i], tested$p.value, mean(fa_group1[,i]), mean(fa_group2[,i])))
  }
  
}


# Following will is some work I've already done for the lda 

groups <- data.frame(t(c(rep(1,6),rep(2,6))))
colnames(groups) <- colnames(combined_nopart)

combined_groupings <- combined_nopart[order(rowSums(combined_nopart),decreasing = TRUE),]
combined_groupings <- combined_groupings[1:100,]

combined_groupings <- rbind(responder =groups,combined_groupings)
combined_groupings <- data.frame(t(combined_groupings ))
combined_groupings <- combined_groupings*100000000 #this LDA protocol doesn't like the small values

lda_results <- lda(formula=responder ~ ., data=combined_groupings)

train <- sample(1:100, 50)

plda = predict(object = lda_results, # predictions
               newdata = combined_groupings[,-train ])


