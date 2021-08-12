
library(nlme)
# library(lme4)
library(lmerTest)
library(ggplot2)
library(scales)
library(reshape)

source("shared/helper.R")
source("shared/read_data_oct.R")

rna_new = read.csv("../data/2016_OF_RNA_NORM_FINAL.tsv", sep="\t")
# rna_new = read.csv("../data/2016_OVERFEED_RNA_VOOM.tsv", sep="\t")

colnames(rna_new) <- sub("^X$", "GeneName", colnames(rna_new) )
colnames(rna_new) <- sub("^X", "", colnames(rna_new) )
colnames(rna_new) <- sub("\\.", "-", colnames(rna_new) )
rna_data = rna_new[, -1 ]
rownames(rna_data) = rna_new$GeneName


pdf("../results/rna_density.pdf")
plot(density(rna_new[, 2] ), ylim=c(0, 0.2) )
for (idx in 2:ncol(rna_new) )
lines(density(rna_new[, idx] ) )
dev.off()

# Variance decomposition of RNA
# {{{

# Variance decomposition using only initial variables
# {{{
tt = sapply(1:nrow(rna_data), function(idx) {
  print(idx)
  tmp = rna_data[idx, ]
  tmp2 = data.frame(logRNA=as.numeric(tmp), sample=names(tmp))
  dfm = merge(tmp2, newmapping, by.x="sample", by.y="ID")
  # dfm[ (dfm$SSPG == 0.0 ), ]$SSPG <- mean(subset(dfm, SSPG != 0)$SSPG)
  lmodel = tryCatch({ lm (logRNA ~  1 + factor(patient_id) + sample_timepoint + deltaBMI + SSPG, dfm) }, error = function(e) {} )
  sq = summary(aov(lmodel))[[1]][, 'Sum Sq']
  rr = c(summary(lmodel)$r.squared, sq / sum(sq) * 100 ) 
  # if (length(rr) != 7 ) {return(na_return)}
  return(rr)
})
# 
# Create and store explanatory variables
varexpl = data.frame(t(tt))
colnames(varexpl) <- c("Rsq",  "Personal", "Timepoint", "deltaBMI", "SSPG", "Other")
varexpl$GeneName = rna_new$GeneName

for (clidx in 1:(ncol(varexpl)) )
{ print(sum(varexpl[, clidx]/ nrow(varexpl) ) ) }
# [1] 54.6829  # Personal    
# [1] 1.946821 # Timepoint   
# [1] 1.990333 # Weightdelta 
# [1] 1.748586 # SSPG        
# [1] 39.63136 # Unexplained 

write.csv(varexpl, "../results/rna_variable_explain_4var.csv", row.names=F)
varexpl = read.csv("../results/rna_variable_explain_4var.csv" )

for (clidx in 1:(ncol(varexpl)) )
{ print(median(varexpl[, clidx], na.rm=T ) ) }
# [1] 0.6118478
# [1] 55.26408
# [1] 1.370843
# [1] 1.188688
# [1] 0.9117657
# [1] 38.81522


# Only 5.9 % of variance is in any of the variables (on average)

## }}}

# varexpl = varexpl_6var [,-c(1, ncol(varexpl)-1,ncol(varexpl))]
varexpl_q = varexpl[,-c(1, ncol(varexpl))]
png("../results/variance_decomposition_3var_rna.png", width=4*8000, height=2000)
makeBarplotm(7.5, 3, varexpl_q[, c(2:4, 1, 5)]) 
dev.off()

# Ternary plot http://stackoverflow.com/questions/10879361/ternary-plot-and-filled-contour
varexpl_new = data.frame(Perturbation=varexpl$Timepoint + varexpl$SSPG + varexpl$deltaBMI, Other=varexpl$Other, Personal=varexpl$Personal)
## ternary plot
# varexpl_new = data.frame(Experiment=varexpl$Timepoint + varexpl$SSPG + varexpl$deltaBMI, 
#                          Other=varexpl$Other, 
#                          Personal=varexpl$Person)
varexpl_new$id <- 1:nrow(varexpl_new)
library(ggtern)
p = ggtern(data = varexpl_new,aes(x = Perturbation,y = Personal, z = Other)) + 
  # geom_density_tern(aes(fill=..level.., alpha=abs(..level..)),bins=500, binwidth=100) +
  geom_density_tern(aes(colour=..level..),bins=4500) +
  scale_fill_gradient(low="yellow",high="red") + 
  scale_color_gradient(low="yellow",high="red")   +
  # theme_bw()  +
  theme(legend.justification=c(0,1), legend.position=c(0,1)) + 
  guides(fill = guide_colorbar(order=1),
         alpha= guide_legend(order=2),
         color="none") + 
  labs(  title= "RNA",
         fill = "Value, V",alpha="|V - 0|")

pdf("../results/vd_tern_rna_3var.pdf", width=14)
print(p)
p = ggtern(data = varexpl_new,aes(x = Perturbation,y = Personal, z = Other)) + 
  geom_point(size=0.25)+ labs(title= "RNA") 
print(p)
dev.off()

# }}}


