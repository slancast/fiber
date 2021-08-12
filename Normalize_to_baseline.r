library(MASS)

pcl_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/relab-pathabundance.pcl", sep="\t", header=TRUE, row.names = 1)
bugs_list_df = read.csv("/Users/SLancaster/Desktop/Projects/Fiber/Microbiome/Data/Meta_genus_trimmed_NMDS_metadata.txt", sep="\t", header=TRUE, row.names = 1)

pcl_metadata_rows = 63
fiber_subset = "Arabinoxylan"
pcl_df <- data.frame(t(pcl_df))
pcl_df <- pcl_df[which(pcl_df$Fiber==fiber_subset),] #subsetting by fiber
pcl_df <- data.frame(t(pcl_df))
pcl_metadata = data.frame(head(pcl_df,pcl_metadata_rows))
bugs_list_df2 <- bugs_list_df[which(bugs_list_df$Fiber==fiber_subset),] #subsetting by fiber
bugs_list_df2 <- data.frame(t(bugs_list_df2))
bugs_list_df2 <- tail(bugs_list_df2, -pcl_metadata_rows) #Now getting rid of the metadata
bugs_list_df2 <- as.matrix(bugs_list_df2)
class(bugs_list_df2) <- "numeric"
bugs_list_df2 <- bugs_list_df2/100

pcl_df2 <- tail(pcl_df, -pcl_metadata_rows) #Now getting rid of the metadata
colnames(bugs_list_df2) <- colnames(pcl_df2)
pcl_df2 <- rbind(bugs_list_df2,pcl_df2)
pcl_df2 <- as.matrix(pcl_df2)
class(pcl_df2) <- "numeric"

pcl_df3 <- as.data.frame(pcl_metadata["Dose",])
pcl_df3 <- as.data.frame(t(pcl_df3))
partic <- as.data.frame(t(pcl_metadata["Participant",]))
pcl_df3 <- cbind(partic, pcl_df3)
pcl_df3 <- cbind(pcl_df3, data.frame(t(pcl_df2)))
pcl_df3 <- cbind(pcl_df3, data.frame(t(bugs_list_df2)))
pcl_df4 <- aggregate(pcl_df3,list(Participant=pcl_df3$Participant, Dose=pcl_df3$Dose),  mean)
pcl_df5 <- pcl_df4[,-3]
pcl_df5 <- pcl_df5[,-3]
baseline_means <- pcl_df5[ which(pcl_df5$Dose=='Baseline'),]
baseline_means <- data.frame(baseline_means)


pclfile = data.frame(t(pclfile))[,-2] #Creating a dataframe and making the files the same dimensions

match_df <- baseline_means[match(pclfile$Participant, baseline_means$Participant),na.omit(match(colnames(pclfile),colnames(baseline_means)))]
pclfile <- as.matrix(pclfile)
class(pclfile) <- "numeric"
pclfile <- pclfile[,-"Dose"]
pclfile <- pclfile[,-"Participant"]
match_df <- as.matrix(match_df)
class(match_df) <- "numeric"
match_df <- match_df[,-"Dose"]
match_df <- match_df[,-"Participant"]

normalized_df = pclfile - match_df

normalzied_pclfile <- rbind(metadata, t(normalized_df))