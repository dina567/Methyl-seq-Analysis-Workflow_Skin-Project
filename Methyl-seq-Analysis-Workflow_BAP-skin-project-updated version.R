# ------ Load packages ------- 
library(dplyr)
library(DSS)
library(bsseq)
# install.packages("gtools") ## Uncomment if not already installed
library(gtools)

packageVersion("ChIPseeker")
packageVersion("DSS")
citation("DSS")
# ------ Read in Data ------ 
getwd()
setwd("C:/Rutgers/Konglab/Projects/HCT116-Delphinidin/HCT_methyl-seq_oct-2019/output")

data <- read.table(file = "C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/raw data/combined_r5s5d50x300_anno-2022-04-04.tsv", sep = "\t", header = TRUE)
data <- read.table(file = "HCT_combined_def_anno.tsv", sep = "\t", header = TRUE)

data[is.na(data)] <- 0
colnames(data)

data.copy <- data

## becuase start position contains duplicated value !!!
any(duplicated(data$start)) #TRUE

# solution 1: assigned start position from 1 as a gene ID (DMLtest Does not work, even if no duplicated)
# 
# data[, 2] <- seq(1:length(data$start))
# any(duplicated(data$start)) #FALSE

# solution 2: assigned mean of the start and end position from 1 as a gene ID (Does not work)

# data <- mutate(data, pos = (data$start+data$end)/2)
# any(duplicated(data$pos)) #TRUE
# any(duplicated(data))
# colnames(data)
# data <- data[,c("seqnames", "start", "end", "pos")]
# colnames(data)
# duplicate <- data[which(duplicated(data$pos)),]

# solution 3: plus 1 for all the duplicated start position 

any(duplicated(data))
colnames(data)
duplicate <- data[which(duplicated(data$start)),]
non_duplicate <- data[which(!duplicated(data$start)),]
duplicate$start <- duplicate$start +1
data <- rbind (non_duplicate, duplicate)
any(duplicated(data$start)) #FALSE

# (optional) rearrange the order of start and chr
data <- arrange(data, start)
data <- arrange(data, seqnames)

data.copy <- data

#(optional) Download the new table with newly assigned start position (if needed)
write.table(data, "C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/combined_r5s5d50x300_anno-2022-04-04_ordered_2022-08-05-no duplicated start pos.tsv",
           sep = "\t")

## (optional) sort chr column by order
# data$chr <- as.vector(sub("chr", "", data$chr))
# data$chr <- as.vector(data$chr)
# mixedsort(data$chr)
# 
# data <- data[order(data$chr), ]
# 
# data[mixedsort(data$chr),]
# 
# data <- data[order(data$chr), ]
# data$chr <- paste("chr", data$chr, sep = "")

# ------ Read in Design Table ------
design <- read.csv("C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/design table_Methyl-seq-07-16-22.csv", sep = ",", header = TRUE)
colnames(design)
design <- data.frame(design[,3:4])
design

# ------ Prepare Data Subsets with Different Conditions ------
# ------Subset the data into separated lists--------
# the data contain 4 variables (chr, pos, N, X)
# sample ID: RW01, RW02, ...., to RW22 

# subset using select() from dplyr package or grepl()
# (1) method 1: select(data,contains("RW01"))
# (2) method 2: grepl()

# subset data for delphinidin
C_del <- cbind(data[,1:2], data[,grepl("HCT.D", colnames(data))])

# subset data for BAP UA in vivo study
C_5w_1 <- cbind(data[,1:2], data[,grepl("RW01", colnames(data))])
C_5w_2 <- cbind(data[,1:2],data[,grepl("RW02", colnames(data))])
C_5w_3 <- cbind(data[,1:2],data[,grepl("RW03", colnames(data))])
C_20w_1 <- cbind(data[,1:2],data[,grepl("RW04", colnames(data))])
C_20w_2 <- cbind(data[,1:2],data[,grepl("RW05", colnames(data))])
C_20w_3 <- cbind(data[,1:2],data[,grepl("RW06", colnames(data))])
BAP_5w_1 <- cbind(data[,1:2],data[,grepl("RW07", colnames(data))])
BAP_5w_2 <- cbind(data[,1:2],data[,grepl("RW08", colnames(data))])
BAP_5w_3 <- cbind(data[,1:2],data[,grepl("RW09", colnames(data))])
BAP_20w_1 <- cbind(data[,1:2],data[,grepl("RW10", colnames(data))])
BAP_20w_2 <- cbind(data[,1:2],data[,grepl("RW11", colnames(data))])
BAP_20w_3 <- cbind(data[,1:2],data[,grepl("RW12", colnames(data))])
BAP_26w_1 <- cbind(data[,1:2],data[,grepl("RW13", colnames(data))])
BAP_26w_2 <- cbind(data[,1:2],data[,grepl("RW14", colnames(data))])
UA_5w_1 <- cbind(data[,1:2],data[,grepl("RW15", colnames(data))])
UA_5w_2 <- cbind(data[,1:2],data[,grepl("RW16", colnames(data))])
UA_5w_3 <- cbind(data[,1:2],data[,grepl("RW17", colnames(data))])
UA_20w_1 <- cbind(data[,1:2],data[,grepl("RW18", colnames(data))])
UA_20w_2 <- cbind(data[,1:2],data[,grepl("RW19", colnames(data))])
UA_20w_3 <- cbind(data[,1:2],data[,grepl("RW20", colnames(data))])
UA_26w_1 <- cbind(data[,1:2],data[,grepl("RW21", colnames(data))])
UA_26w_2 <- cbind(data[,1:2],data[,grepl("RW22", colnames(data))])
# C_5w_1_allcol <- mutate(data[,!grepl("RW", colnames(data))],data[,grepl("RW01", colnames(data))])

# Combined all the separated data into a list 
Control_list <- list(C_5w_1, C_5w_2, C_5w_3, C_20w_1, C_20w_2, C_20w_3)
BAP_list <- list(BAP_5w_1, BAP_5w_2, BAP_5w_3, BAP_20w_1, BAP_20w_2, BAP_20w_3)
BAP_list_26wk <- list(BAP_5w_1, BAP_5w_2, BAP_5w_3, BAP_20w_1, BAP_20w_2, BAP_20w_3, BAP_26w_1, BAP_26w_2)
UA_list <- list(UA_5w_1, UA_5w_2, UA_5w_3, UA_20w_1, UA_20w_2, UA_20w_3)
UA_list_26wk <- list(UA_5w_1, UA_5w_2, UA_5w_3, UA_20w_1, UA_20w_2, UA_20w_3, UA_26w_1, UA_26w_2)
data_list <- list(C_5w_1, C_5w_2, C_5w_3, C_20w_1, C_20w_2, C_20w_3, BAP_5w_1, BAP_5w_2, BAP_5w_3, BAP_20w_1, BAP_20w_2, BAP_20w_3, UA_5w_1, UA_5w_2, UA_5w_3, UA_20w_1, UA_20w_2, UA_20w_3)
data_list_26wk <- list(C_5w_1, C_5w_2, C_5w_3, C_20w_1, C_20w_2, C_20w_3, BAP_5w_1, BAP_5w_2, BAP_5w_3, BAP_20w_1, BAP_20w_2, BAP_20w_3, BAP_26w_1, BAP_26w_2, UA_5w_1, UA_5w_2, UA_5w_3, UA_20w_1, UA_20w_2, UA_20w_3, UA_26w_1, UA_26w_2)
length(data_list)
length(data_list_26wk)

# Assign the same column names to all the separated list
for (i in 1:length(data_list)) {
  colnames(data_list[[i]]) <- c("chr", "pos", "N", "X") 
  i = i + 1
}

colnames(data_list[[18]])

for (i in 1:length(data_list_26wk)) {
  colnames(data_list_26wk[[i]]) <- c("chr", "pos", "N", "X") 
  i = i + 1
}
colnames(data_list_26wk[[22]])

# ------ Run DSS ------
Exp = makeBSseqData(data_list_26wk,
                    c("C_5w_1", "C_5w_2", "C_5w_3", 
                      "C_20w_1", "C_20w_2", "C_20w_3", 
                      "BAP_5w_1", "BAP_5w_2", "BAP_5w_3", 
                      "BAP_20w_1", "BAP_20w_2", "BAP_20w_3", 
                      "BAP_26w_1", "BAP_26w_2", 
                      "UA_5w_1", "UA_5w_2", "UA_5w_3", 
                      "UA_20w_1", "UA_20w_2", "UA_20w_3", 
                      "UA_26w_1", "UA_26w_2"))

head(data_list_26wk)
warnings() # even if the pos is ordered, still get CG positions in chromosome chr1...chrY is not ordered. Reorder CG sites.

# Exp_5w = makeBSseqData(data_list_26wk[[c(1:3,7:9, 15:17)]],
#                        c("C_5w_1", "C_5w_2", "C_5w_3", 
#                          "BAP_5w_1", "BAP_5w_2", "BAP_5w_3",
#                          "UA_5w_1", "UA_5w_2", "UA_5w_3"))
# 
# Exp_20w = makeBSseqData(data_list_26wk[[c(4:6,10:12, 18:20)]],
#                        c("C_20w_1", "C_20w_2", "C_20w_3", 
#                          "BAP_20w_1", "BAP_20w_2", "BAP_20w_3", 
#                          "UA_20w_1", "UA_20w_2", "UA_20w_3"))
# Exp_26w = makeBSseqData(data_list_26wk[[c(4:6, 13:14, 21:22)]],
#                         c("C_20w_1", "C_20w_2", "C_20w_3",
#                           "BAP_26w_1", "BAP_26w_2", 
#                           "UA_26w_1", "UA_26w_2"))

# 1. BAP+TPA vs C

dmlTest_BAPvsC_all = DMLtest(Exp, 
                             group1=c("C_5w_1", "C_5w_2", "C_5w_3", "C_20w_1", "C_20w_2", "C_20w_3"), 
                             group2=c("BAP_5w_1", "BAP_5w_2", "BAP_5w_3", "BAP_20w_1", "BAP_20w_2", "BAP_20w_3","BAP_26w_1", "BAP_26w_2"),
                             smoothing = TRUE, 
                             smoothing.span = 500)
head(dmlTest_BAPvsC_all)

dmlTest_BAPvsC_5w = DMLtest(Exp, 
                            group1=c("C_5w_1", "C_5w_2", "C_5w_3"), 
                            group2=c("BAP_5w_1", "BAP_5w_2", "BAP_5w_3"),
                            smoothing = TRUE, 
                            smoothing.span = 500)
head(dmlTest_BAPvsC_5w)
write.csv(dmlTest_BAPvsC_5w, "DSS_result_BAPvsC_5w.csv")

dmlTest_BAPvsC_20w = DMLtest(Exp, 
                            group1=c("C_20w_1", "C_20w_2", "C_20w_3"), 
                            group2=c("BAP_20w_1", "BAP_20w_2", "BAP_20w_3"),
                            smoothing = TRUE, 
                            smoothing.span = 500)
head(dmlTest_BAPvsC_20w)
write.csv(dmlTest_BAPvsC_20w, "DSS_result_BAPvsC_20w.csv")

dmlTest_BAPvsC_26w = DMLtest(Exp, 
                             group1=c("C_20w_1", "C_20w_2", "C_20w_3"), 
                             group2=c("BAP_26w_1", "BAP_26w_2"),
                             smoothing = TRUE, 
                             smoothing.span = 500)
head(dmlTest_BAPvsC_26w)
write.csv(dmlTest_BAPvsC_26w, "DSS_result_BAPvsC_26w.csv")

# 2. UA vs BAP+TPA

dmlTest_UAvsBAP_all = DMLtest(Exp, 
                             group1=c("BAP_5w_1", "BAP_5w_2", "BAP_5w_3", "BAP_20w_1", "BAP_20w_2", "BAP_20w_3","BAP_26w_1", "BAP_26w_2"), 
                             group2=c("UA_5w_1", "UA_5w_2", "UA_5w_3", "UA_20w_1", "UA_20w_2", "UA_20w_3", "UA_26w_1", "UA_26w_2"),
                             smoothing = TRUE, 
                             smoothing.span = 500)
head(dmlTest_UAvsBAP_all)
write.csv(dmlTest_UAvsBAP_all, "DSS_result_UAvsBAP_all.csv")

dmlTest_UAvsBAP_5w = DMLtest(Exp, 
                             group1=c("BAP_5w_1", "BAP_5w_2", "BAP_5w_3"), 
                             group2=c("UA_5w_1", "UA_5w_2", "UA_5w_3"),
                             smoothing = TRUE, 
                             smoothing.span = 500)
head(dmlTest_UAvsBAP_5w)
write.csv(dmlTest_UAvsBAP_5w, "DSS_result_UAvsBAP_5w.csv")

dmlTest_UAvsBAP_20w = DMLtest(Exp, 
                              group1=c("BAP_20w_1", "BAP_20w_2", "BAP_20w_3"), 
                              group2=c("UA_20w_1", "UA_20w_2", "UA_20w_3"),
                              smoothing = TRUE, 
                              smoothing.span = 500)
head(dmlTest_UAvsBAP_20w)
write.csv(dmlTest_UAvsBAP_20w, "DSS_result_UAvsBAP_20w.csv")

dmlTest_UAvsBAP_26w = DMLtest(Exp, 
                              group1=c("BAP_26w_1", "BAP_26w_2"), 
                              group2=c("UA_26w_1", "UA_26w_2"),
                              smoothing = TRUE, 
                              smoothing.span = 500)
head(dmlTest_UAvsBAP_26w)
write.csv(dmlTest_UAvsBAP_26w, "DSS_result_UAvsBAP_26w.csv")

## Note for equal.disp:
# if equal dispersion is assumed, put equal.disp = TRUE in DMLtest 
# A flag to indicate whether the dispersion in two groups are deemed equal. 
# Default is FALSE, and the dispersion shrinkages are performed on two conditions independently
# When there is no biological replicate in one or both treatment groups, users can either 
#(1) specify equal.disp=TRUE, 
# which assumes both groups have the same dispersion, then the data from two groups are combined and used as replicates to estimate dispersion; or 
#(2) specify smoothing=TRUE,
# which uses the smoothed means (methylation levels) to estimate dispersions via a shrinkage estimator. 
# This smoothing procedure uses data from neighboring CpG sites as "pseudo-replicate" for estimating biological variance


## ------Explanation of parameters/values in DML result-------
# chr Chromosome number.
# pos Genomic coordinates.
# mu1, mu2 Mean methylations of two groups.
# diff Difference of mean methylations of two groups. "diff=mu1-mu2". !!! (be careful!!!)
# diff.se Standard error of the methylation difference.
# stat Wald statistics.
# pval P-values. This is obtained from normal distribution.
# phi1, phi2 Estimated dispersions in two groups.
# pval P-values. This is obtained from normal distribution.
# fdr False discovery rate.
# postprob.overThreshold: The posterior probability of the difference in methylation greater than delta. This columns is only available when delta>0.

# ------ Fliter and Visualize results of DMLs or DMRs -------
# callDML or callDMR 
# By default, the test is based on the null hypothesis that the difference in methylation levels is 0. 
# Users can specify a threshold for difference by delta = __
# When delta is specified, the function will compute the posterior probability that the difference of the means is greater than delta. 
# the threshold for p-value here actually refers to the threshold for 1-posterior probability, or the local FDR. Here we use the same parameter name for the sake of the consistence of function syntax.

## method 1 -------
# combined the gene names and annotation from original data with the DML analysis results--------

# 1. BAP+TPA vs C
# (1) combined data
dmlTest_BAPvsC_all_copy <- dmlTest_BAPvsC_all

any(duplicated(dmlTest_BAPvsC_all))
any(duplicated(dmlTest_BAPvsC_all$pos))
data <- arrange(data, start)
dmlTest_BAPvsC_all <- dmlTest_BAPvsC_all[order(dmlTest_BAPvsC_all$pos),]

geneID_BAPvsC_all <- data[which(data$start %in% dmlTest_BAPvsC_all$pos),]
colnames(geneID_BAPvsC_all)

DML_geneID_BAPvsC_all <- cbind(geneID_BAPvsC_all[,c("geneId", "annotation", "seqnames","start","end")], 
                               dmlTest_BAPvsC_all)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_BAPvsC_all$identical = ifelse(DML_geneID_BAPvsC_all$start == DML_geneID_BAPvsC_all$pos, 'T','F')
any(DML_geneID_BAPvsC_all$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_BAPvsC_all)
DML_geneID_BAPvsC_all <- DML_geneID_BAPvsC_all[,-c(6, 7,17)]
DML_geneID_BAPvsC_all_0.1 <- filter(DML_geneID_BAPvsC_all, fdr <= 0.1)
DML_geneID_BAPvsC_all_0.05 <- filter(DML_geneID_BAPvsC_all, fdr <= 0.05)
write.csv(DML_geneID_BAPvsC_all, "DSS_result_BAPvsC_all_anno.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_all_0.1, "DSS_result_BAPvsC_all_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_all_0.05, "DSS_result_BAPvsC_all_anno_fdr0.05.csv", row.names = T)

## 5w
# (1) combined data
dmlTest_BAPvsC_5w_copy <- dmlTest_BAPvsC_5w

any(duplicated(dmlTest_BAPvsC_5w))
any(duplicated(dmlTest_BAPvsC_5w$pos))
data <- arrange(data, start)
dmlTest_BAPvsC_5w <- dmlTest_BAPvsC_5w[order(dmlTest_BAPvsC_5w$pos),]

geneID_BAPvsC_5w <- data[which(data$start %in% dmlTest_BAPvsC_5w$pos),]
colnames(geneID_BAPvsC_5w)

DML_geneID_BAPvsC_5w <- cbind(geneID_BAPvsC_5w[,c("geneId", "annotation", "seqnames","start","end")], 
                              dmlTest_BAPvsC_5w)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_BAPvsC_5w$identical = ifelse(DML_geneID_BAPvsC_5w$start == DML_geneID_BAPvsC_5w$pos, 'T','F')
any(DML_geneID_BAPvsC_5w$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_BAPvsC_5w)

DML_geneID_BAPvsC_5w <- DML_geneID_BAPvsC_5w[,-c(6, 7,17)]
DML_geneID_BAPvsC_5w_0.1 <- filter(DML_geneID_BAPvsC_5w, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_BAPvsC_5w_0.05 <- filter(DML_geneID_BAPvsC_5w, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_BAPvsC_5w, "DSS_result_BAPvsC_5w_anno.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_5w_0.1, "DSS_result_BAPvsC_5w_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_5w_0.05, "DSS_result_BAPvsC_5w_anno_fdr0.05.csv", row.names = T)

## 20w
# (1) combined data
dmlTest_BAPvsC_20w_copy <- dmlTest_BAPvsC_20w

any(duplicated(dmlTest_BAPvsC_20w))
any(duplicated(dmlTest_BAPvsC_20w$pos))
data <- arrange(data, start)
dmlTest_BAPvsC_20w <- dmlTest_BAPvsC_20w[order(dmlTest_BAPvsC_20w$pos),]

geneID_BAPvsC_20w <- data[which(data$start %in% dmlTest_BAPvsC_20w$pos),]
colnames(geneID_BAPvsC_20w)

DML_geneID_BAPvsC_20w <- cbind(geneID_BAPvsC_20w[,c("geneId", "annotation", "seqnames","start","end")], 
                               dmlTest_BAPvsC_20w)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_BAPvsC_20w$identical = ifelse(DML_geneID_BAPvsC_20w$start == DML_geneID_BAPvsC_20w$pos, 'T','F')
any(DML_geneID_BAPvsC_20w$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_BAPvsC_20w)

DML_geneID_BAPvsC_20w <- DML_geneID_BAPvsC_20w[,-c(6, 7,17)]
DML_geneID_BAPvsC_20w_0.1 <- filter(DML_geneID_BAPvsC_20w, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_BAPvsC_20w_0.05 <- filter(DML_geneID_BAPvsC_20w, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_BAPvsC_20w, "DSS_result_BAPvsC_20w_anno.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_20w_0.1, "DSS_result_BAPvsC_20w_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_20w_0.05, "DSS_result_BAPvsC_20w_anno_fdr0.05.csv", row.names = T)

## 26w
# (1) combined data
dmlTest_BAPvsC_26w_copy <- dmlTest_BAPvsC_26w

any(duplicated(dmlTest_BAPvsC_26w))
any(duplicated(dmlTest_BAPvsC_26w$pos))
data <- arrange(data, start)
dmlTest_BAPvsC_26w <- dmlTest_BAPvsC_26w[order(dmlTest_BAPvsC_26w$pos),]

geneID_BAPvsC_26w <- data[which(data$start %in% dmlTest_BAPvsC_26w$pos),]
colnames(geneID_BAPvsC_26w)

DML_geneID_BAPvsC_26w <- cbind(geneID_BAPvsC_26w[,c("geneId", "annotation", "seqnames","start","end")], 
                               dmlTest_BAPvsC_26w)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_BAPvsC_26w$identical = ifelse(DML_geneID_BAPvsC_26w$start == DML_geneID_BAPvsC_26w$pos, 'T','F')
any(DML_geneID_BAPvsC_26w$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_BAPvsC_26w)

DML_geneID_BAPvsC_26w <- DML_geneID_BAPvsC_26w[,-c(6, 7,17)]
DML_geneID_BAPvsC_26w_0.1 <- filter(DML_geneID_BAPvsC_26w, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_BAPvsC_26w_0.05 <- filter(DML_geneID_BAPvsC_26w, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_BAPvsC_26w, "DSS_result_BAPvsC_26w_anno.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_26w_0.1, "DSS_result_BAPvsC_26w_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_BAPvsC_26w_0.05, "DSS_result_BAPvsC_26w_anno_fdr0.05.csv", row.names = T)

####

#2. UA vs BAP+TPA
# (1) combined data
dmlTest_UAvsBAP_all_copy <- dmlTest_UAvsBAP_all

any(duplicated(dmlTest_UAvsBAP_all))
any(duplicated(dmlTest_UAvsBAP_all$pos))
data <- arrange(data, start)
dmlTest_UAvsBAP_all <- dmlTest_UAvsBAP_all[order(dmlTest_UAvsBAP_all$pos),]

geneID_UAvsBAP_all <- data[which(data$start %in% dmlTest_UAvsBAP_all$pos),]
colnames(geneID_UAvsBAP_all)

DML_geneID_UAvsBAP_all <- cbind(geneID_UAvsBAP_all[,c("geneId", "annotation", "seqnames","start","end")], 
                               dmlTest_UAvsBAP_all)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_UAvsBAP_all$identical = ifelse(DML_geneID_UAvsBAP_all$start == DML_geneID_UAvsBAP_all$pos, 'T','F')
any(DML_geneID_UAvsBAP_all$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_UAvsBAP_all)

DML_geneID_UAvsBAP_all <- DML_geneID_UAvsBAP_all[,-c(6, 7,17)]
DML_geneID_UAvsBAP_all_0.1 <- filter(DML_geneID_UAvsBAP_all, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_UAvsBAP_all_0.05 <- filter(DML_geneID_UAvsBAP_all, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_UAvsBAP_all, "DSS_result_UAvsBAP_all_anno.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_all_0.1, "DSS_result_UAvsBAP_all_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_all_0.05, "DSS_result_UAvsBAP_all_anno_fdr0.05.csv", row.names = T)
## 5w
# (1) combined data
dmlTest_UAvsBA_5w_copy <- dmlTest_UAvsBAP_5w

any(duplicated(dmlTest_UAvsBAP_5w))
any(duplicated(dmlTest_UAvsBAP_5w$pos))
data <- arrange(data, start)
dmlTest_UAvsBAP_5w <- dmlTest_UAvsBAP_5w[order(dmlTest_UAvsBAP_5w$pos),]

geneID_UAvsBAP_5w <- data[which(data$start %in% dmlTest_UAvsBAP_5w$pos),]
colnames(geneID_UAvsBAP_5w)

DML_geneID_UAvsBAP_5w <- cbind(geneID_UAvsBAP_5w[,c("geneId", "annotation", "seqnames","start","end")], 
                              dmlTest_UAvsBAP_5w)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_UAvsBAP_5w$identical = ifelse(DML_geneID_UAvsBAP_5w$start == DML_geneID_UAvsBAP_5w$pos, 'T','F')
any(DML_geneID_UAvsBAP_5w$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_UAvsBAP_5w)

DML_geneID_UAvsBAP_5w <- DML_geneID_UAvsBAP_5w[,-c(6, 7,17)]
DML_geneID_UAvsBAP_5w_0.1 <- filter(DML_geneID_UAvsBAP_5w, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_UAvsBAP_5w_0.05 <- filter(DML_geneID_UAvsBAP_5w, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_UAvsBAP_5w, "DSS_result_UAvsBAP_5w_anno.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_5w_0.1, "DSS_result_UAvsBAP_5w_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_5w_0.05, "DSS_result_UAvsBAP_5w_anno_fdr0.05.csv", row.names = T)

## 20w
# (1) combined data
dmlTest_UAvsBAP_20w_copy <- dmlTest_UAvsBAP_20w

any(duplicated(dmlTest_UAvsBAP_20w))
any(duplicated(dmlTest_UAvsBAP_20w$pos))
data <- arrange(data, start)
dmlTest_UAvsBAP_20w <- dmlTest_UAvsBAP_20w[order(dmlTest_UAvsBAP_20w$pos),]

geneID_UAvsBAP_20w <- data[which(data$start %in% dmlTest_UAvsBAP_20w$pos),]
colnames(geneID_UAvsBAP_20w)

DML_geneID_UAvsBAP_20w <- cbind(geneID_UAvsBAP_20w[,c("geneId", "annotation", "seqnames","start","end")], 
                               dmlTest_UAvsBAP_20w)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_UAvsBAP_20w$identical = ifelse(DML_geneID_UAvsBAP_20w$start == DML_geneID_UAvsBAP_20w$pos, 'T','F')
any(DML_geneID_UAvsBAP_20w$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_UAvsBAP_20w)
DML_geneID_UAvsBAP_20w <- DML_geneID_UAvsBAP_20w[,-c(6, 7,17)]
DML_geneID_UAvsBAP_20w_0.1 <- filter(DML_geneID_UAvsBAP_20w, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_UAvsBAP_20w_0.05 <- filter(DML_geneID_UAvsBAP_20w, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_UAvsBAP_20w, "DSS_result_UAvsBAP_20w_anno.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_20w_0.1, "DSS_result_UAvsBAP_20w_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_20w_0.05, "DSS_result_UAvsBAP_20w_anno_fdr0.05.csv", row.names = T)

## 26w
# (1) combined data
dmlTest_UAvsBAP_26w_copy <- dmlTest_UAvsBAP_26w

any(duplicated(dmlTest_UAvsBAP_26w))
any(duplicated(dmlTest_UAvsBAP_26w$pos))
data <- arrange(data, start)
dmlTest_UAvsBAP_26w <- dmlTest_UAvsBAP_26w[order(dmlTest_UAvsBAP_26w$pos),]

geneID_UAvsBAP_26w <- data[which(data$start %in% dmlTest_UAvsBAP_26w$pos),]
colnames(geneID_UAvsBAP_26w)

DML_geneID_UAvsBAP_26w <- cbind(geneID_UAvsBAP_26w[,c("geneId", "annotation", "seqnames","start","end")], 
                               dmlTest_UAvsBAP_26w)

# (2) check if all the "start" values match the "pos" values from DML analysis
DML_geneID_UAvsBAP_26w$identical = ifelse(DML_geneID_UAvsBAP_26w$start == DML_geneID_UAvsBAP_26w$pos, 'T','F')
any(DML_geneID_UAvsBAP_26w$identical == FALSE) 

# (3) save the DML result with gene/annotation
colnames(DML_geneID_UAvsBAP_26w)

DML_geneID_UAvsBAP_26w <- DML_geneID_UAvsBAP_26w[,-c(6, 7,17)]
DML_geneID_UAvsBAP_26w_0.1 <- filter(DML_geneID_UAvsBAP_26w, fdr <= 0.1) # filter data by fdr <=0.1
DML_geneID_UAvsBAP_26w_0.05 <- filter(DML_geneID_UAvsBAP_26w, fdr <= 0.05) # filter data by fdr <=0.05
write.csv(DML_geneID_UAvsBAP_26w, "DSS_result_UAvsBAP_26w_anno.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_26w_0.1, "DSS_result_UAvsBAP_26w_anno_fdr0.1.csv", row.names = T)
write.csv(DML_geneID_UAvsBAP_26w_0.05, "DSS_result_UAvsBAP_26w_anno_fdr0.05.csv", row.names = T)

# methylation difference heatmap-------

library(ggfortify)

setwd("C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/DSS results/non annotated")

w5_TPA <- read.csv("DSS_result_BAPvsC_5w_anno.csv")
w20_TPA <- read.csv("DSS_result_BAPvsC_20w_anno.csv")
w26_TPA <- read.csv("DSS_result_BAPvsC_26w_anno.csv")
w5_UA <- read.csv("DSS_result_UAvsBAP_5w_anno.csv")
w20_UA <- read.csv("DSS_result_UAvsBAP_20w_anno.csv")
w26_UA <-read.csv("DSS_result_UAvsBAP_26w_anno.csv")

w5_TPA <- read.csv("DSS_result_BAPvsC_5w.csv")
w20_TPA <- read.csv("DSS_result_BAPvsC_20w.csv")
w26_TPA <- read.csv("DSS_result_BAPvsC_26w.csv")
w5_UA <- read.csv("DSS_result_UAvsBAP_5w.csv")
w20_UA <- read.csv("DSS_result_UAvsBAP_20w.csv")
w26_UA <-read.csv("DSS_result_UAvsBAP_26w.csv")

library(dplyr)
colnames(w5_TPA)
w5_TPA <- w5_TPA[,c(3,6,12)] %>% filter(fdr <= 0.1 & diff <= -0.1|diff >= 0.1)
w20_TPA <- w20_TPA[,c(3,6,12)] %>% filter(fdr <= 0.1& diff <= -0.1|diff >= 0.1)
w26_TPA <- w26_TPA[,c(3,6,12)] %>% filter(fdr <= 0.1& diff <= -0.1|diff >= 0.1)
w5_UA <- w5_UA[,c(3,6,12)] %>% filter(fdr <= 0.1& diff <= -0.1|diff >= 0.1)
w20_UA <- w20_UA[,c(3,6,12)] %>% filter(fdr <= 0.1& diff <= -0.1|diff >= 0.1)
w26_UA <- w26_UA[,c(3,6,12)] %>% filter(fdr <= 0.1& diff <= -0.1|diff >= 0.1)

#
w5_TPA <- w5_TPA[,c(3,6,12)] %>% filter(fdr <= 0.1)
w20_TPA <- w20_TPA[,c(3,6,12)] %>% filter(fdr <= 0.1)
w26_TPA <- w26_TPA[,c(3,6,12)] %>% filter(fdr <= 0.1)
w5_UA <- w5_UA[,c(3,6,12)] %>% filter(fdr <= 0.1)
w20_UA <- w20_UA[,c(3,6,12)] %>% filter(fdr <= 0.1)
w26_UA <- w26_UA[,c(3,6,12)] %>% filter(fdr <= 0.1)

#
w5_TPA <- w5_TPA[,c(3,6,12)]
w20_TPA <- w20_TPA[,c(3,6,12)]
w26_TPA <- w26_TPA[,c(3,6,12)]
w5_UA <- w5_UA[,c(3,6,12)]
w20_UA <- w20_UA[,c(3,6,12)]
w26_UA <- w26_UA[,c(3,6,12)]
colnames(w5_TPA) <- c("pos", "diff_w5_TPA", "fdr_w5_TPA")
colnames(w20_TPA) <- c("pos", "diff_w20_TPA", "fdr_w20_TPA")
colnames(w26_TPA) <- c("pos", "diff_w26_TPA", "fdr_w26_TPA")
colnames(w5_UA) <- c("pos", "diff_w5_UA", "fdr_w5_UA")
colnames(w20_UA) <- c("pos", "diff_w20_UA", "fdr_w20_UA")
colnames(w26_UA) <- c("pos", "diff_w26_UA", "fdr_w26_UA")

df <- merge(w5_TPA, w20_TPA, by = "pos")
df <- merge(df, w26_TPA, by = "pos")
df <- merge(df, w5_UA, by = "pos")
df <- merge(df, w20_UA, by = "pos")
df <- merge(df, w26_UA, by = "pos")

head(df)
df <- df[,c(2, 4, 6, 8, 10, 12)]
colnames(df)
anno <- c("BAP+TPA vs. C","BAP+TPA vs. C", "BAP+TPA vs. C", "UA vs. BAP+TPA", "UA vs. BAP+TPA",  "UA vs. BAP+TPA")
anno <- as.data.frame(anno) 

rownames(anno) <- colnames(df)
anno$wk <- c("5", "20", "26","5", "20", "26")
colnames(anno) <- c("trt", "wk")
anno$wk <- factor(anno$wk, levels = c("5", "20", "26"))

p <- pheatmap(df, scale = "row", show_rownames = FALSE, 
              cluster_cols = F,
              treeheight_row = 30,
              annotation = anno, 
              show_colnames=F
              )


ggsave(p, width = 7, height = 5, filename = "C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/DSS results/methyldiff_heatmap_all.tiff")



# Method 2: use intersect() --------
# # data$start <- as.character(data$start)
# # dmlTest_BAPvsC_all$pos <- as.character(dmlTest_BAPvsC_all$pos)
# 
# Inter <- intersect(dmlTest_BAPvsC_all$pos, data$start)
# 
# geneID_BAPvsC_all <- data[which(data$start %in% Inter),]
# length(geneID_BAPvsC_all$start)
# 
# geneID_BAPvsC_all <- arrange(geneID_BAPvsC_all, start)
# 
# DML_BAPvsC_all <- dmlTest_BAPvsC_all[which(dmlTest_BAPvsC_all$pos %in% Inter),]
# DML_BAPvsC_all <- as.data.frame(DML_BAPvsC_all)
# length(DML_BAPvsC_all$pos)
# 
# DML_BAPvsC_all <- DML_BAPvsC_all[order(DML_BAPvsC_all$pos),]
# 
# DML_geneID_BAPvsC_all <- cbind(geneID_BAPvsC_all[,c( "seqnames","start","end","annotation","geneId")], DML_BAPvsC_all)
# 

# 2. callDMR
# DMRs are sorted by areaStat, which is defined in bsseq as the sum of the test statistics of all CpG sites within the DMR.
# Definition from bsseq
# areaStat:  areaStat which is the sum of the t-statistics in each CpG. This is kind of the area of the DMR, except that it is weighted by the number of CpGs and not by genomic length.
# This quantity does not have a direct biological meaning. 
# In ranking DMRs, it is not clear whether the length (number of CpG sites) or the height (methylation differences between two groups) is more important. areaStat is a combination of the two quantities. 
# It is an ad hoc way to rank DMRs, but larger areaStat is more likely to be a DMR.
# more info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8171293/

# Similarly, users can specify a threshold for difference. To filter by difference greater than 0.1, do:  delta = 0.1

# (1) BAP+TPA vs C
# filter by p value of 0.05
dmrs_BAPvsC_all = callDMR(dmlTest_BAPvsC_all, p.threshold = 0.05)
dmrs_BAPvsC_5w = callDMR(dmlTest_BAPvsC_5w, p.threshold = 0.05)
dmrs_BAPvsC_20w = callDMR(dmlTest_BAPvsC_20w, p.threshold = 0.05)
dmrs_BAPvsC_26w = callDMR(dmlTest_BAPvsC_26w, p.threshold = 0.05)

write.csv(dmrs_BAPvsC_all, "DMR_result_BAPvsC_all_p0.05.csv")
write.csv(dmrs_BAPvsC_5w, "DMR_result_BAPvsC_5w_p0.05.csv")
write.csv(dmrs_BAPvsC_20w, "DMR_result_BAPvsC_20w_p0.05.csv")
write.csv(dmrs_BAPvsC_26w, "DMR_result_BAPvsC_26w_p0.05.csv")

# filter by p value of 0.05 & methylation diff of 0.1
dmrs_BAPvsC_all_d0.1 = callDMR(dmlTest_BAPvsC_all, delta=0.1, p.threshold=0.05)
dmrs_BAPvsC_5w_d0.1 = callDMR(dmlTest_BAPvsC_5w, delta=0.1, p.threshold=0.05)
dmrs_BAPvsC_20w_d0.1 = callDMR(dmlTest_BAPvsC_20w, delta=0.1, p.threshold=0.05)
dmrs_BAPvsC_26w_d0.1 = callDMR(dmlTest_BAPvsC_26w, delta=0.1, p.threshold=0.05)

write.csv(dmrs_BAPvsC_all_d0.1, "DMR_result_BAPvsC_all_p0.05_d0.1.csv")
write.csv(dmrs_BAPvsC_5w_d0.1, "DMR_result_BAPvsC_5w_p0.05_d0.1.csv")
write.csv(dmrs_BAPvsC_20w_d0.1, "DMR_result_BAPvsC_20w_p0.05_d0.1.csv")
write.csv(dmrs_BAPvsC_26w_d0.1, "DMR_result_BAPvsC_26w_p0.05_d0.1.csv")
# (2) UA vs BAP+TPA
# filter by p value of 0.05
dmrs_UAvsBAP_all = callDMR(dmlTest_UAvsBAP_all, p.threshold = 0.05)
dmrs_UAvsBAP_5w = callDMR(dmlTest_UAvsBAP_5w, p.threshold = 0.05)
dmrs_UAvsBAP_20w = callDMR(dmlTest_UAvsBAP_20w, p.threshold = 0.05)
dmrs_UAvsBAP_26w = callDMR(dmlTest_UAvsBAP_26w, p.threshold = 0.05)

write.csv(dmrs_UAvsBAP_all, "DMR_result_UAvsBAP_all_p0.05.csv")
write.csv(dmrs_UAvsBAP_5w, "DMR_result_UAvsBAP_5w_p0.05.csv")
write.csv(dmrs_UAvsBAP_20w, "DMR_result_UAvsBAP_20w_p0.05.csv")
write.csv(dmrs_UAvsBAP_26w, "DMR_result_UAvsBAP_26w_p0.05.csv")


# filter by p value of 0.05 & methylation diff of 0.1
dmrs_UAvsBAP_all_d0.1 = callDMR(dmlTest_UAvsBAP_all, delta=0.1, p.threshold=0.05)
dmrs_UAvsBAP_5w_d0.1 = callDMR(dmlTest_UAvsBAP_5w, delta=0.1, p.threshold=0.05)
dmrs_UAvsBAP_20w_d0.1 = callDMR(dmlTest_UAvsBAP_20w, delta=0.1, p.threshold=0.05)
dmrs_UAvsBAP_26w_d0.1 = callDMR(dmlTest_UAvsBAP_26w, delta=0.1, p.threshold=0.05)

write.csv(dmrs_UAvsBAP_all_d0.1, "DMR_result_UAvsBAP_all_p0.05_d0.1.csv")
write.csv(dmrs_UAvsBAP_5w_d0.1, "DMR_result_UAvsBAP_5w_p0.05_d0.1.csv")
write.csv(dmrs_UAvsBAP_20w_d0.1, "DMR_result_UAvsBAP_20w_p0.05_d0.1.csv")
write.csv(dmrs_UAvsBAP_26w_d0.1, "DMR_result_UAvsBAP_26w_p0.05_d0.1.csv")


# ------Explanation of parameters/values in DMR result --------
# A data frame for DMRs. Each row is for a DMR. Rows are sorted by "areaStat", which is the sum
# of test statistics of all CpG sites in the region. The columns are:
# chr: Chromosome number.
# start, end: Genomic coordinates.
# length: Length of the DMR, in bps.
# nCG: Number of CpG sites contained in the DMR.
# meanMethy1, meanMethy2: Average methylation levels in two conditions.
# diff.Methy: The difference in the methylation levels between two conditions. 
  ##diff.Methy=meanMethy1-meanMethy2. (be careful!!)
# areaStat" The sum of the test statistics of all CpG sites within the DMR.


# -----Visualize DMRs using showOneDMR function ------
# This function provides more information than the plotRegion function in bsseq. 
# It plots the methylation percentages as well as the coverage depths at each CpG sites, instead of just the smoothed curve. So the coverage depth information will be available in the figure.
# more info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8171293/
# It displays the methylation levels as well as sequencing coverage depth information at each CpG site. The pink-shaded area is the detected DMR.
# Coverage (or depth) in DNA sequencing is the number of unique reads that include a given nucleotide in the reconstructed sequence.

par("mar")
par(mar=c(0.1,0.1,0.1,0.1))
showOneDMR(dmrs_BAPvsC_all_d0.1[1,], Exp)
colnames(Exp)

# check a specific DMR methylation % and the number of unique reads (one DMR at a time)

showOneDMR(dmrs_BAPvsC_5w[1,], Exp[, c(1:3, 7:9)])
showOneDMR(dmrs_BAPvsC_20w[1,], Exp[, c(4:6,10:12)])
showOneDMR(dmrs_BAPvsC_26w[1,], Exp[, c(4:6, 13:14)])

showOneDMR(dmrs_UAvsBAP_all[1,], Exp)
showOneDMR(dmrs_UAvsBAP_5w[1,], Exp[, c(7:9, 15:17)])
showOneDMR(dmrs_UAvsBAP_20w[1,], Exp[, c(10:12, 18:20)])
showOneDMR(dmrs_UAvsBAP_26w[1,], Exp[, c(13:14, 21:22)])


#---------- Analysis for data with multifactor design----------

# for results from multifactor design, {delta} is NOT supported. 
# This is because in multifactor design, the estimated coefficients in the regression are based on a GLM framework (loosely speaking), thus they don’t have clear meaning of methylation level differences. 
# So when the input DMLresult is from {DMLtest.multiFactor}, {delta} cannot be specified.

# Fit a linear model using DMLfit.multiFactor function, include trt, wk, and trt by wk interaction. 
# Similar to in a multiple regression, the model only needs to be fit once, and then the parameters can be tested based on the model fitting results.

## for combined data set ##
## (1) For interaction model, use:-------
DMLfit_int = DMLfit.multiFactor(Exp, design = design, formula = ~trt+wk+trt:wk)
colnames(DMLfit_int$X)


## (2) linear model without interaction , use:----------

DMLfit = DMLfit.multiFactor(Exp, design = design, formula = ~trt+wk)
head(DMLfit)

colnames(DMLfit$X)
DMLtest.trt = DMLtest.multiFactor(DMLfit, coef="wk")
ix=sort(DMLtest.trt[,"pvals"], index.return=TRUE)$ix
head(DMLtest.trt[ix,])
callDMR(DMLtest.trt, p.threshold=0.05)

## 5wks  ##

Exp_5w <- Exp[, c(1:3, 7:9)]
# (1) For interaction model, use:
# DMLfit_int = DMLfit.multiFactor(Exp_5w, design = design, formula = ~trt+wk+trt:wk)
# colnames(DMLfit_int$X)

# (2) linear model without interaction , use:
DMLfit = DMLfit.multiFactor(Exp_5w, design = design, formula = ~trt+wk)
head(DMLfit)

colnames(DMLfit$X)
DMLtest.trt = DMLtest.multiFactor(DMLfit, coef="wk")
ix=sort(DMLtest.trt[,"pvals"], index.return=TRUE)$ix
head(DMLtest.trt[ix,])
callDMR(DMLtest.cell, p.threshold=0.05)


# convert methylation level FC between comparison 1 and 2 to log2FC 
c1 <- read.table(file = "C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/DSS results/annotated/DSS_result_BAPvsC_5w_anno.csv", sep = ",", header = TRUE)
c2 <- read.table(file = "C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/DSS results/annotated/DSS_result_UAvsBAP_5w_anno.csv", sep = ",", header = TRUE)

c1.copy <- c1
c2.copy <- c2

colnames(c1)
c1 <- c1[, c("geneId", "start", "mu1", "mu2","pval", "fdr")]
c2 <- c2[, c("geneId", "start","mu1", "mu2","pval", "fdr")]

df <- merge(c1, c2, by = "start", all.y = F, all.x = F)

# df <- df %>% filter(fdr.x <=0.1 & fdr.y <= 0.1)

# remove the rows with NAs in mu1 and mu2 columns 
colnames(data)
df <- df[complete.cases(df[ , c(3,4,8,9)]), ]

# double check if gene names from two datasets matches and remove the gene name from second dataset
identical(df[,c("geneId.x")], df[,c("geneId.y")])
df <- as.data.frame(df)
df <- df[, -7]

# calculate the log2FC
# the methylation is inversely related to gene expression so we transformed log2FC by adding a minus sign

df$FC_c1 <- -log2(df$mu2.x/df$mu1.x)
df$FC_c2 <- -log2(df$mu2.y/df$mu1.y)

# rearrange order and save the file

df <- df[,c(2,1, 3, 4, 11,5, 6,7,8,12,9,10)]
write.csv(df, 
          file = "C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/DSS results/annotated/methyl_minus log2FC_5w-2.csv",
          row.names = FALSE)

# -----plot venn diagram ---------

setwd("C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/DSS results/annotated/")

res_TPA.vs.C_5wk <- read.csv("DSS_result_BAPvsC_5w_anno_fdr0.1.csv", header = T, sep = ",")
res_TPA.vs.C_20wk <- read.csv("DSS_result_BAPvsC_20w_anno_fdr0.1.csv", header = T, sep = ",")
res_TPA.vs.C_26wk <- read.csv("DSS_result_BAPvsC_26w_anno_fdr0.1.csv", header = T, sep = ",")

res_TPA.vs.C_5wk$diff <- -res_TPA.vs.C_5wk$diff
res_TPA.vs.C_20wk$diff <- -res_TPA.vs.C_20wk$diff
res_TPA.vs.C_26wk$diff <- -res_TPA.vs.C_26wk$diff

res_UA.vs.TPA_5wk <- read.csv("DSS_result_UAvsBAP_5w_anno_fdr0.1.csv", header = T, sep = ",")
res_UA.vs.TPA_20wk <- read.csv("DSS_result_UAvsBAP_20w_anno_fdr0.1.csv", header = T, sep = ",")
res_UA.vs.TPA_26wk <- read.csv("DSS_result_UAvsBAP_26w_anno_fdr0.1.csv", header = T, sep = ",")

res_UA.vs.TPA_5wk$diff <- -res_UA.vs.TPA_5wk$diff
res_UA.vs.TPA_20wk$diff <- -res_UA.vs.TPA_20wk$diff
res_UA.vs.TPA_26wk$diff <- -res_UA.vs.TPA_26wk$diff

# (1-1) Data preparation by geneId
## single geneId mapped to multipl DMRs so not working
# res_TPA.vs.C_5wk_up <- subset(res_TPA.vs.C_5wk, diff > 0.1) [, c("geneId")]
# res_TPA.vs.C_5wk_down <- subset(res_TPA.vs.C_5wk, diff < -0.1)[, c("geneId")]
# res_TPA.vs.C_20wk_up <- subset(res_TPA.vs.C_20wk,  diff > 0.1)[, c("geneId")]
# res_TPA.vs.C_20wk_down <- subset(res_TPA.vs.C_20wk, diff < -0.1)[, c("geneId")]
# res_TPA.vs.C_26wk_up <- subset(res_TPA.vs.C_26wk,  diff > 0.1)[, c("geneId")]
# res_TPA.vs.C_26wk_down <- subset(res_TPA.vs.C_26wk, diff < -0.1)[, c("geneId")]
# 
# res_UA.vs.TPA_5wk_up <- subset(res_UA.vs.TPA_5wk,  diff > 0.1)[, c("geneId")]
# res_UA.vs.TPA_5wk_down <- subset(res_UA.vs.TPA_5wk, diff < -0.1)[, c("geneId")]
# res_UA.vs.TPA_20wk_up <- subset(res_UA.vs.TPA_20wk,  diff > 0.1)[, c("geneId")]
# res_UA.vs.TPA_20wk_down <- subset(res_UA.vs.TPA_20wk,diff < -0.1)[, c("geneId")]
# res_UA.vs.TPA_26wk_up <- subset(res_UA.vs.TPA_26wk,  diff > 0.1)[, c("geneId")]
# res_UA.vs.TPA_26wk_down <- subset(res_UA.vs.TPA_26wk, diff < -0.1)[, c("geneId")]

## (1-2) Data preparation by start position
res_TPA.vs.C_5wk_up <- subset(res_TPA.vs.C_5wk, diff > 0.1) [, c("start")]
res_TPA.vs.C_5wk_down <- subset(res_TPA.vs.C_5wk, diff < -0.1)[, c("start")]
res_TPA.vs.C_20wk_up <- subset(res_TPA.vs.C_20wk,  diff > 0.1)[, c("start")]
res_TPA.vs.C_20wk_down <- subset(res_TPA.vs.C_20wk, diff < -0.1)[, c("start")]
res_TPA.vs.C_26wk_up <- subset(res_TPA.vs.C_26wk,  diff > 0.1)[, c("start")]
res_TPA.vs.C_26wk_down <- subset(res_TPA.vs.C_26wk, diff < -0.1)[, c("start")]

res_UA.vs.TPA_5wk_up <- subset(res_UA.vs.TPA_5wk,  diff > 0.1)[, c("start")]
res_UA.vs.TPA_5wk_down <- subset(res_UA.vs.TPA_5wk, diff < -0.1)[, c("start")]
res_UA.vs.TPA_20wk_up <- subset(res_UA.vs.TPA_20wk,  diff > 0.1)[, c("start")]
res_UA.vs.TPA_20wk_down <- subset(res_UA.vs.TPA_20wk,diff < -0.1)[, c("start")]
res_UA.vs.TPA_26wk_up <- subset(res_UA.vs.TPA_26wk,  diff > 0.1)[, c("start")]
res_UA.vs.TPA_26wk_down <- subset(res_UA.vs.TPA_26wk, diff < -0.1)[, c("start")]

## (2) Visualize in Venn Diagram
ggvenn(
  list("BAP+TPA vs. Control" = res_TPA.vs.C_5wk_up,
       "BAP+TPA+UA vs. BAP+TPA" = res_UA.vs.TPA_5wk_down),
  fill_color = c("pink", "light green"),
  stroke_size = 0, set_name_size = 8, text_size = 8) +
  ggtitle("5-wk Trt") +
  theme(plot.title = element_text(size = 28, hjust = 0.5,vjust = 2))
ggsave(file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/venn/venn_5wk_UAdown.tiff")

ggvenn(
  list("BAP+TPA vs. Control" = res_TPA.vs.C_5wk_down,
       "BAP+TPA+UA vs. BAP+TPA" = res_UA.vs.TPA_5wk_up),
  fill_color = c("light green", "pink"),
  stroke_size = 0, set_name_size = 8, text_size = 8) +
  ggtitle("5-wk Trt") +
  theme(plot.title = element_text(size = 28, hjust = 0.5,vjust = 2))
ggsave(file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/venn/venn_5wk_UAup.tiff")


#20wks#
ggvenn(
  list("BAP+TPA vs. Control" = res_TPA.vs.C_20wk_up,
       "BAP+TPA+UA vs. BAP+TPA" = res_UA.vs.TPA_20wk_down),
  fill_color = c("pink", "light green"),
  stroke_size = 0, set_name_size = 8, text_size = 8) +
  ggtitle("20-wk Trt") +
  theme(plot.title = element_text(size = 28, hjust = 0.5,vjust = 2))
ggsave(file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/venn/venn_20wk_UAdown.tiff")
ggvenn(
  list("BAP+TPA vs. Control" = res_TPA.vs.C_20wk_down,
       "BAP+TPA+UA vs. BAP+TPA" = res_UA.vs.TPA_20wk_up),
  fill_color = c("light green", "pink"),
  stroke_size = 0, set_name_size = 8, text_size = 8) +
  ggtitle("20-wk Trt") +
  theme(plot.title = element_text(size = 28, hjust = 0.5,vjust = 2))
ggsave(file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/venn/venn_20wk_UAup.tiff")


#26wks#
ggvenn(
  list("BAP+TPA vs. Control" = res_TPA.vs.C_26wk_up,
       "BAP+TPA+UA vs. BAP+TPA" = res_UA.vs.TPA_26wk_down),
  fill_color = c("pink", "light green"),
  stroke_size = 0, set_name_size = 8, text_size = 8) +
  ggtitle("26-wk Trt") +
  theme(plot.title = element_text(size = 28, hjust = 0.5,vjust = 2))
ggsave(file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/venn/venn_26wk_UAdown.tiff")

ggvenn(
  list("BAP+TPA vs. Control" = res_TPA.vs.C_26wk_down,
       "BAP+TPA+UA vs. BAP+TPA" = res_UA.vs.TPA_26wk_up),
  fill_color = c("light green", "pink"),
  stroke_size = 0, set_name_size = 8, text_size = 8) +
  ggtitle("26-wk Trt") +
  theme(plot.title = element_text(size = 28, hjust = 0.5,vjust = 2))
ggsave(file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/venn/venn_26wk_UAup.tiff")



# ------- plot heatmap of genes with inverse diff with filters------- 
## For two comparison groups (TPAvsC versus UAvsTPA) 
## (1) data has been pre-filtered by padj <= 0.1
## (2) find the inverse correlated subset of gens
## 5wk -UA down
colnames(res_TPA.vs.C_5wk)
subset_5w_TPAup <- subset(res_TPA.vs.C_5wk, diff > 0.1) [, c("geneId", "start", "diff")]
subset_5w_UAdown <- subset(res_UA.vs.TPA_5wk, diff < -0.1) [, c("geneId", "start", "diff")]
DEG_5w_TPAup_UAdown <- merge(subset_5w_TPAup, subset_5w_UAdown, by ="start")

## 5wk -UA up
subset_5w_TPAdown <- subset(res_TPA.vs.C_5wk, diff < -0.1) [, c("geneId", "start", "diff")]
subset_5w_UAup <- subset(res_UA.vs.TPA_5wk, diff > 0.1) [, c("geneId", "start", "diff")]
DEG_5w_TPAdown_UAup <- merge(subset_5w_TPAdown, subset_5w_UAup, by ="start")

## combined above gene lists
head(DEG_5w)
DEG_5w <- rbind(DEG_5w_TPAup_UAdown, DEG_5w_TPAdown_UAup)
DEG_5w <- DEG_5w[, -c(1,4)]
DEG_5w <- DEG_5w[order(abs(DEG_5w[,2] -DEG_5w[,3]), decreasing = T), ]

u <- match(unique(DEG_5w$geneId.x),DEG_5w$geneId.x)

DEG_5w <- DEG_5w[u,]
DEG_5w_20 <- head(DEG_5w[order(abs(DEG_5w[,2] -DEG_5w[,3]), decreasing = T), ], 20)
rownames(DEG_5w_20) <- DEG_5w_20[,1]
DEG_5w_20 <- DEG_5w_20[,-1]
DEG_5w_20[,1] <- DEG_5w_20[,1]*100
DEG_5w_20[,2] <- DEG_5w_20[,2]*100

## 20wk -UA down

subset_20w_TPAup <- subset(res_TPA.vs.C_20wk, diff > 0.1) [, c("geneId", "start", "diff")]
subset_20w_UAdown <- subset(res_UA.vs.TPA_20wk, diff < -0.1) [, c("geneId", "start", "diff")]
DEG_20w_TPAup_UAdown <- merge(subset_20w_TPAup, subset_20w_UAdown, by ="start")

## 20wk -UA up
subset_20w_TPAdown <- subset(res_TPA.vs.C_20wk, diff < -0.1) [, c("geneId", "start", "diff")]
subset_20w_UAup <- subset(res_UA.vs.TPA_20wk, diff > 0.1) [, c("geneId", "start", "diff")]
DEG_20w_TPAdown_UAup <- merge(subset_20w_TPAdown, subset_20w_UAup, by ="start")

## combined above gene lists
DEG_20w <- rbind(DEG_20w_TPAup_UAdown, DEG_20w_TPAdown_UAup)
DEG_20w <- DEG_20w[, -c(1,4)]
DEG_20w <- DEG_20w[order(abs(DEG_20w[,2] -DEG_20w[,3]), decreasing = T), ]

u <- match(unique(DEG_20w$geneId.x),DEG_20w$geneId.x)

DEG_20w <- DEG_20w[u,]
DEG_20w_20 <- head(DEG_20w[order(abs(DEG_20w[,2] -DEG_20w[,3]), decreasing = T), ], 20)
rownames(DEG_20w_20) <- DEG_20w_20[,1]
DEG_20w_20 <- DEG_20w_20[,-1]
DEG_20w_20[,1] <- DEG_20w_20[,1]*100
DEG_20w_20[,2] <- DEG_20w_20[,2]*100

## 26wk -UA down

subset_26w_TPAup <- subset(res_TPA.vs.C_26wk, diff > 0.1) [, c("geneId", "start", "diff")]
subset_26w_UAdown <- subset(res_UA.vs.TPA_26wk, diff < -0.1) [, c("geneId", "start", "diff")]
DEG_26w_TPAup_UAdown <- merge(subset_26w_TPAup, subset_26w_UAdown, by ="start")

## 26wk -UA up
subset_26w_TPAdown <- subset(res_TPA.vs.C_26wk, diff < -0.1) [, c("geneId", "start", "diff")]
subset_26w_UAup <- subset(res_UA.vs.TPA_26wk, diff > 0.1) [, c("geneId", "start", "diff")]
DEG_26w_TPAdown_UAup <- merge(subset_26w_TPAdown, subset_26w_UAup, by ="start")

## combined above gene lists
DEG_26w <- rbind(DEG_26w_TPAup_UAdown, DEG_26w_TPAdown_UAup)
DEG_26w <- DEG_26w[, -c(1,4)]
DEG_26w <- DEG_26w[order(abs(DEG_26w[,2] -DEG_26w[,3]), decreasing = T), ]

u <- match(unique(DEG_26w$geneId.x),DEG_26w$geneId.x)

DEG_26w <- DEG_26w[u,]
DEG_26w_20 <- head(DEG_26w[order(abs(DEG_26w[,2] -DEG_26w[,3]), decreasing = T), ], 20)
rownames(DEG_26w_20) <- DEG_26w_20[,1]
DEG_26w_20 <- DEG_26w_20[,-1]
DEG_26w_20[,1] <- DEG_26w_20[,1]*100
DEG_26w_20[,2] <- DEG_26w_20[,2]*100

# (3) visualize in heatmaps

trt <- c("BAP+TPA vs. C", "BAP+TPA+UA vs. BAP+TPA")
trt <- as.data.frame(trt)
trt
colnames(DEG_5w_20) <- c("BAP+TPA vs. C", "BAP+TPA+UA vs. BAP+TPA")
colnames(DEG_20w_20) <- c("BAP+TPA vs. C", "BAP+TPA+UA vs. BAP+TPA")
colnames(DEG_26w_20) <- c("BAP+TPA vs. C", "BAP+TPA+UA vs. BAP+TPA")

rownames(trt) <- c("BAP+TPA vs. C", "BAP+TPA+UA vs. BAP+TPA")

library (RColorBrewer)
display.brewer.all()

breakslist <- seq(-50, 50, length.out = 1000)
DMR_5w_heatmap <- pheatmap(DEG_5w_20, cluster_cols = F, treeheight_row = 0,
                           show_rownames = T, show_colnames=F, 
                           main = "5-wk", fontsize = 6, fontsize_row = 7,
                           cellwidth=20, cellheight = 20, annotation_col = trt,
                           color =  colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(length(breakslist)),
                           breaks = breakslist)

DMR_20w_heatmap <- pheatmap(DEG_20w_20, cluster_cols = F, treeheight_row = 0,
                           show_rownames = T, show_colnames=F, 
                           main = "20-wk", fontsize = 6, fontsize_row = 7,
                           cellwidth=20, cellheight = 20, annotation_col = trt,
                           color =  colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(length(breakslist)),
                           breaks = breakslist)
DMR_26w_heatmap <- pheatmap(DEG_26w_20, cluster_cols = F, treeheight_row = 0,
                           show_rownames = T, show_colnames=F, 
                           main = "26-wk", fontsize = 6, fontsize_row = 7,
                           cellwidth=20, cellheight = 20, annotation_col = trt,
                           color =  colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(length(breakslist)),
                           breaks = breakslist)

# (4) save plots
ggsave(width = 3.5, height = 6.4, DMR_5w_heatmap, file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/heatmap/Top20DMR_5w_heatmap.tiff")
ggsave(width = 3.5, height = 6.4, DMR_20w_heatmap, file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/heatmap/Top20DMR_20w_heatmap.tiff")
ggsave(width = 3.5, height = 6.4, DMR_26w_heatmap, file="C:/Rutgers/Github/Methyl-seq-Analysis-Workflow_BAP-UA-skin-project/heatmap/Top20DMR_26w_heatmap.tiff")
