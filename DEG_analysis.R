library(edgeR)
library(DESeq2)
library(limma)
library(DEGseq)
library(EBSeq)
library(dplyr)


########################################
############### edgeR ##################
########################################

edgeR_DEG <- function(count_data,x) {
  group <- factor(c(rep("test",3),rep("control",3)))
  lrt <- NULL
  keep <- rowSums(count_data) >= x
  count_data <- count_data[keep, ]
  y <- DGEList(counts = count_data, genes = rownames(count_data), group = group)
  y <- estimateDisp(y)
  y <- calcNormFactors(y)
  et <- exactTest(y)
  DEG <- summary(de <- decideTestsDGE(et, adjust.method = "fdr", p.value = 0.05, lfc = 1))
  down <- DEG[1, ]
  NotDEG <- DEG[2, ]
  up <- DEG[3, ]
  lrt <- topTags(et, n = nrow(y$counts))
  lrt <- as.data.frame(lrt)
  lrt$genes <- sapply(strsplit(lrt$genes, split = "|", fixed = TRUE), "[", 1)
  lrt[which(lrt$logFC >= 1 & lrt$FDR < 0.05), 'sig'] <- 'Upregulate'
  lrt[which(lrt$logFC <= -1 & lrt$FDR < 0.05), 'sig'] <- 'Downregulate'
  lrt[which(abs(lrt$logFC) < 1 | lrt$FDR >= 0.05), 'sig'] <- 'NonDE'
  colnames(lrt)[colnames(lrt) == "genes"] <- "ensembl_id"
  return(lrt)
}

########################################
############### limma ##################
########################################

limma_DEG <- function(count_data,x) {
  group <- factor(c(rep("test",3),rep("control",3)))
  design <- model.matrix(~group)
  colnames(design) <- levels(group)
  keep <- rowSums(count_data) >= x
  count_data <- count_data[keep, ]
  dge <- DGEList(count_data)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=F, normalize="quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  lrt <- topTable(fit, coef=ncol(design),n=Inf)
  up <- nrow(subset(lrt, adj.P.Val < 0.05 & logFC >= 1))
  down <- nrow(subset(lrt, adj.P.Val < 0.05 & logFC <= -1))
  row_names <- rownames(lrt)
  rownames(lrt) <- sapply(strsplit(row_names, split = "|", fixed = TRUE), "[", 1)
  lrt[which(lrt$logFC >= 1 & lrt$adj.P.Val < 0.05), 'sig'] <- 'Upregulate'
  lrt[which(lrt$logFC <= -1 & lrt$adj.P.Val < 0.05), 'sig'] <- 'Downregulate'
  lrt[which(abs(lrt$logFC) < 1 | lrt$adj.P.Val >= 0.05), 'sig'] <- 'NonDE'
  lrt$ensembl_id <- row.names(lrt)
  return(lrt)
}

########################################
############### DESeq2 #################
########################################

DEseq_DEG <- function(count_data,x) {
  condition <- factor(c(rep("test",3), rep("control",3)))
  count_data[,c(1:6)] <- as.data.frame(lapply(count_data[,c(1:6)], as.integer))
  coldata <- data.frame(row.names = colnames(count_data), condition)
  dds <- DESeqDataSetFromMatrix(countData = count_data,colData = coldata, design = ~ condition)
  keep <- rowSums(counts(dds)) >= x
  dds <- dds[keep,]
  n <- nrow(dds)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[!is.na(res$padj), ]
  res <- as.data.frame(res)
  row_names <- row.names(res)
  rownames(res) <- sapply(strsplit(row_names, split = "|", fixed = TRUE), "[", 1)
  res[which(res$log2FoldChange >= 1 & res$padj < 0.05), 'sig'] <- 'Upregulate'
  res[which(res$log2FoldChange <= -1 & res$padj < 0.05), 'sig'] <- 'Downregulate'
  res[which(abs(res$log2FoldChange) < 1 | res$padj >= 0.05), 'sig'] <- 'NonDE'
  res$ensembl_id <- row.names(res)
  return(res)
}

########################################
############### EBSeq  #################
########################################

EBseq_DEG <- function(count_data,x) {
  condition=as.factor(rep(c("test","control"),each=3))
  keep <- rowSums(count_data) >= x
  count_data <- count_data[keep,]
  n <- nrow(count_data)
  Sizes=MedianNorm(count_data)
  count_data <- as.matrix(count_data)
  EBOut = EBTest(count_data,Conditions = condition,sizeFactors=Sizes, maxround=5)
  EBDERes=GetDEResults(EBOut, FDR=0.05)
  GeneFC=PostFC(EBOut)
  result <- merge(EBDERes$PPMat,EBDERes$Status,by = "row.names")
  result <- merge(result,GeneFC,by.x = "Row.names",by.y = "row.names")
  DEG <- result[result$y == "DE" | result$y == "EE", ]
  DEG$logFC <- log2(DEG$RealFC)
  DEG[which(DEG$logFC >= 1 & DEG$y == "DE"), 'sig'] <- 'Upregulate'
  DEG[which(DEG$logFC <= -1 & DEG$y == "DE"), 'sig'] <- 'Downregulate'
  DEG[which(abs(DEG$logFC) < 1 | DEG$y == "EE"), 'sig'] <- 'NonDE'
  colnames(DEG)[colnames(DEG) == "Row.names"] <- "ensembl_id"
  return(DEG)
}

########################################
############### DEGseq #################
########################################

DEGseq_DEG <- function(count_data,x,tmp_path) {
  keep <- rowSums(count_data[, 2:7]) >= x
  count_data <- count_data[keep,]
  DEGexp(geneExpMatrix1=count_data, geneCol1=1, expCol1=c(2,3,4), groupLabel1="test", 
         geneExpMatrix2=count_data, geneCol2=1, expCol2=c(5,6,7), groupLabel2="control", 
         method="MARS", qValue=0.05, outputDir=tmp_path, thresholdKind=3)
  DEG <- read.table(paste0(tmp_path,"/output_score.txt"),sep="\t",header=1,check.names=FALSE)
  DEG <- DEG[, c("GeneNames","log2(Fold_change) normalized","q-value(Benjamini et al. 1995)")]
  colnames(DEG) <- c("GeneNames","logFC","q_value")
  DEG[which(DEG$logFC >= 1 & DEG$q_value < 0.05), 'sig'] <- 'Upregulate'
  DEG[which(DEG$logFC <= -1 & DEG$q_value < 0.05), 'sig'] <- 'Downregulate'
  DEG[which(abs(DEG$logFC) < 1 | DEG$q_value >= 0.05), 'sig'] <- 'NonDE'
  colnames(DEG)[colnames(DEG) == "GeneNames"] <- "ensembl_id"
  return(DEG)
}

###################################################################
###    calculating MCC based on the Quartet reference datasets  ###
###################################################################

MCC_calculation_1 <- function(DEG_data,truth_1) {
  truth_1 <- truth_1[, c("ensembl_id", "DEG")]
  df <- merge(truth_1,DEG_data,by = "ensembl_id")
  df$combined <- paste(df$DEG, df$sig, sep = "_")
  TP <- sum(df$combined == "Upregulate_Upregulate" | df$combined == "Downregulate_Downregulate")
  FP1 <- sum(df$combined == "Upregulate_Downregulate" | df$combined == "Downregulate_Upregulate")
  FP2 <- sum(df$combined == "NonDE_Downregulate" | df$combined == "NonDE_Upregulate")
  TN <- sum(df$combined == "NonDE_NonDE")
  TPR <- TP / (5036 - FP1)  #for pc gene, modified as 3337
  FN <- 5036 - TP - FP1
  precision <- TP / (TP + FP1 + FP2)
  FPR <- (FP1 + FP2) / (FP1 + FP2 + TN)
  mcc <- (TP * TN - (FP1 + FP2) * FN) / sqrt((TP + (FP1 + FP2)) * (TP + FN) * (TN + (FP1 + FP2)) * (TN + FN))
  F1 = 2 * (precision * TPR) / (precision + TPR)
  return(mcc)
}



###################################################################
###     calculating MCC based on the Quartet TaqMan datasets    ###
###################################################################

MCC_calculation_2 <- function(DEG_data,truth_2){
  truth_2 <- truth_2[, c("ensembl_id", "DEG")] 
  df <- merge(truth_2,DEG_data,by = "ensembl_id")
  df$combined <- paste(df$DEG, df$sig, sep = "_")
  TP <- sum(df$combined == "Upregulate_Upregulate" | df$combined == "Downregulate_Downregulate")
  FP1 <- sum(df$combined == "Upregulate_Downregulate" | df$combined == "Downregulate_Upregulate")
  FP2 <- sum(df$combined == "NonDE_Downregulate" | df$combined == "NonDE_Upregulate")
  TN <- sum(df$combined == "NonDE_NonDE")
  TPR <- TP/(197-FP1)
  precision <- TP/(TP+FP1+FP2)
  FPR <- (FP1+FP2)/(FP1+FP2 + TN)
  FN <- 197 - TP - FP1 
  mcc <- (TP * TN - (FP1+FP2) * FN) / sqrt((TP + (FP1+FP2)) * (TP + FN) * (TN + (FP1+FP2)) * (TN + FN))
  F1 = 2 * (precision * TPR) / (precision + TPR)
  return(mcc)
}

###################################################################
###      calculating MCC based on the MAQC  TaqMan datasets     ###
###################################################################

MCC_calculation_3 <- function(DEG_data,truth_3){
  truth_3 <- truth_3[, c("ensembl_id", "DEG")]  
  df <- merge(truth_3,DEG_data,by = "ensembl_id")
  df$combined <- paste(df$DEG, df$sig, sep = "_")
  TP <- sum(df$combined == "Upregulate_Upregulate" | df$combined == "Downregulate_Downregulate")
  FP1 <- sum(df$combined == "Upregulate_Downregulate" | df$combined == "Downregulate_Upregulate")
  FP2 <- sum(df$combined == "NonDE_Downregulate" | df$combined == "NonDE_Upregulate")
  TN <- sum(df$combined == "NonDE_NonDE")
  TPR <- TP/(536-FP1)
  precision <- TP/(TP+FP1+FP2)
  FPR <- (FP1+FP2)/(FP1+FP2 + TN)
  FN <- 536 - TP - FP1 
  mcc <- (TP * TN - (FP1+FP2) * FN) / sqrt((TP + (FP1+FP2)) * (TP + FN) * (TN + (FP1+FP2)) * (TN + FN))
  F1 = 2 * (precision * TPR) / (precision + TPR)
  return(mcc)
}

### For calculating three highest MCC at the same times, this will output a list, including MCC based on Quartet reference datasets, Quartet TaqMan datasets, and MAQC TaqMan datasets, respectively ###
function_list <- list(edgeR_DEG, limma_DEG, DEseq_DEG, EBseq_DEG, DEGseq_DEG)

MCC_all <- function(count_data,truth_1,truth_2,truth_3,i){
  MCC_1_list <- list()
  MCC_2_list <- list()
  MCC_3_list <- list()
  
  maqc <- count_data[,c("A_1","A_2","A_3","B_1","B_2","B_3")]
  m8_vs_d6 <- count_data[,c("M8_1", "M8_2", "M8_3", "D6_1", "D6_2","D6_3")]
  f7_vs_d6 <- count_data[,c("F7_1", "F7_2", "F7_3", "D6_1", "D6_2","D6_3")]
  d5_vs_d6 <- count_data[,c("D5_1", "D5_2", "D5_3", "D6_1", "D6_2","D6_3")]
  
  if (i == 5){
    raw_data$gene_id <- row.names(raw_data)
    maqc <- count_data[,c("gene_id","A_1","A_2","A_3","B_1","B_2","B_3")]
    m8_vs_d6 <- count_data[,c("gene_id","M8_1", "M8_2", "M8_3", "D6_1", "D6_2","D6_3")]
    f7_vs_d6 <- count_data[,c("gene_id","F7_1", "F7_2", "F7_3", "D6_1", "D6_2","D6_3")]
    d5_vs_d6 <- count_data[,c("gene_id","D5_1", "D5_2", "D5_3", "D6_1", "D6_2","D6_3")]
    
  }
  ### Using different thresholds to filter low-expression genes, then select the highest MCC value ###
  for (x in c(seq(0, 60, by = 1))) {
    maqc_DEG <- function_list[[i]](maqc,x)
    m8_DEG <- function_list[[i]](m8_vs_d6,x)
    f7_DEG <- function_list[[i]](f7_vs_d6,x)
    d5_DEG <- function_list[[i]](d5_vs_d6,x)
    if (i == 5){
      maqc_DEG <- function_list[[i]](maqc,x,tmp_path)
      m8_DEG <- function_list[[i]](m8_vs_d6,x,tmp_path)
      f7_DEG <- function_list[[i]](f7_vs_d6,x,tmp_path)
      d5_DEG <- function_list[[i]](d5_vs_d6,x,tmp_path)
    }
    m8_DEG$ensembl_id <- paste("M8_", m8_DEG$ensembl_id, sep = "")
    f7_DEG$ensembl_id <- paste("F7_", f7_DEG$ensembl_id, sep = "")
    d5_DEG$ensembl_id <- paste("D5_", d5_DEG$ensembl_id, sep = "")
    quartet_DEG <- bind_rows(m8_DEG, f7_DEG, d5_DEG)
    MCC_1 <- MCC_calculation_1(quartet_DEG,truth)
    MCC_2 <- MCC_calculation_2(quartet_DEG,truth_2)
    MCC_3 <- MCC_calculation_3(maqc_DEG,truth_3)
    cat(MCC_1,MCC_2,MCC_3)
    MCC_1_list <- append(MCC_1_list, MCC_1)
    MCC_2_list <- append(MCC_2_list, MCC_2)
    MCC_3_list <- append(MCC_3_list, MCC_3)
  }
  res <- c(max(unlist(MCC_1_list)), max(unlist(MCC_2_list)), max(unlist(MCC_3_list))) 
  return(res)
}





### start ###

###  1. Read an read count matrix containing 25 columns, the first column is gene names and other 18 columns are samples, each column representing a sample and each row representing a gene  ####
###     Column names must be "ensembl_id","A_1","A_2","A_3","B_1","B_2","B_3","M8_1", "M8_2", "M8_3", "D5_1", "D5_2","D5_3","F7_1", "F7_2", "F7_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###

count_data <- read.table("/your path/gene_count.csv", sep="\t", header = T, row.names = 1)

### 2. Read the Quartet reference datasets, this could be downloaded from https://github.com/wangduo-ux/Asssessment-of-RNA-seq-performance.git or https://chinese-quartet.org/#/dashboard ###

truth_1 <- read.table("/your path/quartet_reference_datasets.csv", sep="\t", header = T)

### 3. Read the Quartet TaqMan datasets, this could be downloaded from https://github.com/wangduo-ux/Asssessment-of-RNA-seq-performance.git ###

truth_2 <-read.table("/your path/quartet_taqman.csv", sep="\t", header = T)


### 4. Read the MAQC TaqMan datasets, this could be downloaded from https://github.com/wangduo-ux/Asssessment-of-RNA-seq-performance.git ###

truth_3 <-read.table("/your path/maqc_taqman.csv", sep="\t", header = T)

tmp_path <- "/your path/"    ### provide a path to save tmp file when running DEGseq ###

###  5. calculating optimal MCC after applying a series of thresholds for filtering low-expression genes  ###
MCC <- MCC_all(count_data,truth_1,truth_2,truth_3,i)     ### set a paramater i to select DEG analysis tools, "1" for edgeR, "2" for limma, "3" for DEseq2, "4" for EBSeq, and "5" for DEGseq  ###


