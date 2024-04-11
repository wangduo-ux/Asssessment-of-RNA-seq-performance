
library(data.table)


### function for SNR calculation  ###

snrdb_function<-function(exprMat,group){
  
  library(data.table)
  
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),retx=T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  
  signoise_db <- 10*log10(signoise)
  return(signoise_db)
  
}


### calculating SNR17 ###

snr17_calculation <- function(expr){
  expr <- log2(expr + 0.01)
  expr <- expr[!apply(expr, 1, function(x) length(unique(x))) == 1, ]
  sample_group <- data.frame(
    sample = c("M8", "M8", "M8", "F7", "F7", "F7", "D5", "D5", "D5", "T1", "T1", "T1", "T2", "T2", "T2", "D6", "D6","D6"),
    library = c("M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3", "D5_1", 'D5_2', 'D5_3', "T1_1", "T1_2", 'T1_3', "T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3")
  )
  snr17 <- {}
  for (j in 1:18){
    
    expr <- expr[, -j]
    sample_group <- sample_group[-j, ]
    expr <- log2(expr + 0.01)
    dat_z <- t(apply(expr,1,function(x){(x-mean(x))/(sd(x))}))  ### expression matrix should be z-scaled before using snrdb_function ###
    sample <- sample_group$sample[match(colnames(dat_z),sample_group$library)]
    SNR <- round(snrdb_function(dat_z,as.factor(sample)),1)
    
  }
  snr17 <- append(snr17, SNR)
  return(snr17)
}



### claculating SNR18 ###
snr18_calculation <- function(expr){
  expr <- log2(expr + 0.01)
  expr <- expr[!apply(expr, 1, function(x) length(unique(x))) == 1, ]
  sample_group <- data.frame(
    sample = c("M8", "M8", "M8", "F7", "F7", "F7", "D5", "D5", "D5", "T1", "T1", "T1", "T2", "T2", "T2", "D6", "D6","D6"),
    library = c("M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3", "D5_1", 'D5_2', 'D5_3', "T1_1", "T1_2", 'T1_3', "T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3")
  )
  
  dat_z <- t(apply(expr,1,function(x){(x-mean(x))/(sd(x))}))  ### expression matrix should be z-scaled before using snrdb_function ###

  sample<-sample_group$sample[match(colnames(dat_z),sample_group$library)]

  SNR18 <- round(snrdb_function(dat_z,as.factor(sample)),1)  ### this will output a SNR of 18 samples ###
  return(SNR18)
}



### start ###

### 1. Read an expression matrix containing 19 columns, the first column is gene names and the other 18 colomuns are samples, each column representing a sample and each row representing a gene ###
###    Column names must be "ensembl_id","M8_1", "M8_2", "M8_3", "F7_1", "F7_2", "F7_3","D5_1", "D5_2","D5_3", "T1_1", "T1_2","T1_3","T2_1", "T2_2", "T2_3", "D6_1", "D6_2","D6_3" ###
expr<-read.table("/your path/quartet.expr.csv",row.names=1,sep = '\t', header = TRUE,check.names=FALSE)

SNR18 <- snr18_calculation(expr) 
SNR17 <- snr17_calculation(expr)  ###  this will return a list including 18 SNR17 values, calculated from any 17 out of 18 samples ###

