setwd("H:/Bradley/cellrangerGRCh38LGL7flows/Seurat/phenograph")
##Because the data is too large, I used the biowulf to calculate see scaleData.R. The result file is around 36GB with colSums>10000
rm(list=ls())
load("scaleDataAll.RData")
library(data.table)
dataOutAll<-fread("dataOutAll_meta_expression_tSNE_cluster.csv", data.table=FALSE)
rownames(dataOutAll)<-dataOutAll[,1]
dataOutAll<-dataOutAll[,-1]
#colnames(dataOutAll)
dim(scaledDataAll); dim(dataOutAll);length(intersect(rownames(dataOutAll), colnames(scaledDataAll)));
dataOutput<-matrix(1, nrow=dim(scaledDataAll)[1], ncol=32*6)
rownames(dataOutput)<-rownames(scaledDataAll)
colnames(dataOutput)<-c(paste("S", 1:32, sep="_"), paste("S_sd", 1:32, sep="_"), paste("SCD8", 1:32, sep="_")
, paste("SCD8_sd", 1:32, sep="_"), paste("SCD4", 1:32, sep="_"), paste("SCD4_sd", 1:32, sep="_"))
for(sample in 1:32){
  cat(sample); cat("\n")
  cell<-rownames(dataOutAll[dataOutAll[, "orig.ident"]==paste("S", sample, sep=""), ])
  dataOutput[, paste("S", sample, sep="_")]<-apply(scaledDataAll[, cell], 1, mean)
  dataOutput[, paste("S_sd", sample, sep="_")]<-apply(scaledDataAll[, cell], 1, sd)
  #CD8
  cell<-rownames(dataOutAll[dataOutAll[, "orig.ident"]==paste("S", sample, sep="") & dataOutAll[, "cluster_phenographmeta_combat"]%in%c(0, 8, 4,5,6), ])
  dataOutput[, paste("SCD8", sample, sep="_")]<-apply(scaledDataAll[, cell], 1, mean)
  dataOutput[, paste("SCD8_sd", sample, sep="_")]<-apply(scaledDataAll[, cell], 1, sd)
  #CD4
  cell<-rownames(dataOutAll[dataOutAll[, "orig.ident"]==paste("S", sample, sep="") & dataOutAll[, "cluster_phenographmeta_combat"]%in%c(1,2,3,7,9), ])
  dataOutput[, paste("SCD4", sample, sep="_")]<-apply(scaledDataAll[, cell], 1, mean)
  dataOutput[, paste("SCD4_sd", sample, sep="_")]<-apply(scaledDataAll[, cell], 1, sd)
}
write.csv(dataOutput, file="dataOutputAverageSample_AliceCD4CD8.csv")