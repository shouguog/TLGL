setwd("H:/Bradley/Integration/DetailedAnalysis_round1/clonaltypesummary")
library(stringr)
#Let us filter combine the data of 32
rm(list=ls())
patientInfo<-read.csv("../../../180418LGLExpList.csv", row.names = 1)

dataOutput<-matrix(1, nrow=32, ncol=11)
rownames(dataOutput)<-rownames(patientInfo)
colnames(dataOutput)<-c("allCells", "withA", "withB", "withAB", "with1A", "with2A","with3A", "with1B", "with2B","with3B","with1A1B")

for(sample in rownames(patientInfo)){
  ###Let us read the data
  dataFiltered1<-read.csv(file=paste("../../../cellrangerVDJ/SNOLANE1_",gsub("S", "", sample), "_VDJT/outs/clonotypes.csv", sep=""), header=T)
  dataOutput[sample,1]<-sum(dataFiltered1$frequency)
  dataFilteredTRA<-dataFiltered1[grepl("TRA:", dataFiltered1$cdr3s_nt),]
  dataOutput[sample,2]<-sum(dataFilteredTRA$frequency)
  dataFilteredTRB<-dataFiltered1[grepl("TRB:", dataFiltered1$cdr3s_nt),]
  dataOutput[sample,3]<-sum(dataFilteredTRB$frequency)
  dataFilteredTRAB<-dataFiltered1[grepl("TRA:", dataFiltered1$cdr3s_nt) & grepl("TRB:", dataFiltered1$cdr3s_nt),]
  dataOutput[sample,4]<-sum(dataFilteredTRAB$frequency)
  
  dataFiltereduniqueTRA<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRA:")==1,]
  dataOutput[sample,5]<-sum(dataFiltereduniqueTRA$frequency)
  dataFiltereddoubleTRA<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRA:")==2,]
  dataOutput[sample,6]<-sum(dataFiltereddoubleTRA$frequency)
  dataFilteredtripleTRA<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRA:")==3,]
  dataOutput[sample,7]<-sum(dataFilteredtripleTRA$frequency)
  
  dataFiltereduniqueTRB<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRB:")==1,]
  dataOutput[sample,8]<-sum(dataFiltereduniqueTRB$frequency)
  dataFiltereddoubleTRB<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRB:")==2,]
  dataOutput[sample,9]<-sum(dataFiltereddoubleTRB$frequency)
  dataFilteredtripleTRB<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRB:")==3,]
  dataOutput[sample,10]<-sum(dataFilteredtripleTRB$frequency)

  dataFiltereduniqueTRAB<-dataFiltered1[str_count(dataFiltered1$cdr3s_nt, "TRA:")==1 & str_count(dataFiltered1$cdr3s_nt, "TRB:")==1,]
  dataOutput[sample,11]<-sum(dataFiltereduniqueTRAB$frequency)
  write.csv(dataFiltereduniqueTRAB, paste("paired1A1Bdata/", sample, ".csv", sep=""))
}
write.csv(cbind(patientInfo,dataOutput), "clonoCapture.csv")
