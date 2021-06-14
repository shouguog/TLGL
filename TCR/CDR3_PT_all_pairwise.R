setwd("H:/Bradley/Integration/DetailedAnalysis_round1/circlize")
rm(list=ls())
intTopNum<-2000
###Let us read the data
dataFiltered1<-read.csv(file="../../../cellrangerVDJ/SNOLANE1_1_VDJT/outs/filtered_contig_annotations.csv", header=T)
dataFiltered1<-dataFiltered1[dataFiltered1[, "raw_consensus_id"]!="None",]
dataFiltered2<-read.csv(file="../../../cellrangerVDJ/SNOLANE1_3_VDJT/outs/filtered_contig_annotations.csv", header=T)
dataFiltered2<-dataFiltered2[dataFiltered2[, "raw_consensus_id"]!="None",]
dataFiltered3<-read.csv(file="../../../cellrangerVDJ/SNOLANE1_5_VDJT/outs/filtered_contig_annotations.csv", header=T)
dataFiltered3<-dataFiltered3[dataFiltered3[, "raw_consensus_id"]!="None",]
dataFiltered4<-read.csv(file="../../../cellrangerVDJ/SNOLANE1_9_VDJT/outs/filtered_contig_annotations.csv", header=T)
dataFiltered4<-dataFiltered4[dataFiltered4[, "raw_consensus_id"]!="None",]
dataFiltered5<-read.csv(file="../../../cellrangerVDJ/SNOLANE1_11_VDJT/outs/filtered_contig_annotations.csv", header=T)
dataFiltered5<-dataFiltered5[dataFiltered5[, "raw_consensus_id"]!="None",]
dataFiltered6<-read.csv(file="../../../cellrangerVDJ/SNOLANE1_13_VDJT/outs/filtered_contig_annotations.csv", header=T)
dataFiltered6<-dataFiltered6[dataFiltered6[, "raw_consensus_id"]!="None",]

####check CDR3
dataFiltered_CDR3<-list()
dataFiltered1_CDR3<-table(dataFiltered1[, "cdr3_nt"])
dataFiltered_CDR3[[1]]<-names(sort(dataFiltered1_CDR3, decreasing = T)[1:intTopNum])
dataFiltered2_CDR3<-table(dataFiltered2[, "cdr3_nt"])
dataFiltered_CDR3[[2]]<-names(sort(dataFiltered2_CDR3, decreasing = T)[1:intTopNum])
dataFiltered3_CDR3<-table(dataFiltered3[, "cdr3_nt"])
dataFiltered_CDR3[[3]]<-names(sort(dataFiltered3_CDR3, decreasing = T)[1:intTopNum])
dataFiltered4_CDR3<-table(dataFiltered4[, "cdr3_nt"])
dataFiltered_CDR3[[4]]<-names(sort(dataFiltered4_CDR3, decreasing = T)[1:intTopNum])
dataFiltered5_CDR3<-table(dataFiltered5[, "cdr3_nt"])
dataFiltered_CDR3[[5]]<-names(sort(dataFiltered5_CDR3, decreasing = T)[1:intTopNum])
dataFiltered6_CDR3<-table(dataFiltered6[, "cdr3_nt"])
dataFiltered_CDR3[[6]]<-names(sort(dataFiltered6_CDR3, decreasing = T)[1:intTopNum])

cat("S1\tS2\toverlap\n")
for(ii in 1:5){
  for(jj in (ii+1):6){
    cat(ii);cat("\t");cat(jj);cat("\t")
    cat(length(intersect(dataFiltered_CDR3[[ii]], dataFiltered_CDR3[[jj]])));cat("\n")
  }
}

#S1	S2	overlap
#1	2	1
#1	3	1
#1	4	3
#1	5	1
#1	6	0
#2	3	7
#2	4	2
#2	5	1
#2	6	3
#3	4	3
#3	5	6
#3	6	2
#4	5	2
#4	6	0
#5	6	1


