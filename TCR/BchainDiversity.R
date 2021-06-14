setwd("H:/Bradley/Integration/diversity/Bchain")
library(tcR)
library(vegan)
rm(list=ls())
data<-read.csv(paste("../../cloneDefination/Bchain/Bchain_", 1, "_summary.csv", sep=""))
indexResult<-c(inverse.simpson(data[,2]/sum(data[,2]), TRUE, 0),
               diversity(data[,2]/sum(data[,2])),
               gini(data[,2]/sum(data[,2]), TRUE, 0),  
               gini.simpson(data[,2]/sum(data[,2]), TRUE, 0),
               diversity(data[,2]),
               diversity(data[,2], index="simpson"),
               diversity(data[,2], index="invsimpson")
)
dataOut<-matrix(indexResult, nrow=1)
colnames(dataOut)<-c("inverse.simpson_tCR","diversity_tCR","gini_tCR","gini.simpson_tCR","shannon_vega","simpson_vegan","invsimpson_vegan")
for(ii in 2:32){
  dataItem<-read.csv(paste("../../cloneDefination/Bchain/Bchain_", ii, "_summary.csv", sep=""))
  indexResult<-c(inverse.simpson(dataItem[,2]/sum(dataItem[,2]), TRUE, 0),
                 diversity(dataItem[,2]/sum(dataItem[,2])),
                 gini(dataItem[,2]/sum(dataItem[,2]), TRUE, 0),  
                 gini.simpson(dataItem[,2]/sum(dataItem[,2]), TRUE, 0),
                 diversity(dataItem[,2]),
                 diversity(dataItem[,2], index="simpson"),
                 diversity(dataItem[,2], index="invsimpson")
  )
  dataOut<-rbind(dataOut,matrix(indexResult, nrow=1))
}

patientInfo<-read.csv("../../../180418LGLExpList.csv", row.names = 1)
PInfo<-paste(patientInfo[,"PatientID"], patientInfo[,"Treatment"], patientInfo[,"Age_Sex"], sep="#")
rownames(dataOut)<-PInfo
write.csv(dataOut, file="diversityBchain_all.csv")
plot(dataOut[, "shannon_vega"], dataOut[,"simpson_vegan"])
