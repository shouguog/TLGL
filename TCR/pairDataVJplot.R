setwd("H:/Bradley/Integration/DetailedAnalysis_round1/clonaltypesummary/paired1A1Bdata")
#Let us filter combine the data of 32
rm(list=ls())
library(latticeExtra)
patientInfo<-read.csv("../../../../180418LGLExpList.csv", row.names = 1)
patientInfoAll<-paste(patientInfo[,1],patientInfo[,2],patientInfo[,3], sep="_")
names(patientInfoAll)<-rownames(patientInfo)
for(sample in paste("S", 1:32, sep="")){
  dataFiltered1A1B<-read.csv(paste("filtered_contig_annotations_", sample, ".csv", sep=""), header=T)
  dataPlot<-as.matrix(table(dataFiltered1A1B[, c("v_gene", "j_gene")]))
  dataPlot3<-as.data.frame(matrix(1, nrow=dim(dataPlot)[1]*dim(dataPlot)[2], ncol=3))
  colnames(dataPlot3)<-c("x", "y", "z")
  for(ii in 1:dim(dataPlot)[1]){
    for(jj in 1:dim(dataPlot)[2]){
      dataPlot3[(ii-1)*dim(dataPlot)[2]+jj, 1]<-rownames(dataPlot)[ii]
      dataPlot3[(ii-1)*dim(dataPlot)[2]+jj, 2]<-as.integer(dataPlot[ii,jj])
      dataPlot3[(ii-1)*dim(dataPlot)[2]+jj, 3]<-colnames(dataPlot)[jj]
    }
  }
  dataPlot3$x<-as.factor(dataPlot3$x)
  dataPlot3$y<-as.integer(dataPlot3$y)
  dataPlot3$z<-as.factor(dataPlot3$z)
  #png(paste("result/", sample,  patientInfoAll[sample], ".png", sep=""), width=3000, height=3000, res=50)
  #print(cloud(y~x+z, dataPlot3, panel.3d.cloud=panel.3dbars, col.facet='red', 
  #    xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
  #    par.settings = list(axis.line = list(col = "transparent")), aspect=c(0.3,0.3)))
  #dev.off()
  dataPlot3<-dataPlot3[order(dataPlot3$y, decreasing = TRUE),]
  write.csv(dataPlot3, file=paste("result/", sample,  patientInfoAll[sample], "_barPlot.csv", sep=""))
}

#d <- read.table(text=' x   y     z
#t1   5   high
#t1   2   low
#t1   4   med
#t2   8   high
#t2   1   low
#t2   3   med
#t3  50   high
#t3  12   med
#t3  35   low', header=TRUE)

#library(latticeExtra)
#png("extralattice.png", width=2000, height=600, res=200)
#cloud(y~x+z, d, panel.3d.cloud=panel.3dbars, col.facet='red', 
#      xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
#      par.settings = list(axis.line = list(col = "transparent")), aspect=c(0.3,0.3))
#dev.off()



