setwd("H:/Bradley/Integration/powerlaw/pairedchain")
###The following code is used to check the power law identification
library("poweRlaw")
rm(list=ls())
####Let us read the patient information
patientInfo<-read.csv("../../../180418LGLExpList.csv", header=T)
patientInfo<-paste(patientInfo[,1], patientInfo[,2], patientInfo[,3],patientInfo[,4], sep=":")
png("Logpowerlaw_pairedchain_part1.png", width=2000, height=2000, res=100)
par(mfrow=c(4,4))
for(sample in 1:16){
  data<-read.csv(paste("../../cloneDefination/10xgenomics/10xgenomics_", sample, "_ABpair_summary.csv", sep=""), header=T)
  data<-data[data[, "x"]>0,]
  ##create object
  m_pl = dislnorm$new(data[, "x"])
  ###estimate Xmin
  est = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  ####Plot fitting plots
  plot(m_pl, main=patientInfo[sample])
  lines(m_pl, col=2)
  #bs_p = bootstrap_p(m_pl)
  #bs_p$p
}
dev.off()

png("Logpowerlaw_pairedchain_part2.png", width=2000, height=2000, res=100)
par(mfrow=c(4,4))
for(sample in 17:32){
  data<-read.csv(paste("../../cloneDefination/10xgenomics/10xgenomics_", sample, "_ABpair_summary.csv", sep=""), header=T)
  data<-data[data[, "x"]>0,]
  ##create object
  m_pl = dislnorm$new(data[, "x"])
  ###estimate Xmin
  est = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  ####Plot fitting plots
  plot(m_pl, main=patientInfo[sample])
  lines(m_pl, col=2)
  #bs_p = bootstrap_p(m_pl)
  #bs_p$p
}
dev.off()

