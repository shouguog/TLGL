#setwd("H:/Bradley/Integration/DetailedAnalysis_round1/DiffusionMap")
rm(list=ls())
seed1<-111
dataDM<-read.csv(paste("../sva_diffusion_map_", seed1, "_1200cells.csv", sep=""), header=T, row.names = 1)
load("../../../datasets/allData.RData")
dataDM$cluster<-dataOutAll_meta_expression_tSNE_cluster[rownames(dataDM), "cluster_phenograph_combat"]
dataDM$metacluster<-dataOutAll_meta_expression_tSNE_cluster[rownames(dataDM), "cluster_phenographmeta_combat"]
coldefine=c("red", "orange", "yellow", "green","violet", "black", "blue", "pink", "grey", "cyan4")
dataDM_CD48<-dataDM####Shouguo gao[dataDM$metacluster %in% c(4,5,6),]
library(data.table)
###read AUC
AUCscore<-t(as.data.frame(fread("../../../../cellrangerGRCh38LGL7flows/Seurat/AUCellotherGSE/testAUC_Sample_all_otherGSE.csv")))
colnames(AUCscore)<-AUCscore[1,]; AUCscore<-AUCscore[-1,]
AUCscore_DM_CD48<-AUCscore[rownames(dataDM_CD48),]
dataDM_CD48<-cbind(dataDM_CD48, AUCscore_DM_CD48)
dataDM_CD48_new<-dataDM_CD48
for(ii in 13:26){
  dataDM_CD48_new[, ii]<-as.numeric(as.character(dataDM_CD48[, ii]))
}
rm(Achain_clonotypes,AUCscore,Bchain_clonotypes,coldefine,dataDM,AUCscore_DM_CD48,dataOutAll_meta_expression_tSNE_cluster,
  ExpressionCluster_combat_secondround,ExpressionCluster_limma_secondround,filtered_contig_annotations_all,LGLExpList180418,
  tenXgenomics_clonotypes,tenXgenomics_clonotypes_ABpair, dataDM_CD48)

corValue<-cor(dataDM_CD48_new[, c(1:10, 13:26)])
library(reshape2)
output<-setNames(melt(corValue), c('rowname', 'colname', 'correlation'))
output<-output[order(abs(output$correlation), decreasing = T),]
output<-output[output$correlation!=1,]  #remove self correlation
output[output$rowname=="V1",]
output<-output[output$rowname %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),]  #remove self correlation
output<-output[!(output$colname %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")),]  #remove self correlation
outputV1V2<-output[output$rowname %in% c("V1","V2"),]  #remove self correlation

dataDM_CD48_new<-dataDM_CD48_new[order(dataDM_CD48_new$V1),]
nbin<-300
DM1<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
for(ii in 1:length(DM1)){
  DM1[ii]<-mean(dataDM_CD48_new$V1[1:nbin+(ii-1)*nbin])
}

#for(cellset in colnames(dataDM_CD48_new)[13:26]){
#  png(paste("DM", seed1, "check_1200_CD48_DM1_", cellset, "_CD8.png", sep=""), width=2000, height=2000, res=120)
#  featureTrend<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
#  for(ii in 1:length(DM1)){
#    featureTrend[ii]<-mean(dataDM_CD48_new[1:nbin+(ii-1)*nbin, cellset])
#  }
#  plot(DM1, featureTrend, pch=19, xlab="DM1", ylab=cellset)
#  dev.off()
#}

colors<-c("hotpink","violet","turquoise3","black")
png(paste("DM", seed1, "check_1200_CD48_DM1_allset_CD8_newcolor.png", sep=""), width=1000, height=1000, res=200)
featureTrend<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
for(ii in 1:length(DM1)){
  featureTrend[ii]<-mean(dataDM_CD48_new[1:nbin+(ii-1)*nbin, "GSE93777_CD8_naivePvalue"])
}
featureTrend<-(featureTrend-min(featureTrend))/(max(featureTrend)-min(featureTrend))*0.9+0.05
colornum<-1
plot(DM1, featureTrend, pch=19, xlab="DM1", ylab="GSE93777_CD8_naivePvalue", col=colors[colornum], ylim=c(0,1), type="line", lwd=3)
#for(cellset in c("GSE93777_CD8_effectPvalue","GSE93777_CD8_centralPvalue","GSE93777_CD8_CD45ROPvalue")){
for(cellset in c("GSE93777_CD8_effectPvalue","GSE93777_CD8_centralPvalue")){
  featureTrend<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
  for(ii in 1:length(DM1)){
    featureTrend[ii]<-mean(dataDM_CD48_new[1:nbin+(ii-1)*nbin, cellset])
  }
  featureTrend<-(featureTrend-min(featureTrend))/(max(featureTrend)-min(featureTrend))*0.9+0.05
  colornum<-colornum+1
  lines(DM1, featureTrend, pch=19, xlab="DM1", ylab=cellset, col=colors[colornum], lwd=3)
}
#legend(-0.3, 0.7, legend=c("CD8_naive", "CD8_effect","CD8_central","CD8_CD45RO"),col=colors, lty=1:2, cex=0.8)
legend(-0.3, 0.7, legend=c("CD8_naive", "CD8_effect","CD8_central"),col=colors[1:3], lty=1:2, cex=0.8)
dev.off()
########plot DM2
dataDM_CD48_new<-dataDM_CD48_new[order(dataDM_CD48_new$V2),]
DM2<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
for(ii in 1:length(DM2)){
  DM2[ii]<-mean(dataDM_CD48_new$V2[1:nbin+(ii-1)*nbin])
}
colors<-c("hotpink","violet","turquoise3","black")
png(paste("DM", seed1, "check_1200_CD48_DM2_allset_CD8_newcolor.png", sep=""), width=1000, height=1000, res=200)
featureTrend<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
for(ii in 1:length(DM2)){
  featureTrend[ii]<-mean(dataDM_CD48_new[1:nbin+(ii-1)*nbin, "GSE93777_CD8_naivePvalue"])
}
featureTrend<-(featureTrend-min(featureTrend))/(max(featureTrend)-min(featureTrend))*0.9+0.05
colornum<-1
plot(DM2, featureTrend, pch=19, xlab="DM2", ylab="GSE93777_CD8_naivePvalue", col=colors[colornum], ylim=c(0,1), type="line", lwd=3)
#for(cellset in c("GSE93777_CD8_effectPvalue","GSE93777_CD8_centralPvalue","GSE93777_CD8_CD45ROPvalue")){
for(cellset in c("GSE93777_CD8_effectPvalue","GSE93777_CD8_centralPvalue")){
  featureTrend<-rep(0, floor(dim(dataDM_CD48_new)[1]/nbin))
  for(ii in 1:length(DM1)){
    featureTrend[ii]<-mean(dataDM_CD48_new[1:nbin+(ii-1)*nbin, cellset])
  }
  featureTrend<-(featureTrend-min(featureTrend))/(max(featureTrend)-min(featureTrend))*0.9+0.05
  colornum<-colornum+1
  lines(DM2, featureTrend, pch=19, xlab="DM2", ylab=cellset, col=colors[colornum], lwd=3)
}
#legend(-0.06, 0.7, legend=c("CD8_naive", "CD8_effect","CD8_central","CD8_CD45RO"),col=colors, lty=1:2, cex=0.8)
#legend(-0.1, 0.7, legend=c("CD8_naive", "CD8_effect","CD8_central"),col=colors[1:3], lty=1:2, cex=0.8)
dev.off()

