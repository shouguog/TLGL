setwd("H:/Bradley/Integration/DetailedAnalysis_round1/DiffusionMap")

library(Matrix)
library(sva)
library(Rtsne)
library(scran)
library(limma)
library(diffusionMap)
rm(list=ls())
#load("../../datasets/scran/sce_sva_limma_correct.Robj")
#rm(list=c("allSamples.limma", "countData","expressionAll","metaData","samples","sceObj"))
#save(allSamples.combat, file="allSamples.combat.RData")

load("allSamples.combat.RData")
cat("load data\n")
dataExpression=as.matrix(t(allSamples.combat))
numExtract<-200
cells<-rownames(dataExpression)
###For sample 1
set.seed(666)
cells1<-cells[!grepl("_", cells)]
cellselected<-cells1[sample(1:length(cells1), length(cells1))[1:numExtract]]
for(ii in 2:32){
  set.seed(666)
  cells1<-cells[grepl(paste("S", ii, "_", sep=""), cells)]
  cellselected<-c(cellselected,cells1[sample(1:length(cells1), length(cells1))[1:numExtract]])
}
cat("cellPerm:");cat(length(cellselected)); cat("\n");
D = 1-cor(t(dataExpression[cellselected,]), method="pearson")
D = as.dist(D)
cat("calculate distance\n")
dmap = diffuse(D,neigen=10) # compute diffusion map
cat("Finished mapping\n")
resultMapping<-dmap$X
rownames(resultMapping)<-cellselected
write.csv(resultMapping, file="sva_diffusion_map_200.csv")


###The information for clonal TCR
TCRinfomration<-read.csv("../../datasets/filtered_contig_annotations_all.csv")
length(unique(TCRinfomration$expressionBarcode))
length(cells)
length(intersect(unique(TCRinfomration$expressionBarcode), cells))
