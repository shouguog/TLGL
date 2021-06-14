####We did this because we found that Seurat can not remove the effect of PtID if we combine P1_p4 and C1-C3
####Seurat can not make cells distribute into different clusters. may be right, may be not
####Seurat for >2 groups is still underdevelopment. MNNCorrect in scran can do this for several datasets 
####Seurat 2.0 result can be found in the according folder
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat/scran")
library(Matrix)
library(scran)
rm(list=ls())
?mnnCorrect
load("sce.Robj")
expressionAll<-exprs(sceObj)   #dim(expressionAll)
### MNN batch correction
expressionList<-list()
for(ii in 1:32){
  metaDataSample<-rownames(metaData[metaData[,"orig.ident"]==paste("S", ii, sep=""),])
  cells<-intersect(metaDataSample, colnames(expressionAll))
  expressionList[[ii]]<-expressionAll[,cells]
  cat(dim(expressionList[[ii]])[1]); cat("\n");
}

LGLXmnn <- mnnCorrect(expressionList[[1]],expressionList[[2]],expressionList[[3]],expressionList[[4]],expressionList[[5]],
expressionList[[6]],expressionList[[7]],expressionList[[8]],expressionList[[9]],expressionList[[10]],
expressionList[[11]],expressionList[[12]],expressionList[[13]],expressionList[[14]],expressionList[[15]],
expressionList[[16]],expressionList[[17]],expressionList[[18]],expressionList[[19]],expressionList[[20]],
expressionList[[21]],expressionList[[22]],expressionList[[23]],expressionList[[24]],expressionList[[25]],
expressionList[[26]],expressionList[[27]],expressionList[[28]],expressionList[[29]],expressionList[[30]],
expressionList[[31]],expressionList[[32]],k=40, sigma=0.5,svd.dim=0) # batch correction is throwing an error because of non-numeric values?
save(sceObj, metaData,LGLXmnn, file="sce_mnncorrect.Robj")


