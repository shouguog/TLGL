####We did this because we found that Seurat can not remove the effect of PtID if we combine P1_p4 and C1-C3
####Seurat can not make cells distribute into different clusters. may be right, may be not
####Seurat for >2 groups is still underdevelopment. MNNCorrect in scran can do this for several datasets 
####Seurat 2.0 result can be found in the according folder
library(Matrix)
library(scran)
library(Seurat)
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat/scran")
rm(list=ls())
load("../S.combined_step5.Robj")
countData<-S.combined@raw.data[S.combined@var.genes,]
metaData<-S.combined@meta.data
rm(S.combined)
sceObj<-calculateQCMetrics(newSCESet(countData=data.frame(as.matrix(countData))))
save(metaData, countData, sceObj, file="sce.Robj")
