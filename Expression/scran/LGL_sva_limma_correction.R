####We did this because we found that Seurat can not remove the effect of PtID if we combine P1_p4 and C1-C3
####Seurat can not make cells distribute into different clusters. may be right, may be not
####Seurat for >2 groups is still underdevelopment. MNNCorrect in scran can do this for several datasets 
####Seurat 2.0 result can be found in the according folder
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat/scran")
library(Matrix)
library(sva)
library(Rtsne)
library(scran)
rm(list=ls())
load("sce.Robj")
###sva batch
expressionAll<-exprs(sceObj)   #dim(expressionAll)
expressionAll<-expressionAll[, rownames(metaData)]##make sure they have same orders
samples<-as.character(metaData[, "orig.ident"])
####correct batch
allSamples.combat <- ComBat(expressionAll,samples,mod=NULL,prior.plots = FALSE)
#all.dists.combat <- as.matrix(dist(t(allSamples.combat)))
#set.seed(0)
#tsne.combat <- Rtsne(all.dists.combat, is_distance=TRUE)
### limma batch correction
library(limma)
allSamples.limma <- removeBatchEffect(expressionAll, factor(samples))
#all.dists.limma <- as.matrix(dist(t(allSamples.limma)))
#set.seed(0)
#tsne.limma <- Rtsne(all.dists.limma, is_distance=TRUE)
save(list=ls(), file="sce_sva_limma_correct.Robj")

