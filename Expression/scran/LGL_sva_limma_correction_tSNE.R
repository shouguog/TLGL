####We did this because we found that Seurat can not remove the effect of PtID if we combine P1_p4 and C1-C3
####Seurat can not make cells distribute into different clusters. may be right, may be not
####Seurat for >2 groups is still underdevelopment. MNNCorrect in scran can do this for several datasets 
####Seurat 2.0 result can be found in the according folder
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat/scran")
library(Matrix)
library(sva)
library(Rtsne)
library(scran)
library(limma)
rm(list=ls())
load("sce_sva_limma_correct.Robj")
####combat
set.seed(0)
tsne.combat <- Rtsne(as.matrix(t(allSamples.combat)), is_distance=FALSE,pca = TRUE)
set.seed(0)
tsne.limma <- Rtsne(as.matrix(t(allSamples.limma)), is_distance=FALSE,pca = TRUE)
save(tsne.combat, tsne.limma, file="sce_sva_limma_tSNE.Robj")


##iris_unique <- unique(iris) # Remove duplicates
##set.seed(42) # Sets seed for reproducibility
##tsne_out <- Rtsne(as.matrix(iris_unique[,1:4])) # Run TSNE
##plot(tsne_out$Y,col=iris_unique$Species) # Plot the result

