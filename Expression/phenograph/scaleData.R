setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat/scran")
rm(list=ls())
library(Seurat)
load("../S.combined.RData")
###Let scale the gene expression
geneCount<-rowSums(S.combined@raw.data)
geneCount10000<-geneCount[geneCount>10000]
S.combined<-ScaleData(S.combined, genes.use = names(geneCount10000))
scaledDataAll<-S.combined@scale.data
save(scaledDataAll, file="scaleDataAll.RData")

