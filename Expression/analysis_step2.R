#setwd("/Users/gaos2/Desktop/mywork/Bradley/cellrangerGRCh38LGL")
#setwd("H:/Bradley/cellrangerGRCh38LGL")
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat")
library(Seurat)
rm(list=ls())
load("S.combined_step2.Robj")
#save.SNN means that you can easily re-run with different resolution values. Here we run with a few different res v
S.combined <- FindClusters(S.combined ,dims.use = 1:25, resolution = 0.5, save.SNN = T)
save(S.combined,file = "S.combined_step3.Robj")
S.combined <- FindClusters(S.combined ,dims.use = 1:25, resolution = 4, save.SNN = T)
save(S.combined,file = "S.combined_step4.Robj")
######Let us output the result with different resolution
S.combined <- SetAllIdent(S.combined, id = "res.0.5")
write.csv(S.combined@data.info, file="data_info_res.0.5.csv")

