#setwd("/Users/gaos2/Desktop/mywork/Bradley/cellrangerGRCh38LGL")
#setwd("H:/Bradley/cellrangerGRCh38LGL")
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat")
library(Seurat)
rm(list=ls())
load("S.combined.RData")
# notice the cell names now have an added identifier
head(x = S.combined@cell.names)
table(S.combined@meta.data$orig.ident)
#####Let us normalize the data
S.combined <- NormalizeData(S.combined)
png("varGenePlot.png")
S.combined <- FindVariableGenes(S.combined, do.plot = T, y.cutoff = 0.4)
dev.off()

####union
vargene<-S.combined@var.genes
S.combined <- ScaleData(object = S.combined, genes.use=vargene, vars.to.regress = c("nUMI"))
save(S.combined,file = "S.combined_step1.Robj")

S.combined <- RunPCA(object = S.combined, pc.genes = S.combined@var.genes, pcs.compute = 25, do.print = TRUE, pcs.print = 1:5,  genes.print = 5)

#select 25 PCs for downstream analysis
S.combined <- RunTSNE(S.combined, dims.use = 1:25, do.fast = T)
save(S.combined,file = "S.combined_step2.Robj")

#save.SNN means that you can easily re-run with different resolution values. Here we run with a few different res v
S.combined <- FindClusters(S.combined ,dims.use = 1:25, resolution = seq(0.5,4,0.5), save.SNN = F, do.sparse = T)

save(S.combined,file = "S.combined_step3.Robj")

######Let us output the result with different resolution
S.combined <- SetAllIdent(S.combined, id = "res.0.5")
write.csv(S.combined@data.info, file="data_info_res.0.5.csv")

