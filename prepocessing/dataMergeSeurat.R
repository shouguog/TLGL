#setwd("/Users/gaos2/Desktop/mywork/Bradley/cellrangerGRCh38LGL")
#setwd("H:/Bradley/cellrangerGRCh38LGL")
setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL")
library(Seurat)
rm(list=ls())
#S1.data <- Read10X(data.dir = "SfivePrimeR2GRCh38_1/outs/filtered_gene_bc_matrices/GRCh38/")
#S1 <- CreateSeuratObject(raw.data = S1.data, project = "S1")
#S2.data <- Read10X(data.dir = "SfivePrimeR2GRCh38_2/outs/filtered_gene_bc_matrices/GRCh38/")
#S2 <- CreateSeuratObject(raw.data = S2.data, project = "S2")
###merge S1  S2
#S.combined <- MergeSeurat(object1 = S1, object2 = S2, add.cell.id1 = "S1", add.cell.id2 = "S2", project = "LGL")

S.combined.data <- Read10X(data.dir = "SfivePrimeR2GRCh38_1/outs/filtered_gene_bc_matrices/GRCh38/")
S.combined <- CreateSeuratObject(raw.data = S.combined.data, project = "S1")

for(sample in 2:11){
  cat("sample"); cat(sample); cat("\n");
  S.data <- Read10X(data.dir = paste("SfivePrimeR2GRCh38_", sample, "/outs/filtered_gene_bc_matrices/GRCh38/", sep=""))
  S.combined <- AddSamples(object = S.combined, new.data = S.data, add.cell.id = paste("S", sample,sep=""))
}
for(sample in 12:32){
  cat("sample"); cat(sample); cat("\n");
  S.data <- Read10X(data.dir = paste("SfivePrimeR2GRCh38_", sample, "/outs/filtered_gene_bc_matrices/GRCh38/", sep=""))
  S.combined <- AddSamples(object = S.combined, new.data = S.data, add.cell.id = paste("S", sample,sep=""))
}

save(S.combined, file="S.combined.RData")
# notice the cell names now have an added identifier
head(x = S.combined@cell.names)
table(S.combined@meta.data$orig.ident)

