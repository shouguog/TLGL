setwd("/data/gaos2/10xgenomics/cellrangerGRCh38LGL7flows/Seurat/plotGeneExpression/DEanalysis")
rm(list=ls())
library(Seurat)
load("../../../../../Bradley/Integration/datasets/allData_withRegev_GSE118165.RData")
rm(Achain_clonotypes,Bchain_clonotypes,ExpressionCluster_combat_secondround,ExpressionCluster_limma_secondround,filtered_contig_annotations_all
   ,LGLExpList180418,tenXgenomics_clonotypes,tenXgenomics_clonotypes_ABpair)
rm(dataGSE118165, dataRegev)

####read the gene expression after sva normalization
load("../S.combined_CD48.Robj")
###We save CD4 and CD8 to save loading time
cellCD8<-rownames(dataOutAll_meta_expression_tSNE_cluster[dataOutAll_meta_expression_tSNE_cluster$cluster_phenographmeta_combat%in%c(4,5,6),])
S.combined.CD8<-SubsetData(S.combined.CD48, cells.use = cellCD8)
S.combined.CD8@raw.data<-S.combined.CD48@raw.data[,cellCD8];
cat("finished")
#Assign cluster number
S.combined.CD8@meta.data$cluster<-dataOutAll_meta_expression_tSNE_cluster[cellCD8,]$cluster_phenograph_combat
rm(S.combined.CD48)
unique(S.combined.CD8@meta.data$cluster)
###DE genes
S.combined.CD8 <- SetAllIdent(S.combined.CD8, id = "cluster")
set.seed(666)
cellCD8<-cellCD8[sample(1:length(cellCD8),150000)]
#S.combined.CD8.sub<-SubsetData(S.combined.CD8, cells.use = cellCD8[1:50000])
#for(ii in unique(S.combined.CD8@meta.data$cluster)){
#  cat(ii); cat("\n")
#  markers <- FindMarkers(object = S.combined.CD8.sub,  ident.1 = ii, logfc.threshold = 0.02, only.pos=TRUE)
#  #head(markers)
#  write.csv(markers, file=paste("CD8_markers_cluster_beforemeta", ii, "_perm1.csv", sep=""))
#}
S.combined.CD8.sub<-SubsetData(S.combined.CD8, cells.use = cellCD8[1:50000+50000])
for(ii in unique(S.combined.CD8@meta.data$cluster)){
  cat(ii); cat("\n")
  markers <- FindMarkers(object = S.combined.CD8.sub,  ident.1 = ii, logfc.threshold = 0.02, only.pos=TRUE)
  #head(markers)
  write.csv(markers, file=paste("CD8_markers_cluster_beforemeta", ii, "_perm2.csv", sep=""))
}

