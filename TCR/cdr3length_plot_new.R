setwd("H:/Bradley/Integration/DetailedAnalysis_round1/CDR3/length")
#Let us filter combine the data of 32
rm(list=ls())
load("../../../datasets/allData.RData")
CD8cells<-dataOutAll_meta_expression_tSNE_cluster[dataOutAll_meta_expression_tSNE_cluster$cluster_phenographmeta_combat%in%c(0,4,5,6,8),c("orig.ident", "CD8A")]
CD4cells<-dataOutAll_meta_expression_tSNE_cluster[dataOutAll_meta_expression_tSNE_cluster$cluster_phenographmeta_combat%in%c(1,2,3,7,9),c("orig.ident", "CD8A")]
cellCountCD4CD8<-data.frame(CD8cells=table(CD8cells$orig.ident), CD4cells=table(CD4cells$orig.ident), cellsAll=table(dataOutAll_meta_expression_tSNE_cluster$orig.ident))
rownames(cellCountCD4CD8)<-names(table(CD4cells$orig.ident))
cellOverCode<-intersect(rownames(CD8cells), filtered_contig_annotations_all[, "expressionBarcode"])
filtered_contig_annotations_all$cdr3length=nchar(filtered_contig_annotations_all$cdr3)

###Let us check CD8 cells
filtered_contig_annotations_all_CD8<-filtered_contig_annotations_all[filtered_contig_annotations_all[, "expressionBarcode"] %in% intersect(rownames(CD8cells), filtered_contig_annotations_all[, "expressionBarcode"]),]
###Let us check CD4 cells
filtered_contig_annotations_all_CD4<-filtered_contig_annotations_all[filtered_contig_annotations_all[, "expressionBarcode"] %in% intersect(rownames(CD4cells), filtered_contig_annotations_all[, "expressionBarcode"]),]


###Let us check CD8 cells
png("cdr3length_CD8_CD4_length_AA.png", width=4000, height=2000, res=150)
par(mfrow=c(4, 8))
for(s in paste("S", 1:32, sep="")){
  cat(s); cat("\n")
  dataPlot<-data.frame(x=8:54*0.5, y=0)
  for(ii in 1:47){
    dataPlot[ii,"y"]<-dim(filtered_contig_annotations_all_CD8[filtered_contig_annotations_all_CD8$sample==s & filtered_contig_annotations_all_CD8$cdr3length==dataPlot[ii,"x"],])[1]
  }
  plot(dataPlot$x, dataPlot$y, xlab="CDR3length", type="l", col="red")
  dataPlot<-data.frame(x=8:54*0.5, y=0)
  for(ii in 1:47){
    dataPlot[ii,"y"]<-dim(filtered_contig_annotations_all_CD4[filtered_contig_annotations_all_CD4$sample==s & filtered_contig_annotations_all_CD4$cdr3length==dataPlot[ii,"x"],])[1]
  }
  lines(dataPlot$x, dataPlot$y, xlab="CDR3length", type="l", col="blue")
}
dev.off()

###Let us check CD8 cells
for(s in paste("S", 1:32, sep="")){
  cat(s); cat("\n")
  png(paste0("samples/", s, "_cdr3length_CD8_CD4_length_AA.png"), width=1000, height=1000, res=250)
  dataPlot<-data.frame(x=8:54*0.5, y=0)
  for(ii in 1:47){
    dataPlot[ii,"y"]<-dim(filtered_contig_annotations_all_CD8[filtered_contig_annotations_all_CD8$sample==s & filtered_contig_annotations_all_CD8$cdr3length==dataPlot[ii,"x"],])[1]
  }
  plot(dataPlot$x, dataPlot$y, xlab="CDR3length", type="l", col="red")
  dataPlot<-data.frame(x=8:54*0.5, y=0)
  for(ii in 1:47){
    dataPlot[ii,"y"]<-dim(filtered_contig_annotations_all_CD4[filtered_contig_annotations_all_CD4$sample==s & filtered_contig_annotations_all_CD4$cdr3length==dataPlot[ii,"x"],])[1]
  }
  lines(dataPlot$x, dataPlot$y, xlab="CDR3length", type="l", col="blue")
  dev.off()
}
