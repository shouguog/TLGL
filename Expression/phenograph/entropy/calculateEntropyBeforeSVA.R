rm(list=ls())
dataCluster<-read.csv("../../data_info_res.4.csv")
tableEntropy<-as.data.frame(table(dataCluster[, c("orig.ident", "res.0.5")]))
tableSize<-as.data.frame(table(dataCluster[, c("orig.ident")]))
rownames(tableSize)<-tableSize$Var1
tableEntropy$ratio<-tableEntropy$Freq/tableSize[tableEntropy$orig.ident,"Freq"]

zz<-file("entropy_BeforeSVA_cluster.txt", "w")
cat("cluster\tentropy\n", file=zz)
entro<-c()
for(cluster in unique(tableEntropy$res.0.5)){
  tableEntropyCluster<-tableEntropy[tableEntropy$res.0.5==cluster,]
  ###Normalize to sum 1
  tableEntropyCluster$ratio<-tableEntropyCluster$ratio/sum(tableEntropyCluster$ratio)
  H<-0
  for(ii in 1:32){
    if(tableEntropyCluster[ii,"ratio"]>0){
      H<-H-tableEntropyCluster[ii,"ratio"]*log2(tableEntropyCluster[ii,"ratio"])
    }
  }
  cat(cluster, file=zz); cat("\t", file=zz);cat(H, file=zz);cat("\n", file=zz)
  entro<-c(entro,H)
}
close(zz)

png("entropy_distribution_beforeSVA.png", width=1000, height=1000, res=150)
barplot(sort(entro))
dev.off()
