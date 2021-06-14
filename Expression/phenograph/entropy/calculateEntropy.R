rm(list=ls())
tabelMeta<-read.csv("../tabelNoMetaAllcluster.csv", header=TRUE)
zz<-file("entropy_NoMeta_cluster.txt", "w")
cat("cluster\tentropy\n", file=zz)
entro<-c()
for(cluster in unique(tabelMeta$cluster_phenograph_combat)){
  tabelMetaCluster<-tabelMeta[tabelMeta$cluster_phenograph_combat==cluster,]
  ###Normalize to sum 1
  tabelMetaCluster$ratioAll<-tabelMetaCluster$ratioAll/sum(tabelMetaCluster$ratioAll)
  H<-0
  for(ii in 1:32){
    if(tabelMetaCluster[ii,"ratioAll"]>0){
      H<-H-tabelMetaCluster[ii,"ratioAll"]*log2(tabelMetaCluster[ii,"ratioAll"])
    }
  }
  cat(cluster, file=zz); cat("\t", file=zz);cat(H, file=zz);cat("\n", file=zz)
  entro<-c(entro,H)
}
close(zz)

png("entropy_distribution.png", width=1000, height=1000, res=150)
barplot(sort(entro))
dev.off()

permData<-c()
for(jj in 1:5000){
  x1 <- runif(32, 0, 1)
  x1 <- x1/sum(x1)
  H<-0
  for(ii in 1:32){
    if(x1[ii]>0){
      H<-H-x1[ii]*log2(x1[ii])
    }
  }
  permData<-c(permData, H)
}
png("entropy_permutation_distribution.png", width=2000, height=1000, res=150)
par(mfrow=c(1,2))
barplot(sort(entro), xlab="cluster", ylab="Entropy")
hist(permData, xlab="Entropy", main=paste("5%", round(quantile(permData, seq(0, 1, 0.05))[2],2), sep=":"))
dev.off()

##quantile(permData, seq(0, 1, 0.05))
