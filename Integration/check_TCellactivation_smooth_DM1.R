library(Seurat)
rm(list=ls())
###load all data
load("../../../../cellrangerGRCh38LGL7flows/Seurat/S.combined.RData")
Tactivation<-c("CD69","CCR7","CD27","BTLA","CD40LG","IL2RA","CD3E","CD47","EOMES","GNLY","GZMA","GZMB","PRF1","IFNG"
,"CD8A","CD8B","CD95L","LAMP1","LAG3","CTLA4","HLA-DRA","TNFRSF4","ICOS","TNFRSF9","TNFRSF18")
Tactivation<-intersect(Tactivation, rownames(S.combined@raw.data))

seed1<-111
##read diffusion result
data<-read.csv(paste("../sva_diffusion_map_", seed1, "_1200cells.csv", sep=""), header=T, row.names = 1)
dataExp<-t(as.matrix(S.combined@raw.data[Tactivation,rownames(data)]))
mean(cor(dataExp)[, "CD8A"])
###Let us plot the trend
dataplot<-as.numeric(rowSums(log(dataExp[order(data[,1]),]+1)))
dataplot10<-c()
for(ii in 1:3840){
  dataplot10<-c(dataplot10, sum(dataplot[1:10+ii*10]))
}

png(paste("meanTactivationMean_dm1_1200.png", sep=""), width=500, height=500, res=200)
###let us try to smooth
x<-1:3840
dataPlot<-data.frame(x=x, y=dataplot10)
loess_mod <- loess(y ~ x, dataPlot)
pred <- predict(loess_mod, dataPlot, se=TRUE)
dataPlot$lwl <- pred$fit-1.96*pred$se.fit
dataPlot$upl <- pred$fit+1.96*pred$se.fit
dataPlot<-dataPlot[dataPlot$x<3500,]
library(ggplot2)
print(ggplot(dataPlot, aes(x = x, y = y)) +
        #geom_point() +
        geom_smooth(method = 'loess') +
        geom_line(aes(y = lwl), color = "red") +
        geom_line(aes(y = upl), color = "red")+ theme_bw() +
        xlab("DM1") +
        ylab(("T activation")))
dev.off()

###Let us plot the trend
dataplot<-as.numeric(rowSums(log(dataExp[order(data[,2]),]+1)))
dataplot10<-c()
for(ii in 1:3840){
  dataplot10<-c(dataplot10, sum(dataplot[1:10+ii*10]))
}

png(paste("meanTactivationMean_dm2_1200.png", sep=""), width=500, height=500, res=200)
###let us try to smooth
x<-1:3840
dataPlot<-data.frame(x=x, y=dataplot10)
loess_mod <- loess(y ~ x, dataPlot)
pred <- predict(loess_mod, dataPlot, se=TRUE)
dataPlot$lwl <- pred$fit-1.96*pred$se.fit
dataPlot$upl <- pred$fit+1.96*pred$se.fit
dataPlot<-dataPlot[dataPlot$x<3500,]
library(ggplot2)
print(ggplot(dataPlot, aes(x = x, y = y)) +
        #geom_point() +
        geom_smooth(method = 'loess') +
        geom_line(aes(y = lwl), color = "red") +
        geom_line(aes(y = upl), color = "red")+ theme_bw() +
        xlab("DM2") +
        ylab(("T activation")))
dev.off()


