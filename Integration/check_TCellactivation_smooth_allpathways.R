library(Seurat)
rm(list=ls())
###load all data
load("../../../../cellrangerGRCh38LGL7flows/Seurat/S.combined.RData")
fileName="../danaSignature.gmt"
con=file(fileName,open="r")
line=readLines(con) 
long=length(line)
zz<-file("pathway_dim_cor.txt", "w")
cat("pathway\tdm\tcorR\n", file=zz);
for (i in 1:long){
  linn=line[i]
  genes<-strsplit(linn, "\t")[[1]]
  pathway<-genes[1]
  genes<-genes[-1]
  genes<-intersect(genes, rownames(S.combined@raw.data))
  seed1<-111
  ##read diffusion result
  data<-read.csv(paste("../sva_diffusion_map_", seed1, "_1200cells.csv", sep=""), header=T, row.names = 1)
  dataExp<-t(as.matrix(S.combined@raw.data[genes,rownames(data)]))
  ###Let us plot the trend
  for(dm in 1:10){
    dataplot<-as.numeric(rowSums(log(dataExp[order(data[,dm]),]+1)))
    dataplot10<-c()
    for(ii in 1:3840){
      dataplot10<-c(dataplot10, sum(dataplot[1:10+ii*10]))
    }
    png(paste(pathway, "Mean_dm", dm, "_1200.png", sep=""), width=500, height=500, res=200)
    ###let us try to smooth
    x<-1:3840
    dataPlot<-data.frame(x=x, y=dataplot10)
    loess_mod <- loess(y ~ x, dataPlot)
    pred <- predict(loess_mod, dataPlot, se=TRUE)
    dataPlot$lwl <- pred$fit-1.96*pred$se.fit
    dataPlot$upl <- pred$fit+1.96*pred$se.fit
    dataPlot<-dataPlot[dataPlot$x<3500,]
    print(ggplot(dataPlot, aes(x = x, y = y)) +
        #geom_point() +
        geom_smooth(method = 'loess') +
        geom_line(aes(y = lwl), color = "red") +
        geom_line(aes(y = upl), color = "red")+ theme_bw() +
        xlab(paste0("DM", dm)) +
        ylab(("T activation")))
    dev.off()
    cat(pathway, file=zz);cat("\t", file=zz); cat(dm, file=zz);cat("\t", file=zz); cat(cor(dataPlot$x, dataPlot$y), file=zz);cat("\n", file=zz);
  }
}
close(zz)

