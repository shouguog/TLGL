###R 3.4.1 
library(fgsea)
rm(list=ls())
pathways <- gmtPathways("apoptosis.gmt")
str(head(pathways))
dataFilter<-read.csv("../../resultPreOverTreatment.csv", header=T, row.names = 1)
####Order and put into a list
dataFilter<-dataFilter[order(dataFilter[,"statistic"]),]
ranksOur <- setNames(dataFilter[, "statistic"], rownames(dataFilter))
###And runnig fgsea:
fgseaRes <- fgsea(pathways, ranksOur, minSize=15, maxSize=1000, nperm=1000)
fgseaRes<-fgseaRes[order(fgseaRes[,"pval"]),]
write.csv(fgseaRes[,1:5], file=paste("fGSEA_resultPreOverTreatment.csv", sep=""))

###plot
library(ggplot2)
p1<-plotEnrichment(pathways[["GO:0043066"]],ranksOur) + labs(title="negative regulation of apoptotic process")
p2<-plotEnrichment(pathways[["GO:0042981"]],ranksOur) + labs(title="regulation of apoptotic process")
p3<-plotEnrichment(pathways[["GO:2001234"]],ranksOur) + labs(title="negative regulation of apoptotic signaling pathway")
p4<-plotEnrichment(pathways[["GO:0006915"]],ranksOur) + labs(title="apoptotic process]")
png("typicalpathways_resultPreOverTreatment.png", width=1000, height = 1000, res=200)
library(gridExtra)
print(grid.arrange(p1,p2,p3,p4,nrow=2))
dev.off()

rm(list=ls())
pathways <- gmtPathways("apoptosis.gmt")
str(head(pathways))
dataFilter<-read.csv("../../resultPreOverTreatment_CD4.csv", header=T, row.names = 1)
####Order and put into a list
dataFilter<-dataFilter[order(dataFilter[,"statistic"]),]
ranksOur <- setNames(dataFilter[, "statistic"], rownames(dataFilter))
###And runnig fgsea:
fgseaRes <- fgsea(pathways, ranksOur, minSize=15, maxSize=1000, nperm=1000)
fgseaRes<-fgseaRes[order(fgseaRes[,"pval"]),]
write.csv(fgseaRes[,1:5], file=paste("fGSEA_resultPreOverTreatment_CD4.csv", sep=""))

rm(list=ls())
pathways <- gmtPathways("apoptosis.gmt")
str(head(pathways))
dataFilter<-read.csv("../../resultPreOverTreatment_CD8.csv", header=T, row.names = 1)
####Order and put into a list
dataFilter<-dataFilter[order(dataFilter[,"statistic"]),]
ranksOur <- setNames(dataFilter[, "statistic"], rownames(dataFilter))
###And runnig fgsea:
fgseaRes <- fgsea(pathways, ranksOur, minSize=15, maxSize=1000, nperm=1000)
fgseaRes<-fgseaRes[order(fgseaRes[,"pval"]),]
write.csv(fgseaRes[,1:5], file=paste("fGSEA_resultPreOverTreatment_CD8.csv", sep=""))
