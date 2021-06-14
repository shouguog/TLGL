rm(list=ls())
library(topGO)
library(org.Hs.eg.db)
###let us read the gene mapping
idmapping<-read.table("../../geneexpression/Pre_vs_treatment/genes.tsv")
result<-read.csv("resultPreOverTreatment.csv", row.names=1)
geneNames<-rownames(result)
myInterestingGenes<-rownames(result[result[, "pvalue"]<0.01 & result[, "statistic"]>0,])  #statistic positive means higher in Pre

geneNames<-as.character(idmapping[idmapping[,2]%in%geneNames,1])
myInterestingGenes <- as.character(idmapping[idmapping[,2]%in%myInterestingGenes,1])

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
score(resultFisher)
goscore<-sort(score(resultFisher))
zz <- file("topGO_higherInPre.txt", "w") 
cat("GOid\tgoterm\tcountInteresting\tInterestingGeneList\tcountAll\tpvalue\n", file = zz, sep="")
for (i in 1:500){
	cat(as.vector(names(goscore)[i][1]), file=zz)
	cat("\t", file=zz)
	cat(Term(names(goscore)[i][1]), file=zz)
	cat("\t", file=zz)
	genelist<-intersect(as.vector(genesInTerm(GOdata, names(goscore)[i][1])[[names(goscore)[i][1]]]), myInterestingGenes)
	genelistAll<-intersect(as.vector(genesInTerm(GOdata, names(goscore)[i][1])[[names(goscore)[i][1]]]), geneNames)
	cat(length(genelist), file=zz)
	cat("\t", file=zz)
	###write gene list
	cat(genelist[1], file=zz)
	for(j in 2:length(genelist)){
		cat("|", file=zz)
		cat(genelist[j], file=zz)
	}
	cat("\t", file=zz)
	cat(length(genelistAll), file=zz)
	cat("\t", file=zz)
	cat(goscore[i], file=zz)
	cat("\n", file=zz)
}
close(zz)


rm(list=ls())
library(topGO)
library(org.Hs.eg.db)
###let us read the gene mapping
idmapping<-read.table("../../geneexpression/Pre_vs_treatment/genes.tsv")
result<-read.csv("resultPreOverTreatment.csv", row.names=1)
geneNames<-rownames(result)
myInterestingGenes<-rownames(result[result[, "pvalue"]<0.001 & result[, "statistic"]<0,])

geneNames<-as.character(idmapping[idmapping[,2]%in%geneNames,1])
myInterestingGenes <- as.character(idmapping[idmapping[,2]%in%myInterestingGenes,1])

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
score(resultFisher)
goscore<-sort(score(resultFisher))
zz <- file("topGO_HigherinTreatment.txt", "w") 
cat("GOid\tgoterm\tcountInteresting\tInterestingGeneList\tcountAll\tpvalue\n", file = zz, sep="")
for (i in 1:500){
  cat(as.vector(names(goscore)[i][1]), file=zz)
  cat("\t", file=zz)
  cat(Term(names(goscore)[i][1]), file=zz)
  cat("\t", file=zz)
  genelist<-intersect(as.vector(genesInTerm(GOdata, names(goscore)[i][1])[[names(goscore)[i][1]]]), myInterestingGenes)
  genelistAll<-intersect(as.vector(genesInTerm(GOdata, names(goscore)[i][1])[[names(goscore)[i][1]]]), geneNames)
  cat(length(genelist), file=zz)
  cat("\t", file=zz)
  ###write gene list
  cat(genelist[1], file=zz)
  for(j in 2:length(genelist)){
    cat("|", file=zz)
    cat(genelist[j], file=zz)
  }
  cat("\t", file=zz)
  cat(length(genelistAll), file=zz)
  cat("\t", file=zz)
  cat(goscore[i], file=zz)
  cat("\n", file=zz)
}
close(zz)




