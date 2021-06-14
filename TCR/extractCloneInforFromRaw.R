setwd("H:/Bradley/Integration/cloneDefination")
rm(list=ls())
for(sample in 1:32){
  dataFiltered<-read.csv(file=paste("../../cellrangerVDJ/SNOLANE1_", sample, "_VDJT/outs/filtered_contig_annotations.csv", sep=""), header=T)
  dataFiltered<-dataFiltered[dataFiltered[, "raw_consensus_id"]!="None",]
  dataAnnotation<-read.csv(file=paste("../../cellrangerVDJ/SNOLANE1_", sample, "_VDJT/outs/consensus_annotations.csv", sep=""), header=T, row.names=2)
  dataAnnotationCombined<-dataAnnotation[as.character(dataFiltered[,"raw_consensus_id"]),]
  data<-cbind(dataFiltered, dataAnnotationCombined)
  #####We need to rename the colnames
  colnames(data)<-c("barcode","is_cell","contig_id","high_confidence","Filtered_contig_annotations_length",
  "Filtered_contig_annotations_chain","Filtered_contig_annotations_v_gene","Filtered_contig_annotations_d_gene"
,"Filtered_contig_annotations_j_gene","Filtered_contig_annotations_c_gene","Filtered_contig_annotations_full_length"
,"Filtered_contig_annotations_productive","Filtered_contig_annotations_cdr3","Filtered_contig_annotations_cdr3_nt",
"Filtered_contig_annotations_reads","Filtered_contig_annotations_umis","raw_clonotype_id","raw_consensus_id","clonotype_id",
"Consensus_annotations_length","Consensus_annotations_chain","Consensus_annotations_v_gene","Consensus_annotations_d_gene",
"Consensus_annotations_j_gene","Consensus_annotations_c_gene","Consensus_annotations_full_length","Consensus_annotations_productive",
"Consensus_annotations_cdr3","Consensus_annotations_cdr3_nt","Consensus_annotations_reads","Consensus_annotations_umis")
  write.csv(data, file=paste("combined/combined_S", sample, ".csv", sep=""))
  #####Let us check the relationship between clonetype and consensus
  data_cloneID_cloneconsensus<-data[, c("clonotype_id", "raw_consensus_id","Consensus_annotations_v_gene","Consensus_annotations_d_gene"
                                        ,"Consensus_annotations_j_gene","Consensus_annotations_c_gene", "Filtered_contig_annotations_cdr3_nt")]
  data_cloneID_cloneconsensus_backup<-data_cloneID_cloneconsensus
  data_cloneID_cloneconsensus<-unique(data_cloneID_cloneconsensus)
  ####Let us combine the clonotype
  allClones<-unique(as.character(data_cloneID_cloneconsensus[,"clonotype_id"]))
  allClonesInfo<-c();allClonesCount<-c(); CloneChainCount<-c();
  for(clone in allClones){
    data_cloneID_cloneconsensus_clone=data_cloneID_cloneconsensus[data_cloneID_cloneconsensus[, "clonotype_id"]==clone,]
    TRs<-sort(paste(data_cloneID_cloneconsensus_clone[, "Consensus_annotations_v_gene"],data_cloneID_cloneconsensus_clone[, "Consensus_annotations_d_gene"]
               ,data_cloneID_cloneconsensus_clone[, "Consensus_annotations_j_gene"], data_cloneID_cloneconsensus_clone[, "Consensus_annotations_c_gene"]
               , data_cloneID_cloneconsensus_clone[, "Filtered_contig_annotations_cdr3_nt"], sep=":"))
    TRsString<-"";
    for(tr in TRs){
      TRsString<-paste(TRsString, tr, sep="#")
    }
    allClonesInfo<-c(allClonesInfo, TRsString)
    ##The clone number is calculated with record number divided by chain number
    CloneChainCount<-c(CloneChainCount, length(TRs))
    allClonesCount<-c(allClonesCount, dim(dataAnnotationCombined[dataAnnotationCombined[, "clonotype_id"]==clone,])[1]/length(TRs))
  }
  dataOut<-data.frame(allClonesInfo=allClonesInfo,allClonesCount=allClonesCount, CloneChainCount=CloneChainCount)
  rownames(dataOut)<-allClones
  dataOut<-dataOut[order(dataOut[, "allClonesCount"], decreasing = T),]
  write.csv(dataOut, file=paste("10xgenomics/10xgenomics_", sample, "_summary.csv", sep=""))
  ####Let us write the 10xgenomics cell information
  dataFilteredOut<-cbind(dataFiltered, dataOut[as.character(dataFiltered[, "raw_clonotype_id"]),])
  write.csv(dataFilteredOut, file=paste("10xgenomics/10xgenomics_", sample, ".csv", sep=""))
}##end for(sample)

