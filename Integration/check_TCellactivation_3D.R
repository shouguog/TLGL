setwd("H:/Bradley/Integration/DetailedAnalysis_round1/DiffusionMap")
rm(list=ls())
library(Seurat)
library(scatterplot3d)
rgl_init <- function(new.device = FALSE, bg = "white", width = 1640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 1)
}
library(rgl)



###load all data
load("../../datasets/S.combined.RData")
Tactivation<-c("CD69","CCR7","CD27","BTLA","CD40LG","IL2RA","CD3E","CD47","EOMES","GNLY","GZMA","GZMB","PRF1","IFNG"
,"CD8A","CD8B","CD95L","LAMP1","LAG3","CTLA4","HLA-DRA","TNFRSF4","ICOS","TNFRSF9","TNFRSF18")
Tactivation<-intersect(Tactivation, rownames(S.combined@raw.data))
TCRinfomration<-read.csv("../../datasets/filtered_contig_annotations_all.csv")
TCRinfomrationTRB<-TCRinfomration[TCRinfomration$chain=="TRB",]
##Filter with expression
TCRinfomrationTRBexp<-TCRinfomrationTRB[TCRinfomrationTRB$umis>=2,]


seeds<-c(111)

for(seed1 in seeds){
  ##read diffusion result
  data<-read.csv(paste("sva_diffusion_map_", seed1, "_1200cells.csv", sep=""), header=T, row.names = 1)
  dataExp<-t(as.matrix(S.combined@raw.data[Tactivation,rownames(data)]))
  ###Let us try the average
  colors<-rep("blue", dim(data)[1])
  colors[rowSums(dataExp)>10]<-"red"
  #sum(rownames(data)%in%TCRinfomrationTRB$expressionBarcode)
  colorsTCR<-rep("blue", dim(data)[1])
  colorsTCR[rownames(data)%in%TCRinfomrationTRBexp$expressionBarcode]<-"red"
  
  png(paste("DM", seed1, "allSum10_1200_TactivationTCR2D.png", sep=""), width=2000, height=1000, res=200)
  par(mfrow=c(1,2))
  plot(data[,1], data[,2], pch = 16, col=colors,xlab="D1",ylab="D2", xlim=c(-0.5,0.5), ylim=c(-0.8,0.3), cex=0.3)
  plot(data[,1], data[,2], pch = 16, col=colorsTCR,xlab="D1",ylab="D2", xlim=c(-0.5,0.5), ylim=c(-0.8,0.3), cex=0.3)
  dev.off()
}



for(seed1 in seeds){
  ##read diffusion result
  data<-read.csv(paste("sva_diffusion_map_", seed1, "_1200cells.csv", sep=""), header=T, row.names = 1)
  dataExp<-t(as.matrix(S.combined@raw.data[Tactivation,rownames(data)]))
  ###Let us try the average
  colors<-rep("blue", dim(data)[1])
  colors[rowSums(dataExp)>10]<-"red"
  png(paste("DM", seed1, "allSum10_1200_Tactivation3D.png", sep=""), width=1600, height=1600, res=200)
  scatterplot3d(data[,c(1,2,3)], pch = 16, color=colors,xlab="D1",ylab="D2",zlab="D3", angle=45
                , col.axis="black", col.grid="white", axis=TRUE, xlim=c(-0.5,0.5), ylim=c(-0.8,0.8), zlim=c(-0.3,0.5))
  dev.off()
  ####change the data to rgl to limit 
  for(ii in 1:dim(data)[1]){
    if(data[ii,1]>0.5){data[ii,1]<-0.5}
    if(data[ii,1]<(-0.5)){data[ii,1]<-(-0.5)}
    if(data[ii,2]>0.8){data[ii,2]<-0.8}
    if(data[ii,2]<(-0.8)){data[ii,2]<-(-0.8)}
    if(data[ii,3]>0.5){data[ii,3]<-0.5}
    if(data[ii,3]<(-0.3)){data[ii,3]<-(-0.3)}
  }
  for(theta in 1:10*20){
    for(phi in 1:10*20){
      rgl_init()
      rgl.spheres(data[,1], data[,2], data[,3], r = 0.007, color = colors) 
      rgl.bbox(color=c("#888888","black"), emission="#888888",specular="#888888", shininess=1, alpha=0.8 ) 
      rgl.viewpoint(theta = theta, phi = phi, zoom = 1)
      rgl.snapshot(filename = paste("anglecheck/plot_theta_", theta, "_phi_", phi, ".png", sep=""))
    }
  }
}
