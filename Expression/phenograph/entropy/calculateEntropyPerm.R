rm(list=ls())
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
png("entropy_permutation.png", width=1000, height=1000, res=150)
hist(permData)
dev.off()




