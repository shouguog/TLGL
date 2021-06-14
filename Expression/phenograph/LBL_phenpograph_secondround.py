###This is used to calculate /data/gaos2/phenograph/PhenoGraph-master
### 
import pandas as pd
import numpy as np
dataAll=pd.read_csv("locCluster_combat.csv", index_col=0)
data=dataAll[["tSNE1_mean","tSNE2_mean"]].as_matrix()
import phenograph
communities, graph, Q = phenograph.cluster(data, k=5, min_cluster_size=3)
file=open('community_combat_secondround.txt','w')
for ii in range(len(communities)):
	file.write(dataAll.index[ii])
	file.write("\t")
	file.write(str(ii))
	file.write("\t")
	file.write(str(communities[ii]))
	file.write("\n")
file.close()
dataAll=pd.read_csv("locCluster_limma.csv", index_col=0)
data=dataAll[["tSNE1_mean","tSNE2_mean"]].as_matrix()
import phenograph
communities, graph, Q = phenograph.cluster(data, k=5, min_cluster_size=3)
file=open('community_limma_secondround.txt','w')
for ii in range(len(communities)):
	file.write(dataAll.index[ii])
	file.write("\t")
	file.write(str(ii))
	file.write("\t")
	file.write(str(communities[ii]))
	file.write("\n")
file.close()

