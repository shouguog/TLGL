###This is used to calculate /data/gaos2/phenograph/PhenoGraph-master
### 
import pandas as pd
import numpy as np
dataAll=pd.read_csv("metaInfor_sva_limma_tSNE.csv", index_col=0)
data=dataAll[["tSNE1_combat","tSNE2_combat"]].as_matrix()
#"tSNE1_combat","tSNE2_combat","tSNE1_limma","tSNE2_limma"
import phenograph
communities, graph, Q = phenograph.cluster(data)
file=open('community_combat.txt','w')
for ii in range(len(communities)):
	file.write(dataAll.index[ii])
	file.write("\t")
	file.write(str(ii))
	file.write("\t")
	file.write(str(communities[ii]))
	file.write("\n")
file.close()
data=dataAll[["tSNE1_limma","tSNE2_limma"]].as_matrix()
#"tSNE1_combat","tSNE2_combat","tSNE1_limma","tSNE2_limma"
import phenograph
communities, graph, Q = phenograph.cluster(data)
file=open('community_limma.txt','w')
for ii in range(len(communities)):
	file.write(dataAll.index[ii])
	file.write("\t")
	file.write(str(ii))
	file.write("\t")
	file.write(str(communities[ii]))
	file.write("\n")
file.close()

