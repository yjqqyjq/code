import faiss
import h5py
import glob
from tqdm import tqdm
import numpy as np
from scipy import spatial
f=h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/1.hdf5', 'r')
p=np.array(f["PartType1"]["Coordinates"])
f.close()
index= faiss.IndexFlatL2(3)


index.add(p)
#for i in tqdm(range(0,1000)):
print("start")  
dm=index.range_search(np.zeros((1000,3)), 5**2)[0]
print(dm,len(p))
Tree= spatial.KDTree(p)
for i in tqdm(range(0,1000)):
  dm2=Tree.query_ball_point(np.array([[0,0,0]]), 5)[0]
  