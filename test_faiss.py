import faiss
import h5py
import glob
from tqdm import tqdm
import numpy as np
from scipy import spatial
from sklearn.neighbors import KDTree
import datetime
f=h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/1.hdf5', 'r')
pa=np.array(f["PartType1"]["Coordinates"])
p=[]
for i in range(0,100):
  p.append(pa)
f.close()
p=np.concatenate(p)
print(len(p))#number of data points
index= faiss.IndexFlatL2(3)


index.add(p)
#for i in tqdm(range(0,1000)):
print("start") 
faiss.omp_set_num_threads(5)
#for i in tqdm(range(0,1000)): 
print(datetime.datetime.now())
dm=index.range_search(np.zeros((1000,3)), 1**2)[0]
print(datetime.datetime.now())

Tree= KDTree(p)
print(datetime.datetime.now())
#for i in tqdm(range(0,1000)):
dm2=Tree.query_radius(np.zeros((1000,3)), 1,count_only=True)[0]
#  print(len(dm2))
print(datetime.datetime.now())