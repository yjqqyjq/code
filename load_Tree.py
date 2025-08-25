#load Tree

import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc

dir="/Users/24756376/data/Flamingo/L1000N0900/Trees/"
'''
files= glob.glob(dir+"*dm.joblib")
dm_all=[]
for file in tqdm(files):
   index,D=joblib.load(file)
   dm=[]
   for i in range(len(index)):
      dm.append(np.array([index[i], D[i]]))
   dm_all.append(dm)
dm_all2=[]
for i in range(0,len(dm_all[0])):#halos
   halo=[]
   for j in range(len(dm_all)):#files
        halo.append(dm_all[j][i])
   halo=np.concatenate(halo, axis=1)
   dm_all2.append(halo)


files= glob.glob(dir+"*g.joblib")
g_all=[]
for file in tqdm(files):
   index,D=joblib.load(file)
   g=[]
   for i in range(len(index)):
      g.append(np.array([index[i], D[i]]))
   g_all.append(g)
g_all2=[]
for i in range(0,len(g_all[0])):#halos
   halo=[]
   for j in range(len(g_all)):#files
        halo.append(g_all[j][i])
   halo=np.concatenate(halo, axis=1)
   g_all2.append(halo)  


joblib.dump(dm_all2,"/Users/24756376/data/Flamingo/L1000N0900/dm.joblib")
joblib.dump(g_all2,"/Users/24756376/data/Flamingo/L1000N0900/g.joblib")
'''
f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/halos_ranked.hdf5", 'r')
r200=np.array(f["r200"])
ids=np.array(f['id'])
r200=r200[ids<=0]
f.close()
dm_all2=joblib.load("/Users/24756376/data/Flamingo/L1000N0900/dm.joblib")
g_all2=joblib.load("/Users/24756376/data/Flamingo/L1000N0900/g.joblib")
bins=10**np.linspace(-1.5, 0.5, 21)
bin=10**np.linspace(-1.45, 0.45, 20)
hist_dm=np.zeros((len(dm_all2), len(bins)-1))
hist_g=np.zeros((len(g_all2), len(bins)-1))

for i in tqdm(range(len(dm_all2))):
    hist, _=np.histogram(dm_all2[i][1], bins=bins*r200[i])
    hist_dm[i] = np.cumsum(hist)
    hist, _ = np.histogram(g_all2[i][1], bins=bins*r200[i])
    hist_g[i]=np.cumsum(hist)

import h5py
f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/profile.hdf5", 'w')
f.create_dataset("dm", data=hist_dm/4/np.pi*3/bin**3)
f.create_dataset("g", data=hist_g/4/np.pi*3/bin**3)
f.create_dataset("bins", data=bin)
f.close()


   
