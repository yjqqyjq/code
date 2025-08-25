#Kdtree s
from sklearn.neighbors import KDTree
import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool
dir="/Users/24756376/data/Flamingo/L1000N0900/particles_ranked.hdf5"
data_dir="/Users/24756376/data/Flamingo/L1000N0900/Trees/"
'''
Coord=[]

file=dir
f=h5py.File(file, 'r')
Coord.append(np.array(f["PartType2"]["Coordinates"]))
f.close()

Coord=np.concatenate(Coord, axis=0)
print(len(Coord))    
tree = KDTree(Coord)
joblib.dump(tree, data_dir+"s.pkl")


Tree_s=joblib.load(data_dir+'s.pkl')
index,D=Tree_s.query_radius(
        centers,
        r=6,
        return_distance=True,
        sort_results=True,
        
    )
joblib.dump([index,D],'/Users/24756376/data/Flamingo/L1000N0900/Trees/s.joblib')

'''
f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/halos_ranked.hdf5", 'r')
centers=np.array([f['centers_x'],f['centers_y'],f['centers_z']]).T
ids=np.array(f['id'])
centers = centers[ids <= 0]
r200=np.array(f['r200'])
r200=r200[ids<=0]
print(len(centers))
f.close()
index,D=joblib.load("/Users/24756376/data/Flamingo/L1000N0900/Trees/s.joblib")
s_all=[]


for i in range(len(index)):
   
   s_all.append(np.array([index[i], D[i]]))

bins=10**np.linspace(-1.5, 0.5, 21)
bin=10**np.linspace(-1.45, 0.45, 20)
hist_s=np.zeros((len(s_all), len(bins)-1))


for i in tqdm(range(len(s_all))):
    hist, _ =np.histogram(s_all[i][1], bins=bins*r200[i])
    hist_s[i]= np.cumsum(hist)

import h5py
f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/profile.hdf5", 'a')
del f["s"]
f.create_dataset("s", data=hist_s/4/np.pi*3/np.diff(bins**3))


f.close()

