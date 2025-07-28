import faiss
import h5py
import glob
from tqdm import tqdm
import numpy as np
dir="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(dir+'halos_ranked.hdf5','r')

centers=np.array([f["center_x"],f["center_y"],f["center_z"]]).T
r200=np.array(f["r200"])
ids=np.array(f['id'])
f.close()
r200=r200[ids<=0]

centers=centers[ids<=0]

f=h5py.File(dir+'particles_ranked.hdf5','r')
p=np.array(f["PartType1"]["Coordinates"],dtype=np.float32)

f.close()
for center in tqdm(centers):
    Coord=p-center
    Coord=Coord[(Coord[:,0]<2)*(Coord[:,0]>-2)*(Coord[:,1]<2)*(Coord[:,1]>-2)*(Coord[:,2]<2)*(Coord[:,2]>-2)]
    
    np.histogram((Coord[:,0]**2+Coord[:,1]**2+Coord[:,2]**2), bins=10, range= [0, 2])
   