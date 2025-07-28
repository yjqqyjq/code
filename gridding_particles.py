#gridding the pparticles
import numpy as np
import h5py
import tracemalloc
from tqdm import tqdm
import sys
import itertools
import multiprocessing as mp 
from functools import partial
import datetime
from collections import Counter
import os

print(datetime.datetime.now())
path="/Users/24756376/data/Flamingo/L1000N0900/"

f=h5py.File(path+'cluster_particles0.0_unbound.hdf5','r')
Coord_dm=np.array(f['PartType1']['Coordinates'],dtype=np.float32)
#Coord_g=np.array(f['PartType0']['coordinates'],dtype=np.int32)
#Coord_s=np.array(f['PartType2']['coordinates'],dtype=np.int32)
f.close()
mask=(Coord_dm[:,0]<600)*((Coord_dm[:,0]>500))*(Coord_dm[:,1]<900)*((Coord_dm[:,1]>800))*(Coord_dm[:,2]<850)*(Coord_dm[:,2]>800)
Coord_dm_new=Coord_dm[mask]
f=h5py.File(path+'cluster_particles0.8333333333333333_unbound.hdf5','r')
Coord_dm=np.array(f['PartType1']['Coordinates'],dtype=np.float32)

f.close()
mask=(Coord_dm[:,0]<600)*((Coord_dm[:,0]>500))*(Coord_dm[:,1]<900)*((Coord_dm[:,1]>800))*(Coord_dm[:,2]<850)*(Coord_dm[:,2]>800)
Coord_dm_new2=Coord_dm[mask]
Coord_dm_s=np.concatenate([Coord_dm_new,Coord_dm_new2],axis=0)
Coord_dm_s[:,0]-=550
Coord_dm_s[:,1]-=850
Coord_dm_s[:,2]-=825
f=h5py.File(path+'slice_dm.hdf5','w')
dm=f.create_group("PartType0")
dm.create_dataset("Coordinates",data=Coord_dm_s)
f.close()