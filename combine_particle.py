#combine_particle

import numpy as np
import h5py
import tracemalloc
from pathlib import Path
from tqdm import tqdm
import datetime

from multiprocessing import Pool
from functools import partial
import os
print("start")
print(datetime.datetime.now())
keys_dm=["Coordinates","Velocities","member"]
keys_s=["Coordinates","Velocities","Luminosities","member"]
keys_g=["Coordinates","Velocities","LastAGNFeedbackScaleFactors","xray_lum_erosita_high","xray_lum_erosita_low","member"]
keys_bh=["Coordinates","member"]
dir='/home/jyang/data/Flamingo/L1000N0900/cluster_particles'
f=h5py.File("/home/jyang/data/Flamingo/L1000N0900/cluster_particles.hdf5",'w')
dm=f.create_group("PartType1")
g=f.create_group("PartType0")
s=f.create_group("PartType2")
bh=f.create_group("PartType3")
for key in keys_dm:
    data=[]
    for i in range(0,32):
        subfile=h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(i)+".hdf5",'r')
        data.append(np.array(subfile["PartType1"][key],dtype=subfile["PartType1"][key][0].dtype))
        subfile.close()
    data=np.concatenate(data)
    dm.create_dataset(key,data=data)
print("dm Done")
for key in keys_g:
    data=[]
    for i in range(0,32):
        subfile=h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(i)+".hdf5",'r')
        data.append(np.array(subfile["PartType0"][key],dtype=subfile["PartType0"][key][0].dtype))
        subfile.close()
    data=np.concatenate(data)
    g.create_dataset(key,data=data)
print("gas Done")
for key in keys_s:
    data=[]
    for i in range(0,32):
        subfile=h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(i)+".hdf5",'r')
        data.append(np.array(subfile["PartType2"][key],dtype=subfile["PartType2"][key][0].dtype))
        subfile.close()
    data=np.concatenate(data)
    s.create_dataset(key,data=data)
for key in keys_bh:
    data=[]
    for i in range(0,32):
        subfile=h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(i)+".hdf5",'r')
        data.append(np.array(subfile["PartType0"][key],dtype=subfile["PartType0"][key][0].dtype))
        subfile.close()
    data=np.concatenate(data)
    bh.create_dataset(key,data=data)
f.close()