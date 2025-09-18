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

dir='/Users/24756376/data/Flamingo/L1000N0900/halos/'

for i in tqdm(range(0,100)):
  f=h5py.File(dir+str(-i)+"_no_sat.hdf5",'w')
  dm=f.create_group("PartType1")
  g=f.create_group("PartType0")

  fu=h5py.File(dir+str(-i)+"_1.0_r100_unbound.hdf5",'r')
  fc=h5py.File(dir+str(-i)+".hdf5",'r')
  

  for key in list(fu["PartType1"].keys()):
 
    data=np.concatenate([np.array(fu["PartType1"][key]),np.array(fc["PartType1"][key])])
    dm.create_dataset(key,data=data)
  for key in list(fu["PartType0"].keys()):
    data=np.concatenate([np.array(fu["PartType0"][key]),np.array(fc["PartType0"][key])])
    
    g.create_dataset(key,data=data)
  fu.close()
  fc.close()
  f.close()
    