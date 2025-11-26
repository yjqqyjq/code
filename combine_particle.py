import numpy as np
import h5py
import tracemalloc
from pathlib import Path
from tqdm import tqdm
import datetime
from functools import partial
from multiprocessing import Pool

import os
def load_dm(i,key):
   
    subfile=h5py.File(dir+'cluster14/particles_radius_add'+str(i)+'.hdf5','r')
    sub=np.array(subfile["PartType1"][key],dtype=subfile["PartType1"][key][0].dtype)
        
    subfile.close()
    return sub
def load_g(i,key):
    
    subfile=h5py.File(dir+'cluster14/particles_radius_add'+str(i)+'.hdf5','r')
    if key=="XrayLuminosities":
      sub=np.array(subfile["PartType0"][key],dtype=subfile["PartType0"][key][0].dtype).T
    else:
      sub=np.array(subfile["PartType0"][key],dtype=subfile["PartType0"][key][0].dtype)
    
    subfile.close()
    return sub
def load_s(i,key):
    subfile=h5py.File(dir+'cluster13/particles_radius_add'+str(i)+".hdf5",'r')
    sub=np.array(subfile["PartType2"][key],dtype=subfile["PartType2"][key][0].dtype)
     
    subfile.close()

    return sub

def load_bh(i,key):
    subfile=h5py.File(dir+'/clusters13/cluster_particles'+str(i)+".hdf5",'r')
    sub=np.array(subfile["PartType3"][key],dtype=subfile["PartType3"][key][0].dtype)
     
    subfile.close()
    return sub
if __name__ == '__main__':
  print("start")
  print(datetime.datetime.now())
  keys_dm=[]#["Coordinates","Velocities"]
  keys_s=["Coordinates","Velocities","Luminosities","member"]
  keys_g=["Temperatures","XrayLuminosities"]#"Coordinates"]#["Velocities","Temperatures","ElectronNumberDensities"]#,"MetalMassFractions"]#
  keys_bh=["Coordinates","member"]
  dir='/Users/24756376/data/Flamingo/L1000N1800_NoCool/'
  f=h5py.File(dir+"/particles_radius_add.hdf5",'w')
  dm=f.create_group("PartType1")
  g=f.create_group("PartType0")
  #s=f.create_group("PartType2")
  #bh=f.create_group("PartType3")
  
  m=h5py.File(dir+"/cluster14/mask.hdf5",'r')
#  mask_dm=np.array(m["mask_dm_ub"],dtype=np.int64)
  mask_g=np.array(m["mask_g_ub"],dtype=np.int64)
  m.close()

  for key in keys_dm:
    print(key)
  
    data=[]
    '''
    for i in range(0,64):   
        data.append(load_dm(i,key)) 
#    with Pool(1) as p:
#      data=p.map(partial(load_dm, key=key), range(0,64))
    
    data=np.concatenate(data,dtype=data[0][0].dtype,axis=0)[mask_dm]
    dm.create_dataset(key,data=data)
  print("dm Done")
    '''
  for key in keys_g:
    print(key)
#    with Pool(6) as p:
#      data=p.map(partial(load_g, key=key), range(0,64))
    data=[]
    for i in range(0,64):   
        data.append(load_g(i,key))  
    
               
    data=np.concatenate(data,dtype=data[0][0].dtype,axis=0)[mask_g]
    g.create_dataset(key,data=data)
  print("gas Done")
  f.close()