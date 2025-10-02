#import unyt

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
#dir="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/flamingo_0077.hdf5"

#soap_dir="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
tracemalloc.start()
#dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/colibre_with_SOAP_membership_0127.hdf5"
#soap_dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"


         
def load_particles(k):
    """
    Load particles from a given directory and region.
    """
    
    
    
    data=h5py.File(folder_path+str(k)+".hdf5",'r')
    member=h5py.File(path_member+str(k)+".hdf5",'r')
    
    member_dm=np.array(member["PartType1"]["GroupNr_bound"],dtype=np.int32)
    
    pmask=np.isin(member_dm,sample_input_ids,assume_unique=True)
    
    member_dm=member_dm[pmask]
    f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(k)+'.hdf5', 'w')
    dm=f.create_group("PartType1")
    keys_dm=["Coordinates","Velocities"]
    for key in keys_dm:
        print(key)
        data_w=np.array(data["PartType1"][key],dtype=data["PartType1"][key][0].dtype)[pmask]
        dm.create_dataset(key,data=data_w)
        print(key)
    
    dm.create_dataset("member",data=member_dm)

#    potential_dm=np.array(data.dark_matter.potentials,dtype=np.float32)[pmask]
    
    print(tracemalloc.get_traced_memory())
    f.close()
    print("dm done")
    
    
    
    member_dm=[]
  
  
    member_star=np.array(member["PartType4"]["GroupNr_bound"],dtype=np.int32)
    pmask=np.isin(member_star,sample_input_ids,assume_unique=True)
    member_star=member_star[pmask]

    keys_s=["Coordinates","Velocities","Luminosities"]
    f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(k)+'.hdf5', 'a')
    s= f.create_group("PartType2")
    for key in keys_s:
        if key=="Luminosities":
            data_w=np.array(data["PartType4"][key][:,4],dtype=np.float32)[pmask]
            s.create_dataset(key,data=data_w)
        else:
            data_w=np.array(data["PartType4"][key],dtype=data["PartType4"][key][0].dtype)[pmask]
            s.create_dataset(key,data=data_w)
    s.create_dataset("member",data=member_star)
    f.close()
  
#    potential_star=potential_star[pmask2]

    print("star done")
    member_bh=np.array(member["PartType5"]["GroupNr_bound"],dtype=np.int32)
    pmask=np.isin(member_bh,sample_input_ids)
    Coord_bh=np.array(data["PartType5"]["Coordinates"][pmask])
#    potential_bh=np.array(data.black_holes.potentials,dtype=np.float32)[pmask]
    print("bh done")
    print(tracemalloc.get_traced_memory())
       
    f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(k)+'.hdf5', 'a')

    
    bh=f.create_group("PartType3")
    bh.create_dataset("Coordinates", data=Coord_bh)
#    bh.create_dataset("potentials", data=potential_bh)
    bh.create_dataset("member",data=member_bh)
    f.close()
    
   
    member_star=[]
    member_bh=[]

    member_gas=np.array(member["PartType0"]["GroupNr_bound"],dtype=np.int32)

    pmask=np.isin(member_gas,sample_input_ids,assume_unique=True)

    member_gas=member_gas[pmask]
    member.close()
  
 #   potential_gas=np.array(data.gas.potentials,dtype=np.float32)
   
    
    f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(k)+'.hdf5', 'a')
    g=f.create_group("PartType0")

    keys_g=["Coordinates","Velocities","LastAGNFeedbackScaleFactors","XrayLuminosities"]
#    g.create_dataset("potentials", data=potential_gas)
    for key in keys_g:
        if key=="XrayLuminosities":
          lum_low=np.array(data["PartType0"][key][:,1],dtype=np.float32)[pmask]
          lum_high=np.array(data["PartType0"][key][:,2],dtype=np.float32)[pmask]
          g.create_dataset("xray_lum_erosita_high", data=lum_low)
          g.create_dataset("xray_lum_erosita_low", data=lum_high)
          lum_low=[]
          lum_high=[]
        else:
          data_w=np.array(data["PartType0"][key],dtype=data["PartType0"][key][0].dtype)[pmask]
          g.create_dataset(key,data=data_w)
    g.create_dataset("member",data=member_gas)
    f.close()
    
 

#    potential_gas=potential_gas[pmask2]
    
    
    print(tracemalloc.get_traced_memory())
    print("gas done")

    

    data.close()
    print(datetime.datetime.now())
def append_particles(k):
    """
    Append particles from a given region.
    """
    data=h5py.File(folder_path+str(k)+".hdf5",'r')
    member=h5py.File(path_member+str(k)+".hdf5",'r')
    member_gas=np.array(member["PartType0"]["GroupNr_bound"],dtype=np.int32)

    pmask=np.isin(member_gas,sample_input_ids,assume_unique=True)

    member_gas=[]
    member.close()
  
 #   potential_gas=np.array(data.gas.potentials,dtype=np.float32)
   
    
    f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_14/cluster_particles'+str(k)+'.hdf5', 'a')
    g=f["PartType0"]

    keys_g=["Temperatures"]
#    g.create_dataset("potentials", data=potential_gas)
    for key in keys_g:
      
          
          data_w=np.array(data["PartType0"][key],dtype=data["PartType0"][key][0].dtype)[pmask]
          g.create_dataset(key,data=data_w)
   
    f.close()
    data.close()
    print(tracemalloc.get_traced_memory())


if __name__ == '__main__':
  data_h=h5py.File("/home/jyang/data/Flamingo/L1000N0900/halos_ranked.hdf5")

  sample_input_ids=np.array(data_h["input_ids"],dtype=np.int64)

  data_h.close()
  path_member="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/membership_0077/membership_0077."

  print("halo Done\n")

#de some sequence here

  folder_path ="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/snapshots/flamingo_0077/flamingo_0077."

  print("start")  
  with Pool() as p: # Create a pool with 5 worker processes
        p.map(append_particles, range(0,32))

#    append_particles(region[j],j)

     
#    load_unbound_particles(region[j],j)    
print(tracemalloc.get_traced_memory())
    




print("Done\n")
print(datetime.datetime.now())
#print(len(member_dm[member_dm!=-1]))
    

