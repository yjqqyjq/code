#baryon_peak
import numpy as np
import unyt


import matplotlib.pyplot as plt
import functions as fn

import h5py
from tqdm import tqdm

import subprocess
import sys

# Run the first script
#subprocess.run(["python", "functions.py"]) 
#work with SOAP
#load the data
path="/Users/24756376/data/Flamingo/L1000N0900/"

halo_ids=fn.halo_ids
main_id=halo_ids[halo_ids<=0]
r200=fn.r50[halo_ids<=0]
print(r200[0],fn.r200[halo_ids<=0][0])
bins=np.linspace(0,2.5,26)
bin=(bins[1:]+bins[:-1])/2
density=np.zeros(len(main_id))
f_bar=[]
f_hot=[]
f=h5py.File(path+'neighbour.hdf5','r')

index=np.array(f["index"])
print(index)
f.close()
for i in tqdm(range(len(main_id))):
    
  
    particle_5200=fn.load_regions(path,main_id[i].astype(int),5*r200[i],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":["Temperatures"],"stars":[]},mode="all")
    rdm=np.sqrt(particle_5200[0][0][:,0]**2+particle_5200[0][0][:,1]**2+particle_5200[0][0][:,2]**2)
    rdm=rdm[rdm>0*r200[i]]/r200[i]#in r200
    rg=np.sqrt(particle_5200[1][0][:,0]**2+particle_5200[1][0][:,1]**2+particle_5200[1][0][:,2]**2)
    T=particle_5200[1][1][rg>0*r200[i]]#in K
    rg=rg[rg>0*r200[i]]/r200[i]#in r200
    
    rs=np.sqrt(particle_5200[2][0][:,0]**2+particle_5200[2][0][:,1]**2+particle_5200[2][0][:,2]**2)
    rs=rs[rs>0*r200[i]]/r200[i]#in r200 
    hdm=np.histogram(rdm,bins=bins)[0]
    
    hgas=np.histogram(rg,bins=bins)[0]
    hstar=np.histogram(rs,bins=bins)[0]
    hdm=np.cumsum(hdm)
    hgas=np.cumsum(hgas)
    hstar=np.cumsum(hstar)
  
#    print(hgas)
    fbar=(hgas*8.56+hstar*6.714)/(hgas*8.56+hstar*6.714+hdm*45.2)/(0.0494/0.316)#in cosmic average
    
    f_bar.append(fbar)
    hg_hot=np.histogram(rg[T>10**5],bins=bins)[0]
    hg_hot=np.cumsum(hg_hot)
    fbar_hot=(hg_hot*8.56)/(hstar*6.714+hgas*8.56+hdm*45.2)/(0.0494/0.316)
   
    f_hot.append(fbar_hot)
f=h5py.File(path+'S_dist.hdf5','a')
del f["f_bar_cu"]
del f["f_hot_cu"]
f.create_dataset("f_bar_cu",data=np.array(f_bar))
f.create_dataset("f_hot_cu",data=np.array(f_hot))

f.close()