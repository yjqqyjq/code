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
bins=np.linspace(0.5,2.5,21)
bin=(bins[1:]+bins[:-1])/2
density=np.zeros(len(main_id))
f=h5py.File(path+'neighbour.hdf5','r')

index=np.array(f["index"])
print(index)
f.close()
for i in tqdm(range(len(main_id))):
    if index[i]>0:
        continue
    
  
    particle_5200=fn.load_regions(path,main_id[i].astype(int),5*r200[i],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
    rdm=np.sqrt(particle_5200[0][0][:,0]**2+particle_5200[0][0][:,1]**2+particle_5200[0][0][:,2]**2)
    rdm=rdm[rdm>0.5*r200[i]]/r200[i]#in r200
    rg=np.sqrt(particle_5200[1][0][:,0]**2+particle_5200[1][0][:,1]**2+particle_5200[1][0][:,2]**2)
    rg=rg[rg>0.5*r200[i]]/r200[i]#in r200
    rs=np.sqrt(particle_5200[2][0][:,0]**2+particle_5200[2][0][:,1]**2+particle_5200[2][0][:,2]**2)
    rs=rs[rs>0.5*r200[i]]/r200[i]#in r200 
    
    fbar=(np.histogram(rg,bins=bins)[0]*8.56+np.histogram(rs,bins=bins)[0]*6.714)/(np.histogram(rdm,bins=bins)[0])*45.2/(0.02222/0.7**2/0.316)#in cosmic average
    max=bin[np.argmax(fbar)]
    density[i]=max
f=h5py.File(path+'S_dist.hdf5','a')
del f["bar_peak"]
f.create_dataset("bar_peak",data=density)

f.close()