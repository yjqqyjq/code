#satellite_ratio
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
R_sat=np.zeros(len(main_id))
R_mms=np.zeros(len(main_id))
for i in tqdm(range(len(main_id))):
    particle=fn.load_particles(path,main_id[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
    particle_c=fn.load_particles(path,main_id[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
#    particle_mms=fn.load_particles(path,halo_ids[np.argwhere(halo_ids==-1)+1],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
    particle_r=fn.load_regions(path,main_id[i].astype(int),fn.r100[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
    rdm=particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2
    rg=particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2
  
    mask_dm=~np.isin(particle[0][0][:,0],particle_c[0][0][:,0])*(rdm<=fn.r100[halo_ids<=0][i]**2)
    mask_g=~np.isin(particle[1][0][:,0],particle_c[1][0][:,0])*(rg<=fn.r100[halo_ids<=0][i]**2)
    rdm=rdm[mask_dm]
    rg=rg[mask_g]
    R_sat[i]=(len(rdm)*45.2+len(rg)*8.56)/(len(particle_r[0][0])*45.2+len(particle_r[1][0])*8.56)
f=h5py.File(path+'S_compare.hdf5','a')
#f.create_dataset("S_r50",data=S_r)
f.create_dataset("R_sat_100",data=R_sat)


f.close()