#count particles
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
hdm=[]
hgas=[]

hdm_sat=[]
hgas_sat=[]

hdm_neigh=[]
hgas_neigh=[]

for i in tqdm(range(len(main_id))):
    
  
    particle_5200=fn.load_regions(path,main_id[i].astype(int),5*r200[i],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
    particle_ub=fn.load_regions(path,main_id[i].astype(int),5*r200[i],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="unbound")
    particle_c=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
    particle_m=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
    particle_n=[[np.concatenate([particle_m[0][0],particle_ub[0][0]])],[np.concatenate([particle_m[1][0],particle_ub[1][0]])]]
    particle_s=[[np.concatenate([particle_c[0][0],particle_ub[0][0]])],[np.concatenate([particle_c[1][0],particle_ub[1][0]])]]
    rdm=np.sqrt(particle_5200[0][0][:,0]**2+particle_5200[0][0][:,1]**2+particle_5200[0][0][:,2]**2)
    rdm=rdm/r200[i]#in r200
    rg=np.sqrt(particle_5200[1][0][:,0]**2+particle_5200[1][0][:,1]**2+particle_5200[1][0][:,2]**2)
   
    rg=rg/r200[i]#in r200
    
  
    hdm.append(np.histogram(rdm,bins=bins)[0])
    
    hgas.append(np.histogram(rg,bins=bins)[0])

    rdm=np.sqrt(particle_n[0][0][:,0]**2+particle_n[0][0][:,1]**2+particle_n[0][0][:,2]**2)
    rdm=rdm/r200[i]#in r200
    rg=np.sqrt(particle_n[1][0][:,0]**2+particle_n[1][0][:,1]**2+particle_n[1][0][:,2]**2)
    rg=rg/r200[i]#in r200
    hdm_neigh.append(np.histogram(rdm,bins=bins)[0])
    hgas_neigh.append(np.histogram(rg,bins=bins)[0])
    rdm=np.sqrt(particle_s[0][0][:,0]**2+particle_s[0][0][:,1]**2+particle_s[0][0][:,2]**2)
    rdm=rdm/r200[i]#in r200
    rg=np.sqrt(particle_s[1][0][:,0]**2+particle_s[1][0][:,1]**2+particle_s[1][0][:,2]**2)
    rg=rg/r200[i]#in r200
    hdm_sat.append(np.histogram(rdm,bins=bins)[0])
    hgas_sat.append(np.histogram(rg,bins=bins)[0])
hdm=np.array(hdm)
hgas=np.array(hgas)
hdm_sat=np.array(hdm_sat)
hgas_sat=np.array(hgas_sat)
hdm_neigh=np.array(hdm_neigh)
hgas_neigh=np.array(hgas_neigh)

hgas_min=np.min(hgas, axis=1)
hgas_neigh_min=np.min(hgas_neigh, axis=1)
hgas_sat_min=np.min(hgas_sat, axis=1)
f=h5py.File(path+"particle_counts.hdf5","w")
f.create_dataset("hdm",data=hdm)
f.create_dataset("hgas",data=hgas)
f.create_dataset("hdm_sat",data=hdm_sat)#without neighbour or satellites
f.create_dataset("hgas_sat",data=hgas_sat)
f.create_dataset("hdm_neigh",data=hdm_neigh)#withought neighbour
f.create_dataset("hgas_neigh",data=hgas_neigh)
f.create_dataset("bin",data=bin)
f.create_dataset("hgas_min",data=hgas_min)
f.create_dataset("hgas_neigh_min",data=hgas_neigh_min)
f.create_dataset("hgas_sat_min",data=hgas_sat_min)
f.close()