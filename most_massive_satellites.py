#most_massive_satellites
import numpy as np
import unyt

#import Falco_S as fd
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
main_r200=fn.r50[halo_ids<=0]
mms=np.zeros(len(halo_ids[halo_ids<=0]))
centers=fn.centers
r50=fn.r50
bins=np.linspace(0,2.5,26)
bin=(bins[1:]+bins[:-1])/2
for i in tqdm(range(len(mms))):
#    particle=fn.load_regions(path,main_id[i].astype(int),main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
#    particle_ub=fn.load_regions(path,main_id[i].astype(int),main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="unbound")
    
#    N_dm_r50=len(particle[0][0])
#    N_g_r50=len(particle[1][0])   
    
    arg_mms=int(np.argwhere(halo_ids==-i))+1
    if halo_ids[arg_mms]<0:
#        N_dm_mms=0
#        N_g_mms=0
       mms[i]=0
    else:
#        N_dm_mms=fn.N_dm[arg_mms]
#        N_g_mms=fn.N_g[arg_mms]
       dist=centers[arg_mms]-centers[arg_mms-1]
       
       dist=np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)#(N_dm_mms*45.2+N_g_mms*8.56)/(N_dm_r50*45.2+N_g_r50*8.56) 
       if dist>10:
         print(dist,r50[arg_mms-1])
#    mms[i]=(len(particle_ub[0][0])*45.2+len(particle_ub[1][0])*8.56)/(N_dm_r50*45.2+N_g_r50*8.56) 
f=h5py.File(path+"satellites.hdf5",'a')
f.create_dataset("mms_dist",data=mms)
f.close()