#number_profile
#radial velocity
import numpy as np
import unyt

#import Falco_S as fd
import matplotlib.pyplot as plt
import functions_13 as fn

import h5py
from tqdm import tqdm

import subprocess
import sys

# Run the first script
#subprocess.run(["python", "functions.py"]) 
#work with SOAP
#load the data
path="/Users/24756376/data/Flamingo/L1000N1800_NoCool/"

halo_ids=fn.halo_ids

H0=68#0
bins=10**np.linspace(-1.5,np.log10(2.5),25)
bin=(bins[1:]+bins[:-1])/2
binsv=np.linspace(-3,3,21)
binv=(binsv[1:]+binsv[:-1])/2
main_id=halo_ids[halo_ids<=0][0:10000]
main_r200=fn.r50[halo_ids<=0][0:10000]
m50=fn.m50[halo_ids<=0][0:10000]
Sa=np.zeros(len(main_id))
dm=np.zeros((len(main_id),len(bin)))
gs=np.zeros((len(main_id),len(bin)))
s=np.zeros((len(main_id),len(bin)))
g=np.zeros((len(main_id),len(bin)))
print(bin)
#vdm=np.zeros(20)

for i in tqdm(range((len(main_id)))):
     
      particle=fn.load_regions(path,main_id[i].astype(int),2.5*main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
      particle_b=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
      for j in range(0,2):
          particle[j][0]=np.concatenate([particle_b[j][0],particle[j][0]])
          
      rdm=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      rg=np.sqrt(particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2)/main_r200[i]
#      rs=np.sqrt(particle[2][0][:,0]**2+particle[2][0][:,1]**2+particle[2][0][:,2]**2)/main_r200[i]
   
      hdm=np.histogram(rdm,bins=bins)[0]
      dm[i]=hdm
  
      hg=np.histogram(rg, bins=bins)[0]
      g[i]=hg
#      hs=np.histogram(rs, bins=bins)[0]
#      s=np.zeros(len(bin))
     

     
f=h5py.File(path+'radial_profile_13.5.hdf5','w')
f.create_dataset("Ndm",data=dm)
f.create_dataset("Ng",data=g)

#f.create_dataset("Ns",data=s)

f.close()