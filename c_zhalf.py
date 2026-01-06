import numpy as np
import unyt
import scipy
#import Falco_S as fd
import matplotlib.pyplot as plt
import functions_13 as fn
from scipy.interpolate import interp1d
import h5py
from tqdm import tqdm

import subprocess
import sys

# Run the first script
#subprocess.run(["python", "functions.py"]) 
#work with SOAP
#load the data
path="/Users/24756376/data/Flamingo/L1000N1800/"




#cdm=np.zeros(50000)
#cg=np.zeros(50000)

vdm=np.zeros(20)

f=h5py.File(path+"mass_evolution_12.hdf5","r")
masses=np.array(f["m200"])

f.close()

#masses=masses[0:19]
masses=(masses/masses[-1]).T

main_id=fn.halo_ids[fn.halo_ids<=0]
main_r200=fn.r200[fn.halo_ids<=0]


bins=np.linspace(0,1,31)
z=np.arange(1,2.5,0.05)[::-1]
cg=np.zeros(len(masses))
z_half=np.zeros(len(masses))
for i in tqdm(range(len(masses))):
      
      particle=fn.load_regions(path,main_id[i].astype(int),main_r200[i],dm=0,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
      '''
   
      particle_b=fn.load_particles(path,main_id[i].astype(int),dm=0,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
      if len(particle_b[0][0])==0:
            print("error")
      for j in range(0,2):
          particle[j][0]=np.concatenate([particle_b[j][0],particle[j][0]])
        
      
#      rdm=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      rg=np.sqrt(particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2)/main_r200[i]
      rs=np.sqrt(particle[2][0][:,0]**2+particle[2][0][:,1]**2+particle[2][0][:,2]**2)/main_r200[i]
      rg=rg[rg<1]
      rs=rs[rs<1]     
#      rdm=np.sort(rdm)
      rg=np.sort(rg)
      rs=np.sort(rs)
      hg=np.histogram(rg, bins=bins)[0]
      hs=np.histogram(rs, bins=bins)[0] 
      m=hg+hs*6.174/8.56
      m=np.cumsum(m)/np.sum(m)
      cg[i]=interp1d(m, bins[1:])(0.5)
#      j=i-len(halo_ids[halo_ids<=0])+10000
    
     

     
      
      '''
      j=i
      if np.min(masses[i])==0:
          z_half[j]=-1
      elif len(np.argwhere(masses[i]<0.5))==0:
          z_half[j]=np.log10(1+z[0])
      
      else:
          z_half[j]=scipy.interpolate.interp1d(masses[i],np.log10(1+z))(0.5)
          
    
    
          
#      if z_half[j]==0:
#          print(z[0],arg,z[arg])
    
      '''

f=h5py.File(path+"c_z_half_13.hdf5","a")


#del f["cdm_200_25"]
del f["cg_200"]
#f.create_dataset("cdm_200_25",data=cdm)
f.create_dataset("cg_200",data=cg)
#del f["log_z_half_200"]
#f.create_dataset("log_z_half_200",data=z_half)
f.close()