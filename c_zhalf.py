import numpy as np
import unyt
import scipy
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
path="/Users/24756376/data/Flamingo/L1000N1800/"

f=h5py.File(path+"c_z_half_13.hdf5","a")

f.close()

halo_ids=fn.halo_ids


main_id=halo_ids[halo_ids<=0]
main_r200=fn.r50[halo_ids<=0]
m50=fn.m50[halo_ids<=0]
cdm=np.zeros(50000)
cg=np.zeros(50000)

vdm=np.zeros(20)
f=h5py.File(path+"mass_evolution_13.hdf5","r")
masses=np.array(f["m50"])

f.close()
#masses=masses[0:19]
masses=(masses/masses[-1]).T




z=np.arange(0,1.5,0.1)[::-1]
print(len(z),len(masses[0]),len(masses[1]))
z_half=np.zeros(len(m50))
for i in tqdm(range(50000)):
      
      particle=fn.load_regions(path,main_id[i].astype(int),main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
   
      
      particle_b=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
      if len(particle_b[0][0])==0:
            print("error")
      for j in range(0,2):
          particle[j][0]=np.concatenate([particle_b[j][0],particle[j][0]])
        
      
      rdm=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      rg=np.sqrt(particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2)/main_r200[i]
      rs=np.sqrt(particle[2][0][:,0]**2+particle[2][0][:,1]**2+particle[2][0][:,2]**2)/main_r200[i]
      rdm=np.sort(rdm)
      rg=np.sort(rg)
      rs=np.sort(rs)

#      j=i-len(halo_ids[halo_ids<=0])+10000
    
     
      cdm[i]=np.percentile(rdm,25)
      cg[i]=np.percentile(rg,25)
     
      
      '''
      j=i
      if np.min(masses[i])==0 or np.max(masses[i])>1.01:
          z_half[j]=-1
      elif len(np.argwhere(masses[i]<0.5))==0:
          z_half[j]=z[0]
      
      else:
          z_half[j]=scipy.interpolate.interp1d(masses[i],z)(0.5)
          
          arg=np.argwhere(masses[i]<0.5)[-1]
          if 0.5-masses[i][arg]>masses[i][arg+1]-0.5:
              arg=arg+1
          z_half[j]=z[arg]
     '''
          
#      if z_half[j]==0:
#          print(z[0],arg,z[arg])
    


f=h5py.File(path+"c_z_half_13.hdf5","a")


#del f["cdm_200_25"]
#del f["cg_200_25"]
f.create_dataset("cdm_200_25",data=cdm)
f.create_dataset("cg_200_25",data=cg)
#del f["z_half_200_com"]
#f.create_dataset("z_half_200_1",data=z_half)
f.close()