#metalicity
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
path="/Users/24756376/data/Flamingo/L1000N1800/"

halo_ids=fn.halo_ids


bins=np.linspace(0,2,21)
bin=(bins[1:]+bins[:-1])/2
binsm=np.linspace(-10,-0.5,16)
binm=(binsm[1:]+binsm[:-1])/2
main_id=halo_ids[halo_ids<=0]
main_r200=fn.r50[halo_ids<=0]

dm=np.zeros((len(halo_ids[halo_ids<=0]),len(bin),len(binm)))
g=np.zeros((len(halo_ids[halo_ids<=0]),len(bin),len(binm)))
vdm=np.zeros(20)

for i in tqdm(range((len(halo_ids[halo_ids<=0])))):
      particle=fn.load_regions(path,main_id[i].astype(int),2*main_r200[i],dm=0,g=1,s=0,coordinate=1,extra_entry={"dm":["Velocities"],"gas":["MetalMassFractions"],"stars":[]},mode="all")


      rg=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      
      Mg=particle[0][1]
      Mg[Mg<10**-10]=10**-10
#      print(np.min(Mg[Mg>10**-15]))
  
      
     
  
 
      hg=np.histogram2d(rg, Mg, bins=[bins,10**binsm])[0]
   
      
      
      g[i]=hg

f=h5py.File(path+'2d_Z.hdf5','w')



f.create_dataset("Ng",data=g)
f.create_dataset("bin",data=bin)
f.create_dataset("binm",data=binm)



f.close()