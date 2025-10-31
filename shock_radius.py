#shock radius
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
path="/Users/24756376/data/Flamingo/L1000N1800_NoCool/"

halo_ids=fn.halo_ids


main_id=halo_ids[halo_ids<=0]
main_r50=fn.r50[halo_ids<=0]
mms=np.zeros(len(halo_ids[halo_ids<=0]))
centers=fn.centers
Rs=np.zeros(len(halo_ids[halo_ids<=0]))

for i in tqdm(range(len(mms))):
    particle=fn.load_regions(path,main_id[i].astype(int),2.5*main_r50[i],dm=0,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":["Temperatures","ElectronNumberDensities"],"stars":[]},mode="all")
    bins=np.linspace(0,2.5*main_r50[i],26)
    bin=(bins[1:]+bins[:-1])/2
    r=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)
    T=particle[0][1]*8.6*10**(-5)#k in ev
    n_e=particle[0][2]/2.938/10**67# (cubic meters)
   
    Rbin=np.digitize(r, bins)-1
    K_m=0
    im=0
    for j in range(len(bin)):
        
        ave_T=np.average(T[Rbin==j])
        ave_n_e=np.average(n_e[Rbin==j])
        K=ave_T/(ave_n_e**(2/3))
        
        
        if K>K_m:
            K_m=K
            im=j
    Rs[i]=bin[im] 
print(Rs)   
        
#    mms[i]=(len(particle_ub[0][0])*45.2+len(particle_ub[1][0])*8.56)/(N_dm_r50*45.2+N_g_r50*8.56) 
f=h5py.File(path+"particle_counts.hdf5",'a')
del f["Shock_Radius"]
f.create_dataset("Shock_Radius",data=Rs)
#print(np.array(f["Shock_Radius"]))
f.close()