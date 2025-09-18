import numpy as np
import unyt
import swiftsimio as sw

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
mass=fn.mass
R=np.zeros(len(halo_ids))
for i in tqdm(range(0,len(halo_ids))):
    particle=fn.load_particles(path,halo_ids[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
    R[i]=len(particle[1][0])/len(particle[0][0])
import matplotlib.pyplot as plt
plt.close()

fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(mass[halo_ids<=0],R[halo_ids<=0],color='r',s=0.3,alpha=0.3,label="central")
ax.scatter(mass[halo_ids>0],R[halo_ids>0],color='b',s=0.3,alpha=0.3,label="sat")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
ax.set_xlabel("M")
ax.set_ylabel("f_gas")

#ax.set_yscale("log")
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/f_gas.png")