# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio import load
import swiftgalaxy as sg
import functions as fn
from swiftgalaxy import SWIFTGalaxy, MaskCollection
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
S_m=np.ones(len(halo_ids[halo_ids<=0]))*10
S_c=np.ones(len(halo_ids[halo_ids<=0]))*10
S_r=np.ones(len(halo_ids[halo_ids<=0]))*10
S_r2=np.ones(len(halo_ids[halo_ids<=0]))*10
main_id=halo_ids[halo_ids<=0]
for i in tqdm(range(len(S_m))):
  if np.sum(fn.cross[(halo_ids>i)*(halo_ids<i+1)+(halo_ids==-i)])==0:
    particle=fn.load_particles(path,main_id[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
    particle_c=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
    particle_r=fn.load_regions(path,main_id[i].astype(int),fn.r200[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle_r2=fn.load_regions(path,main_id[i].astype(int),2*fn.r200[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    S_m[i]=fn.dissociation(particle[0][0],particle[1][0])
    S_c[i]=fn.dissociation(particle_c[0][0],particle_c[1][0])
    S_r[i]=fn.dissociation(particle_r[0][0],particle_r[1][0])
    S_r2[i]=fn.dissociation(particle_r2[0][0],particle_r2[1][0])
#  else:
#   S[i]=-10


#analyse the main halo

#calculate the dissociation

#h=np.histogram(S[S>-10],bins=100)
#print((h[1][1:]+h[1][-1:])/2)
f=h5py.File(path+'S_compare.hdf5','w')
f.create_dataset("central",data=S_c)
f.create_dataset("central_sat",data=S_m)
f.create_dataset("r200",data=S_r)
f.create_dataset("2r200",data=S_r2)
f.close()
'''
#plot
import matplotlib.pyplot as plt
plt.close()

fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(S_m,S_c,color='r',s=0.3,alpha=0.3,label="central")
ax.scatter(S_m,S_r,color='b',s=0.3,alpha=0.3,label="r200")
ax.scatter(S_m,S_r2,color='k',s=0.3,alpha=0.3,label="2r200")
ax.legend()
ax.set_xlabel("S_cen_sat")
ax.set_ylabel("S_other")

#ax.set_yscale("log")
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/Dissociation_compare.png")

plt.close()
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.scatter(S[S>-10],Scale[S>-10],s=0.5,alpha=0.3)
ax.set_xlabel("Dissociations")
ax.set_ylabel("fgas")
ax.set_title("M>10^4")
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/fgas.png")
'''