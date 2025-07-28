#mass_T_relation
import numpy as np
import unyt


import functions as fn

import h5py
from tqdm import tqdm
import matplotlib.pyplot as plt


# Run the first script
#subprocess.run(["python", "functions.py"]) 
#work with SOAP
#load the data
path="/Users/24756376/data/Flamingo/L1000N0900/"

halo_ids=fn.halo_ids

M=np.zeros(len(halo_ids[halo_ids<=0]))
Tmean=np.zeros(len(halo_ids[halo_ids<=0]))
main_id=halo_ids[halo_ids<=0]
for i in tqdm(range(0,len(main_id))):
  if np.sum(fn.cross[(halo_ids>i)*(halo_ids<i+1)+(halo_ids==-i)])==0:
    particle_all=fn.load_particles(path,main_id[i],dm=0,g=1,s=0,coordinate=0,extra_entry={"dm":[],"gas":["T"],"stars":[]},mode="cluster")
   
    
    M[i]=np.sum(fn.mass[(halo_ids>=i)*(halo_ids<i+1)+(halo_ids==-i)])
  
    Tmean[i]=np.average(particle_all[0][0])
#    if i<5:
#      print(M[i])


#M=fn.mass[halo_ids<=0]
plt.figure()
ax=plt.subplot(111)
s=ax.scatter(M,Tmean,s=0.5,alpha=0.3)
plt.colorbar(s,label="N_sat")
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.set_xlabel('Mass(central+sattllite)/10^10Msun')
ax.set_ylabel('T_mean')
#ax.set_title('X-ray contribution f_satellite')
plt.savefig('/Users/24756376/plot/Flamingo/L1000N0900/M_T_relation.png')