#Xray Lum
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
Xray=np.zeros(len(halo_ids[halo_ids<=0]))
M=np.zeros(len(halo_ids[halo_ids<=0]))
N_sat=np.zeros(len(halo_ids[halo_ids<=0]))
main_id=halo_ids[halo_ids<=0]
for i in tqdm(range(0,len(main_id))):
  if np.sum(fn.cross[(halo_ids>i)*(halo_ids<i+1)+(halo_ids==-i)])==0:
    particle_all=fn.load_particles(path,main_id[i],dm=0,g=1,s=0,coordinate=0,extra_entry={"dm":[],"gas":["xray_lum_erosita_low"],"stars":[]},mode="cluster")
    particle_cen=fn.load_particles(path,main_id[i],dm=0,g=1,s=0,coordinate=0,extra_entry={"dm":[],"gas":["xray_lum_erosita_low"],"stars":[]},mode="halo")
    
    M[i]=np.sum(fn.mass[(halo_ids>=i)*(halo_ids<i+1)+(halo_ids==-i)])
    N_sat[i]=len(fn.mass[(halo_ids>=i)*(halo_ids<i+1)+(halo_ids==-i)])
    Xray[i]=(np.sum(particle_all[0][0])-np.sum(particle_cen[0][0]))/np.sum(particle_all[0][0])
#    if i<5:
#      print(M[i])
N_sat2=np.where(N_sat>40,40,N_sat)
print(main_id[(Xray>10**-5)*(Xray<10**-4)])
#M=fn.mass[halo_ids<=0]
plt.figure()
ax=plt.subplot(111)
s=ax.scatter(M,Xray,s=0.5,alpha=0.3,c=N_sat2)
plt.colorbar(s,label="N_sat")
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.set_xlabel('Mass(central+sattllite)/10^10Msun')
ax.set_ylabel('f_sat')
ax.set_title('X-ray contribution f_satellite')
plt.savefig('/Users/24756376/plot/Flamingo/L1000N0900/xray_lum_ratio_test.png')