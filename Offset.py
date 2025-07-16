# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
import numpy as np
import unyt
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
O=np.zeros(len(halo_ids[halo_ids<=0]))
#Scale=np.zeros(len(halo_ids[halo_ids<=0]))
main_id=halo_ids[halo_ids<=0]

for i in tqdm(range(len(O))):
  if np.sum(fn.cross[(halo_ids>i)*(halo_ids<i+1)+(halo_ids==-i)])==0:
    center=fn.centers[fn.halo_ids<=0][i]
    particle=fn.load_particles(path,main_id[i],dm=1,g=0,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
    
    particle=[particle[0][0]-center,particle[1][0]-center]
#    r=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)
    
  
    O[i]=fn.offset(particle[0],particle[1])
  else:
    O[i]=-1
 
    

O=O/fn.r200[halo_ids<=0]
#analyse the main halo

print(main_id[O<0.007])

 
#plot
import matplotlib.pyplot as plt
plt.close()

fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(O[O!=-1], bins=10**np.linspace(-4,0.5,100))
ax.set_xlabel("Offset/r200")
ax.set_ylabel("Counts")
ax.set_title("M>10^14,Offsets between DM and Star CoM")
ax.set_yscale("log")
ax.set_xscale("log")
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/Offsets_DM_star.png")