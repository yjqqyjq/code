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
S=np.zeros(len(halo_ids[halo_ids<=0]))
Scale=np.zeros(len(halo_ids[halo_ids<=0]))
main_id=halo_ids[halo_ids<=0]
for i in tqdm(range(len(S))):
  if np.sum(fn.cross[(halo_ids>i)*(halo_ids<i+1)+(halo_ids==-i)])==0:
    particle=fn.load_particles(path,main_id[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
  
  
    S[i]=fn.dissociation(particle[0][0],particle[1][0])
  else:
    S[i]=-10


#analyse the main halo

#calculate the dissociation
'''
S=np.zeros(len(x_dm))exit
for i in  range(0,len(x_dm)):
      S[i]=fn.dissociation(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
      Offset[i]=fn.offset(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
'''   

print(main_id[(S>0.3)*(S<0.4)])
 
#plot
import matplotlib.pyplot as plt
plt.close()

fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(S[S>=-1], bins=100)
ax.set_xlabel("Diossociation")
ax.set_ylabel("Counts")
ax.set_title("M>10^14")
ax.set_yscale("log")
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/Dissociation.png")

'''
plt.close()
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.scatter(S[S>-10],Scale[S>-10],s=0.5,alpha=0.3)
ax.set_xlabel("Dissociations")
ax.set_ylabel("fgas")
ax.set_title("M>10^4")
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/fgas.png")
'''