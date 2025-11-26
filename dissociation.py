# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
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


bins=np.linspace(0,1,11)
bin=(bins[1:]+bins[:-1])/2

main_id=halo_ids[halo_ids<=0]
main_r200=fn.r50[halo_ids<=0]
S_dist=np.zeros((len(halo_ids[halo_ids<=0]),len(bin)))
S_50=np.zeros(len(halo_ids[halo_ids<=0]))
S_rm=np.zeros(len(halo_ids[halo_ids<=0]))
for i in tqdm(range(len(S_50))):
      particle=fn.load_regions(path,main_id[i].astype(int),2.5*main_r200[i],dm=0,g=0,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
#      particle_ub=fn.load_regions(path,main_id[i].astype(int),1*main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="unbound")
   
#      particle_c=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
#      rdm=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      
#      rg=np.sqrt(particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2)/main_r200[i]
#      particle_c=[[particle_c[0][0][rdm<1]],[particle_c[1][0][rg<1]]]

#      particle_m=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
      
#      particle_r=fn.load_regions(path,main_id[i].astype(int),2*main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
#      particle=particle_ub
    
      rs=np.sqrt(np.sum(particle[0][0]**2,axis=1))/main_r200[i]
      rs=np.sort(rs)
      print(np.percentile(rs,90))


#      S_50[i]=fn.dissociation(particle[0][0],particle[1][0],rdmg=1)

      
#      for j in range(10):
#        S_dist[i][j]=fn.dissociation(particle[0][0][(rdm<bins[j+1])],particle[1][0][(rg<bins[j+1])],rdmg=1)#accumulate within r
#        S_nei[i][j]=fn.dissociation(particle_m[0][0][(rdm_m<bins[j+1])],particle_m[1][0][(rg_m<bins[j+1])],rdmg=1)#accumulate within r


   

            

#      S_rm[i]=fn.dissociation(particle_m[0][0],particle_m[1][0],rdmg=1)


#analyse the main halo

#calculate the dissociation

#h=np.histogram(S[S>-10],bins=100)
#print((h[1][1:]+h[1][-1:])/2)

f=h5py.File(path+'S_compare.hdf5','a')

del f["S_50"]
#del f["S_c_ub"]
#f.create_dataset("S_dist",data=S_dist)
#f.create_dataset("S_c_ub",data=S_rm)

#f.create_dataset("S_g",data=S_g)
f.create_dataset("S_50",data=S_50)
#f.create_dataset("S_r100",data=S_100)
#f.create_dataset("S_r200",data=S_200)

f.close()


