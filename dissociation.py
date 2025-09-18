# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
import numpy as np
import unyt

import Falco_S as fd
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
path="/Users/24756376/data/Flamingo/L1000N0900/"

halo_ids=fn.halo_ids
S_f=np.zeros(len(halo_ids[halo_ids<=0]))
S_m=np.ones(len(halo_ids[halo_ids<=0]))
S_c=np.ones(len(halo_ids[halo_ids<=0]))
S_r=np.ones(len(halo_ids[halo_ids<=0]))
S_rm=np.ones(len(halo_ids[halo_ids<=0]))
S_r2=np.ones(len(halo_ids[halo_ids<=0]))
S_c_ub=np.ones(len(halo_ids[halo_ids<=0]))
qdm_m=np.ones(len(halo_ids[halo_ids<=0]))
qdm_c=np.ones(len(halo_ids[halo_ids<=0]))
qdm_r=np.ones(len(halo_ids[halo_ids<=0]))
qdm_rm=np.ones(len(halo_ids[halo_ids<=0]))
main_id=halo_ids[halo_ids<=0]
main_r100=fn.r100[halo_ids<=0]
for i in tqdm(range(len(S_m))):
  

        
#    particle_r2=fn.load_regions(path,main_id[i].astype(int),fn.r50[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle=fn.load_particles(path,main_id[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
    rdm=particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2
    rg=particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2
#    S_f[i]=fd.calculate_S(particle[1][0][rg<main_r200[i]**2], particle[0][0][rdm<main_r200[i]**2], np.ones(len(particle[1][0][rg<main_r200[i]**2])), np.ones(len(particle[0][0][rdm<main_r200[i]**2]))*45.2/8.56, dim = 3)



    particle_c=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
    
    


    particle_r=fn.load_regions(path,main_id[i].astype(int),fn.r100[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
    particle_r2=fn.load_regions(path,main_id[i].astype(int),fn.r200[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
    #unbound within r100
    particle_ub=fn.load_regions(path,main_id[i].astype(int),fn.r100[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="unbound")
    
   
  
    
    S_m[i]=fn.dissociation(particle[0][0],particle[1][0],rdmg=1)#Here, we select the geometry center of the DM and gas CoM, to see the effects of gemoetry
    S_c[i]=fn.dissociation(particle_c[0][0],particle_c[1][0],rdmg=1)
    S_r[i]=fn.dissociation(particle_r[0][0],particle_r[1][0],rdmg=1)
    S_rm[i]=fn.dissociation(particle[0][0][rdm<=main_r100[main_id[i].astype(int)]**2],particle[1][0][rg<=main_r100[main_id[i].astype(int)]**2],rdmg=1)
    S_r2[i]=fn.dissociation(particle_r2[0][0],particle_r2[1][0],rdmg=1)
#    print(len(particle_c[0][0]),len(particle_ub[0][0]))
    S_c_ub[i]=fn.dissociation(np.concatenate([particle_c[0][0],particle_ub[0][0]],axis=0),np.concatenate([particle_c[1][0],particle_ub[1][0]],axis=0),rdmg=1)
    qdm_m[i]=fn.dissociation(particle[0][0],np.array([[0,0,0.001],[0,0,0]]),rdmg=1)
    qdm_c[i]=fn.dissociation(particle_c[0][0],np.array([[0,0,0.001],[0,0,0]]),rdmg=1)
    qdm_r[i]=fn.dissociation(particle_r[0][0],np.array([[0,0,0.001],[0,0,0]]),rdmg=1)
    qdm_rm[i]=fn.dissociation(particle[0][0][rdm<main_r100[i]**2],np.array([[0,0,0.001],[0,0,0]]),rdmg=1)
    

#analyse the main halo

#calculate the dissociation

#h=np.histogram(S[S>-10],bins=100)
#print((h[1][1:]+h[1][-1:])/2)

f=h5py.File(path+'S_compare.hdf5','w')
#f.create_dataset("S_r50",data=S_r)
f.create_dataset("S_cen_unbound",data=S_c_ub)

f.create_dataset("q_central",data=qdm_c)
f.create_dataset("q_central_sat",data=qdm_m)
f.create_dataset("q_r100",data=qdm_r)
f.create_dataset("q_mem_r100",data=qdm_rm)
f.create_dataset("S_central",data=S_c)
f.create_dataset("S_central_sat",data=S_m)
f.create_dataset("S_r100",data=S_r)
f.create_dataset("S_mem_r100",data=S_rm)
f.create_dataset("S_r200",data=S_r2)
f.close()
'''
fig=plt.figure()
ax=fig.add_subplot(111)
ax.hist(S_f,bins=100,histtype='step')
ax.set_yscale('log')
plt.savefig("/Users/24756376/plot/Flamingo/L1000N0900/S_Falco.png")
'''