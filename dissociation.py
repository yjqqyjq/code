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
path="/Users/24756376/data/Flamingo/L1000N0900/"

halo_ids=fn.halo_ids
'''
S_200=np.zeros(len(halo_ids[halo_ids<=0]))
S_100=np.ones(len(halo_ids[halo_ids<=0]))
S_50=np.ones(len(halo_ids[halo_ids<=0]))
S_2=np.ones(len(halo_ids[halo_ids<=0]))
S_4=np.ones(len(halo_ids[halo_ids<=0]))
S_6=np.ones(len(halo_ids[halo_ids<=0]))
S_c_ub=np.ones(len(halo_ids[halo_ids<=0]))
qdm_m=np.ones(len(halo_ids[halo_ids<=0]))
qdm_c=np.ones(len(halo_ids[halo_ids<=0]))
qdm_r=np.ones(len(halo_ids[halo_ids<=0]))
qdm_rm=np.ones(len(halo_ids[halo_ids<=0]))
'''
f=h5py.File(path+'S_dist.hdf5','r')
peak=np.array(f["bar_peak"])

f.close()
bins=np.linspace(0.5,2.5,21)
bin=(bins[1:]+bins[:-1])/2
S_sat=np.zeros((len(halo_ids[halo_ids<=0]),20))
main_id=halo_ids[halo_ids<=0]
main_r200=fn.r50[halo_ids<=0]
for i in tqdm(range(len(S_sat))):
    
      particle=fn.load_regions(path,main_id[i].astype(int),5*main_r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="all")
      rdm=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      rg=np.sqrt(particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2)/main_r200[i]#the center is the center of the halo
      for j in range(20):
        S_sat[i][j]=fn.dissociation(particle[0][0][(rdm<bins[j])],particle[1][0][(rg<bins[j])],rdmg=1)#accumulate within r
      

'''
    particle_200=fn.load_regions(path,main_id[i].astype(int),fn.r200[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle_100=fn.load_regions(path,main_id[i].astype(int),fn.r100[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle_50=fn.load_regions(path,main_id[i].astype(int),fn.r50[fn.halo_ids<=0][main_id[i].astype(int)],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle_2=fn.load_regions(path,main_id[i].astype(int),2,dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle_4=fn.load_regions(path,main_id[i].astype(int),4,dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    particle_6=fn.load_regions(path,main_id[i].astype(int),6,dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})
    S_100[i]=fn.dissociation(particle_100[0][0],particle_100[1][0],rdmg=1)
    S_200[i]=fn.dissociation(particle_200[0][0],particle_200[1][0],rdmg=1)
    S_50[i]=fn.dissociation(particle_50[0][0],particle_50[1][0],rdmg=1)
    S_2[i]=fn.dissociation(particle_2[0][0],particle_2[1][0],rdmg=1)
    S_4[i]=fn.dissociation(particle_4[0][0],particle_4[1][0],rdmg=1)
    S_6[i]=fn.dissociation(particle_6[0][0],particle_6[1][0],rdmg=1)
            
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
    '''

#analyse the main halo

#calculate the dissociation

#h=np.histogram(S[S>-10],bins=100)
#print((h[1][1:]+h[1][-1:])/2)

f=h5py.File(path+'S_compare.hdf5','a')
del f["S_distance"]
f.create_dataset("S_distance",data=S_sat)
#f.create_dataset("S_r50",data=S_50)
#f.create_dataset("S_r100",data=S_100)
#f.create_dataset("S_r200",data=S_200)
#f.create_dataset("S_r2",data=S_2)
#f.create_dataset("S_r4",data=S_4)
#f.create_dataset("S_r6",data=S_6)
f.close()

'''

fig=plt.figure()
ax=fig.add_subplot(111)
ax.hist(S_f,bins=100,histtype='step')
ax.set_yscale('log')
plt.savefig("/Users/24756376/plot/Flamingo/L1000N0900/S_Falco.png")
'''