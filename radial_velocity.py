#radial velocity
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

H0=68#0
bins=np.linspace(0,2.5,21)
bin=(bins[1:]+bins[:-1])/2
binsv=np.linspace(-3,3,21)
binv=(binsv[1:]+binsv[:-1])/2
main_id=halo_ids[halo_ids<=0]
main_r200=fn.r50[halo_ids<=0]
m50=fn.m50[halo_ids<=0]
Sa=np.zeros(len(halo_ids[halo_ids<=0]))
dm=np.zeros((len(halo_ids[halo_ids<=0]),len(bin),len(binv)))
gs=np.zeros((len(halo_ids[halo_ids<=0]),len(bin),len(binv)))
gt=np.zeros((len(halo_ids[halo_ids<=0]),len(bin),len(binv)))
g=np.zeros((len(halo_ids[halo_ids<=0]),len(bin),len(binv)))
print(bin)
#vdm=np.zeros(20)
'''
for i in tqdm(range((len(halo_ids[halo_ids<=0])))):
     
      particle=fn.load_regions(path,main_id[i].astype(int),2.5*main_r200[i],dm=0,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":["Temperatures","ElectronNumberDensities"],"stars":[]})
      S=particle[0][1]*8.6*10**(-5)/((particle[0][2]/2.938/10**67)**(2/3))
      rg=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      S[particle[0][2]==0]=0
      S=S/0.7754/m50[i]**(2/3)
      hg=np.histogram(rg, bins=bins,weights=S)[0]/np.histogram(rg, bins=bins)[0]
      Sa[i]=np.argmax(hg)
print(np.histogram(Sa))    
 
'''
for i in tqdm(range((len(halo_ids[halo_ids<=0])))):
     
      particle=fn.load_regions(path,main_id[i].astype(int),2.5*main_r200[i],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":["Velocities"],"gas":["Velocities"],"stars":["Velocities"]})
    
#      particle_b=fn.load_particles(path,main_id[i].astype(int),dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":["Velocities"],"gas":["Velocities","Temperatures","ElectronNumberDensities"],"stars":[]},mode="cluster")
#      if len(particle[0][0])==0 or len(particle_b[0][0])==0:
#            print("error")
#      for j in range(0,2):
#          particle[j][0]=np.concatenate([particle_b[j][0],particle[j][0]])
#          particle[j][1]=np.concatenate([particle_b[j][1],particle[j][1]])
#      particle[1][2]=np.concatenate([particle_b[1][2],particle[1][2]])
#      particle[1][3]=np.concatenate([particle_b[1][3],particle[1][3]])
      rdm=np.sqrt(particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2)/main_r200[i]
      rg=np.sqrt(particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2)/main_r200[i]
      rs=np.sqrt(particle[2][0][:,0]**2+particle[2][0][:,1]**2+particle[2][0][:,2]**2)/main_r200[i]
#      T=particle[1][2]/np.average(particle[1][2][rg<0.1])/36000/m50[i]**(2/3)#Normalized by T50
   
#      X=particle[1][3][:,0]/np.average(particle[1][3][:,0][rg<0.1])
#      S=particle[1][2]*8.6*10**(-5)/((particle[1][3]/2.938/10**67)**(2/3))
#      S[particle[1][3]==0]=0
#      S=S/0.7754/m50[i]**(2/3)#From Appendix Joey et al, normalized by K50
#      print(np.min(S))
      vmean=np.average(np.concatenate([45.2/8.56*particle[0][1][rdm<0.1],particle[1][1][rg<0.1]],axis=0),axis=0)/(45.2/8.56+1)#mean velocity within 0.1 r50 as COM velocity
#      print(vmean)
#      vdm=np.sum(particle[0][0] * (particle[0][1]-vmean), axis=1)/np.linalg.norm(particle[0][0], axis=1)
#      vdm=(vdm+H0*rdm*main_r200[i])/(np.sqrt(6.558*m50[i]/main_r200[i]))#normolized by V50, corrected by Hubble flow
      vs=np.sum(particle[2][0] * (particle[2][1]-vmean), axis=1)/np.linalg.norm(particle[2][0], axis=1)
      vs=(vs+H0*rs*main_r200[i])/(np.sqrt(6.558*m50[i]/main_r200[i]))#normolized by V50, corrected by Hubble flow
#      vg=np.sum(particle[1][0] * (particle[1][1]-vmean), axis=1)/np.linalg.norm(particle[1][0], axis=1)
#      vg=(vg+H0*rg*main_r200[i])/(np.sqrt(6.558*m50[i]/main_r200[i]))
      
#      print(np.sqrt(6.558*m50[i]/main_r200[i]),main_r200[i])
#      hdm=np.histogram2d(rdm, vdm, bins=[bins,binsv])[0]
#      dm[i]=hdm
      hg=np.histogram2d(rs, vs, bins=[bins,binsv])[0]
      g[i]=hg
#      hg=np.histogram2d(rg, vg, bins=[bins,binsv],weights=S)[0]/(np.histogram2d(rg, vg, bins=[bins,binsv])[0]+10**-10)
      
#      print(np.sum(hg,axis=1))
#      gs[i]=hg
#      hg=np.histogram2d(rg, vg, bins=[bins,binsv],weights=T)[0]/(np.histogram2d(rg, vg, bins=[bins,binsv])[0]+10**-10)
#      gt[i]=hg
#      hg=np.histogram2d(rg, vg, bins=[bins,binsv])[0]
#      g[i]=hg
     

     
f=h5py.File(path+'2d_fgas_H.hdf5','a')


f.create_dataset("Ns",data=g)
#f.create_dataset("bin",data=bin)
#f.create_dataset("binv",data=binv)



f.close()
'''
    v=np.sqrt(np.sum(particle[0][1]**2,axis=1))
  
      h=np.histogram(v,bins=np.linspace(0,4000,21))[0]
   
      vdm+=h
vdm/=1000

fig=plt.figure()
ax=fig.add_subplot(111)
bin=np.linspace(0,4000,21)
ax.plot((bin[1:]+bin[:-1])/2,vdm)
#ax.set_yscale('log')
fig.savefig('/Users/24756376/plot/Flamingo/L1000N1800/vg_distribution.png')
'''