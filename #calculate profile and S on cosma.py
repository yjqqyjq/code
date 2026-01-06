#calculate profile and S on cosma
import numpy as np
import unyt


import matplotlib.pyplot as plt

import h5py
from tqdm import tqdm

import subprocess
import sys

# Run the first script
#subprocess.run(["python", "functions.py"]) 
#work with SOAP
#load the data
def radial_distance(x, y,z):
    return np.sqrt(x**2 + y**2+z**2)

def spherical_harmonic_0(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(5/16/np.pi)*(3*z**2-r**2)/r**2
def spherical_harmonic_1(x, y,z):
    r = radial_distance(x, y,z)  
    return np.sqrt(15/4/np.pi)*x*z/r**2
def spherical_harmonic_2(x, y,z): 
    r = radial_distance(x, y,z)
    return np.sqrt(15/16/np.pi)*(x**2-y**2)/r**2
def spherical_harmonic__1(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(15/4/np.pi)*y*z/r**2 
def spherical_harmonic__2(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(15/4/np.pi)*x*y/r**2 
def quadrupole(x, y,z):
    r = radial_distance(x, y,z)
    
    f0=np.sum(spherical_harmonic_0(x, y,z)*r) / len(r)
    f1=np.sum(spherical_harmonic_1(x, y,z)*r) / len(r)
    f2=np.sum(spherical_harmonic_2(x, y,z)*r) / len(r)
    f_1=np.sum(spherical_harmonic__1(x, y,z)*r) / len(r)
    f_2=np.sum(spherical_harmonic__2(x, y,z)*r) / len(r)
    return np.sqrt(f0**2+f1**2+f2**2+f_1**2+f_2**2)
def dissociation(Coord_dm ,Coord_g,rdmg):
    x_dm, y_dm, z_dm = Coord_dm[:,0], Coord_dm[:,1], Coord_dm[:,2]
    x_g, y_g, z_g = Coord_g[:,0], Coord_g[:,1], Coord_g[:,2]
    #center of mass
    if len(x_dm)==0 or len(x_g)==0:
        return 0
    else:
       
        xdmc,ydmc,zdmc= np.average(x_dm), np.average(y_dm), np.average(z_dm)
    
        xgc,ygc,zgc= np.average(x_g), np.average(y_g), np.average(z_g)
        xc=(xdmc*rdmg+xgc)/(rdmg+1)
        yc=(ydmc*rdmg+ygc)/(rdmg+1)
        zc=(zdmc*rdmg+zgc)/(rdmg+1)
        xdm=np.array(x_dm-xc)   
        ydm=np.array(y_dm-yc)
        zdm=np.array(z_dm-zc)
        xg=np.array(x_g-xc)
        yg=np.array(y_g-yc)
        zg=np.array(z_g-zc)
        n_dm = len(xdm)
        n_g = len(xg)
    
        r_mean_dm = np.average(radial_distance(xdm, ydm,zdm))
        q_dm = quadrupole(xdm, ydm,zdm)
   
        r_mean_g = np.average(radial_distance(xg, yg,zg))
        q_g = quadrupole(xg, yg, zg)
   
        r_max = max(r_mean_dm, r_mean_g)
   
    
        return np.sqrt(4*np.pi/5) * (q_dm - q_g) / r_max

    return Coord
def load_properties(PartType,entry,k):
    f=h5py.File(path+'k.hdf5','r')
    if PartType=='dm':
     index=np.asarray(f["index_dm"],dtype=np.int32)
     N_dm=np.asarray(f["index_dm"])[:,1]
    f.close()
    f=h5py.File(snap+'k.hdf5','r')
    if PartType=='dm':
        data_all=np.asarray(f["PartType1"][entry],dtype=f["PartType1"][entry][0].dtype)[index]
    data=[]
    for i in range(len(N_dm)):
        i_i=np.sum(N_dm[0:i]) 
        i_f=np.sum(N_dm[0:i+1]) 
        data.append(data_all[i_i:i_f])
    return data
def load_regions_all(id,dm,g,s):#id is negative
   arg=int(id)
   
   
   slide=[]
   if dm==1:
     
   
     dm_s=int(np.sum(N_dm_region[0:-arg]))
     dm_e=int(np.sum(N_dm_region[0:-arg+1]))
     slide.append(slice(dm_s,dm_e))
    
   #  print(arg)
   if g==1:
      g_s=int(np.sum(N_g_region[0:-arg]))
      g_e=int(np.sum(N_g_region[0:-arg+1]))
      slide.append(slice(g_s,g_e))
     
   if s==1:
      s_s=int(np.sum(N_s_region[0:-arg]))
      s_e=int(np.sum(N_s_region[0:-arg+1]))
      
      slide.append(slice(s_s,s_e))

   slide=np.array(slide)
   key=[]
   if dm==1:
       key.append('dm')
   if g==1:
       key.append('gas')
   if s==1:
       key.append('stars')
   return key,slide


def load_regions(path,id,radius,dm=0,g=0,s=0):
   if id>0:#correct id input to 
      id=-int(id) 
   else:
      id=int(id)
   center=centers[-int(id)]
  
   keys,slides=load_regions_all(id,dm,g,s)
   f=h5py.File(path+'particles_region_ranked.hdf5','r')
      
  
      
   
   # create the slice and then find all the particles in the slice, that save the space by avoiding loading everything
   dataset=[]
   
   if dm==1:
      
      dataset.append(f['PartType1'])
   

      
   if g==1:
      
      dataset.append(f['PartType0'])
      

   if s==1:
       
      dataset.append(f['PartType2'])
  
  
    #load entry keys

   particles=[]

   for i in range(0,len(dataset)):
      comp=[]
      data=dataset[i]
      Coord=np.array(data['Coordinates'][slides[i]])
#      
      
     
      comp.append(Coord)
      
      
      
      
      
      particles.append(comp)#particle in shape dm, g, s
   f.close()      
   
   return particles
path="/Users/24756376/data/Flamingo/L1000N1800/"

f=h5py.File(path+'halos_ranked.hdf5','r')

N_dm_region=np.array(f['N_dm_region'])
N_g_region=np.array(f['N_g_region'])
N_s_region=np.array(f['N_s_region'])



mass=np.array(f['mass'])
#N_s_region=np.zeros(len(N_dm_region))

m50=np.array(f['m50'])
m200=np.array(f['m200'])
centers=np.array([f["centers"]])
r50=np.array(f["r50"])

r200=np.array(f["r200"])

f.close()


bins=np.linspace(0,1,11)
bin=(bins[1:]+bins[:-1])/2

S_50=np.zeros(len(mass))
profile_dm=np.zeros((len(mass),len(bins)-1))
for i in tqdm(range(len(S_50))):
      particle=load_regions(path,i.astype(int),r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]})


      
      S_50[i]=dissociation(particle[0][0],particle[1][0],rdmg=1)

      


f=h5py.File(path+'S_compare.hdf5','a')


f.close()