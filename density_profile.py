import numpy as np
import unyt
from velociraptor import load
import functions as fn
import h5py
import tqdm
dir="/mnt/su3-pro/flamingo/L0200N0360/"

data_h=load(dir+"VR/halos_0008.properties.0")

mass=data_h.masses.mass_tot
is_mainhalo=data_h.centrals
mask=(is_mainhalo*(mass>10000))
radius=data_h.radii.r_200crit
m200=data_h.masses.mass_200crit
ms=data_h.masses.mass_star
mainhalo_id=np.array(data_h.ids.id[mask])-1
mainhalo_id=mainhalo_id.astype(int)
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles/"
'''
bin=np.zeros(100)
for i in range(100):
  bin[i]=i/100+0.005
rho_g=[]

for id in tqdm.tqdm(mainhalo_id):
  f=h5py.File(path+str(id)+'.hdf5','r')
  Coord_g=f['PartType1']["Coordinates"]
  xyz_g=np.array(Coord_g).T
  r=float(radius[id])
  x_g=xyz_g[0]
  y_g=xyz_g[1]
  z_g=xyz_g[2]
  r_g=fn.radial_distance(x_g,y_g,z_g)/r
  print(np.max(r_g))
  f.close()
  hist=np.histogram(r_g, bins=100)
  rho=hist[0]/bin**2/float(mass[id])
  rho_g.append(rho)
'''
import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlabel("Mh")
ax.set_ylabel("Ms")
ax.set_yscale('log')
ax.set_xscale('log')
#for i in range(len(mainhalo_id)):
#    ax.plot(bin,rho_g[i],color='b')
ax.scatter(ms,ms/m200,s=0.01,marker='o',color='b')
#ax.plot(np.logspace(0,4,100),np.logspace(0,4,100)*0.1,color='r')
fig.savefig("/home/jyang/plot/Flamingo/L0200N0360/mh_ms_relation.png")
