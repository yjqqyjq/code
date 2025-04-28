import numpy as np
import unyt
import swiftsimio as sw
import functions as fn
import h5py
import tqdm
from PIL import Image#To plot image directly from pixle
dir="/mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
#dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"
data_h=sw.load(dir)

host_id=data_h.soap.host_halo_index#central halo=-1\
halo_id=np.arange(0,len(host_id),1)
mass=data_h.bound_subhalo.total_mass
radius=data_h.bound_subhalo.enclose_radius
xc=data_h.bound_subhalo.centre_of_mass[:,0]
yc=data_h.bound_subhalo.centre_of_mass[:,1]
zc=data_h.bound_subhalo.centre_of_mass[:,2]
star_lumz=data_h.bound_subhalo.stellar_luminosity[:,4]
mainhalo_id=halo_id[(host_id==-1)*(mass>10000)]
BCGid=np.zeros(len(mainhalo_id))
BCGoffset=np.zeros(len(mainhalo_id))
i=0
for id in tqdm.tqdm(mainhalo_id):
  id=int(id)
  sub_id=halo_id[host_id==id]
  mask=np.isin(halo_id,sub_id)
  star_lum_sub=star_lumz[mask]
  if len(star_lum_sub)==0:
    BCGoffset[i]=-0.11#the halo doesn't have any subhalo
    BCGid[i]=id
    i=i+1
    continue
  
#  hosthalo is not included in the subhalo
  if star_lumz[id]>=np.max(star_lum_sub):
    BCGoffset[i]=0
    BCGid[i]=id
  else:
    BCG_id=sub_id[star_lum_sub==np.max(star_lum_sub)]
    BCGid[i]=BCG_id
    BCGoffset[i]=fn.radial_distance(xc[BCG_id]-xc[id],yc[BCG_id]-yc[id],
                                  zc[BCG_id]-zc[id])/float(radius[id])
  i=i+1
np.savetxt("/home/jyang/data/Flamingo/L1000N0900/M>1e14_BCG_halo_offset",np.array([mainhalo_id,BCGid,BCGoffset]),header=
           "check whether the BCG is in the main halo, data stored in mainhalo_id, BCG halo_id, the r between to the COM of BCG and the main halo, normalized by radius",comments="#")
print(BCGoffset[BCGoffset>0])
'''
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
title="offset between CoM of BCG halo and the main halo"
import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlabel("BCGoffset/r_halo")
ax.set_ylabel("Counts")
ax.set_title(title)
#for i in range(len(mainhalo_id)):
#    ax.plot(bin,rho_g[i],color='b')
ax.hist(BCGoffset,bins=20)
ax.set_yscale("log")
#ax.plot(np.logspace(0,4,100),np.logspace(0,4,100)*0.1,color='r')
fig.savefig("/home/jyang/plot/Flamingo/L1000N0900/BCG_offset_CoM_halo.png")
