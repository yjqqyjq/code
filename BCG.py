import numpy as np
import unyt
import swiftsimio as sw
import functions as fn
import h5py
import tqdm
from PIL import Image#To plot image directly from pixle
#dir="/mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
#dir="../../../mnt/su3-pro/colibre/L0200N1504/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"




path="/home/jyang/data/Flamingo/L1000N0900/"
#path="/home/jyang/data/Colibre/L0200N1504/"
f=h5py.File(path+'halos.hdf5','r')

halo_id=np.array(f["halos"]["id"])

radius=np.array(f["halos"]["radius"])
host_id=np.array(f["halos"]["hostid"])
mass=np.array(f["halos"]["mass"])
mbp=np.array(f["halos"]["center"])
com_star=np.array(f["halos"]["com_star_100kpc"])
star_lumz=np.array(f["halos"]["lumz_100kpc"])
f.close()

#ignore all the halos without models of star lum and com
rs=np.sqrt(com_star[:,0]**2+com_star[:,1]**2+com_star[:,2]**2)
mask=(star_lumz>0)#*(rs>0)
halo_id=halo_id[mask]
host_id=host_id[mask]
mbp=mbp[mask]

com_star=com_star[mask]
star_lumz=star_lumz[mask]
mass=mass[mask]
radius=radius[mask]
ids=np.arange(0,len(halo_id),1)
mainhalo_id=ids[(host_id==-1)*(mass>10000)]
BCGid=np.zeros(len(mainhalo_id))
BCGoffset=np.zeros(len(mainhalo_id))
Rbri=np.ones(len(mainhalo_id))
'''
mbp=mbp[(host_id!=-1)]
com_star=com_star[(host_id!=-1)]
radius=radius[(host_id!=-1)]
mass=mass[(host_id!=-1)]
d=com_star-mbp
o=fn.radial_distance(d[:,0],d[:,1],d[:,2])
print(np.max(mass))
import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.subplot(111)
sc=ax.scatter(radius, o/radius,s=0.05,alpha=0.5,c=mass,cmap='rainbow',label="subhalos")
b=plt.colorbar(sc)
b.set_label("mass/10^10Msun")
ax.plot(np.arange(0.1,10,0.1),0.01/np.arange(0.1,10,0.1),'k',label="Softing")
ax.legend()
ax.set_xlabel("r_en/Mpc")
ax.set_ylabel("Offsets/r_en")
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig("/home/jyang/plot/Flamingo/L1000N0900/star_CoM_offset_sub.png")
'''
i=0

for id in tqdm.tqdm(mainhalo_id):
  id=int(id)

  sub_id=ids[host_id==halo_id[id]]
  mask=np.isin(ids,sub_id)
  star_lum_sub=star_lumz[mask]
 
  

  if len(sub_id)==0:
#    BCGoffset[i]=fn.radial_distance(com_star[id][0]-mbp[id][0],
#                                     com_star[id][1]-mbp[id][1],com_star[id][2]-mbp[id][2])/radius[id]#the halo doesn't have any subhalo

    BCGid[i]=id
    i=i+1
    continue
  
#  hosthalo is not included in the subhalo
  if star_lumz[id]>=np.max(star_lum_sub):
#    BCGoffset[i]=fn.radial_distance(com_star[id][0]-mbp[id][0],
#                                    com_star[id][1]-mbp[id][1],com_star[id][2]-mbp[id][2])/radius[id]

    BCGid[i]=id
  else:
    
    BCG_id=sub_id[star_lum_sub==np.max(star_lum_sub)][0]
    BCGid[i]=BCG_id
#    BCGoffset[i]=fn.radial_distance(com_star[BCG_id][0]-mbp[id][0],
#                                    com_star[BCG_id][1]-mbp[id][1],com_star[BCG_id][2]-mbp[id][2])/radius[id]
#    Rbri[i]=np.max(star_lum_sub)/star_lumz[id]
    BCGoffset[i]=fn.radial_distance(mbp[BCG_id][0]-mbp[id][0],
                                    mbp[BCG_id][1]-mbp[id][1],mbp[BCG_id][2]-mbp[id][2])/radius[id]
    if BCGoffset[i]>1:
      print(halo_id[id])
  i=i+1

'''
path="/home/jyang/data/Flamingo/L1000N1800/"
f=h5py.File(path+'BCG_exr100kpc.hdf5','w')
s=f.create_group("PartType0")
s.create_dataset("mainid", data=mainhalo_id)
s.create_dataset("BCGid", data=BCGid)
s.create_dataset("BCGoffset", data=BCGoffset)
s.create_dataset("Rbri", data=Rbri)
f.close()

'''
title="D of mbp BCG and central halo, r_star<100kpc"  
import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlabel("Offset/r_en")
ax.set_ylabel("Counts")
ax.set_title(title)
#for i in range(len(mainhalo_id)):
#    ax.plot(bin,rho_g[i],color='b')
ax.hist(BCGoffset,bins=50)
ax.set_yscale("log")
#ax.plot(np.logspace(0,4,100),np.logspace(0,4,100)*0.1,color='r')
fig.savefig("/home/jyang/plot/Flamingo/L1000N0900/offset_mbp_starcom_exr100kpc.png")
