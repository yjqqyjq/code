import numpy as np
import unyt
import swiftsimio as sw
import functions as fn
import h5py
import tqdm
from PIL import Image#To plot image directly from pixle
#dir="/mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
#dir="../../../mnt/su3-pro/colibre/L0200N1504/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"




path="/Users/24756376/data/Flamingo/L1000N1800/"
#path="/home/jyang/data/Colibre/L0200N1504/"
f=h5py.File(path+'halos.hdf5','r')

halo_id=np.array(f["halos"]["id"])

radius=np.array(f["halos"]["r200"])
host_id=np.array(f["halos"]["hostid"])
mass=np.array(f["halos"]["mass"])
m200=np.array(f["halos"]["m200"])
mbp=np.array(f["halos"]["center"])
com_star=np.array(f["halos"]["com_star_100kpc"])
msg=np.array(f["halos"]["mass_gas_bound"])
ms=np.array(f["halos"]["mass_star_100kpc"])
f.close()

#ignore all the halos without models of star lum and com
rs=np.sqrt(com_star[:,0]**2+com_star[:,1]**2+com_star[:,2]**2)
mask=(rs>-1)
halo_id=halo_id[mask]
host_id=host_id[mask]
mbp=mbp[mask]
m200=m200[mask]
#com_star=com_star[mask]
#star_lumz=star_lumz[mask]
mass=mass[mask]
radius=radius[mask]
ms=ms[mask]
print(len(ms[(ms==0)*(host_id==-1)]))
star_lumz=msg
ids=np.arange(0,len(halo_id),1)
mainhalo_id=ids[(host_id==-1)*(mass>10000)]
R=np.ones(len(mainhalo_id))
mostmassid=np.zeros(len(mainhalo_id))
i=0

for id in (mainhalo_id):
  id=int(id)

  sub_id=ids[host_id==halo_id[id]]
  
  star_lum_sub=star_lumz[host_id==halo_id[id]]
  if len(sub_id)==0:
    R[i]=-1
    mostmassid[i]=-1
    i=i+1
  else:
    R[i]=np.max(star_lum_sub)/star_lumz[id]
    mostmassid[i]=int(sub_id[star_lum_sub==np.max(star_lum_sub)][0])
    i=i+1
print(halo_id[(host_id==-1)*(mass>10000)][R>1])

path="/Users/24756376/data/Flamingo/L1000N1800/"
f=h5py.File(path+'compare_most_massive_sat_gas.hdf5','w')


s=f.create_group("PartType0")
s.create_dataset("centralid", data=mainhalo_id)
s.create_dataset("Brightid", data=mostmassid)
s.create_dataset("Rmass", data=R)

f.close()
'''
Brightid=np.zeros(len(mainhalo_id))
Bright2nd_id=np.zeros(len(mainhalo_id))

i=0

for id in (mainhalo_id):
  id=int(id)

  sub_id=ids[host_id==halo_id[id]]
  
  star_lum_sub=star_lumz[host_id==halo_id[id]]
 
 
  

  if len(sub_id)==0:
    

    Brightid[i]=-1
    Bright2nd_id[i]=-1
    i=i+1
    continue
  
#  hosthalo is not included in the subhalo
  if star_lumz[id]>=np.max(star_lum_sub):
    

    Brightid[i]=id
    Bright2nd_id[i]=int(sub_id[star_lum_sub==np.max(star_lum_sub)][0])
#    if np.max(star_lum_sub)/star_lumz[int(Bright2nd_id[i])]!=1:
#      print(np.max(star_lum_sub),star_lumz[int(Bright2nd_id[i])])
      
  else:
    #purt all halos togethers to compare
#    subid=np.append(sub_id,id)
#    starlumsub=np.append(star_lum_sub,star_lumz[id])
    
    Brightid[i]=int(sub_id[star_lum_sub==np.max(star_lum_sub)][0])
    

    Bright2nd_id[i]=id
  
    
#   if Bright2nd_id[i]!=id:
#      print("fatal error, please smash your laptop immediately")
# 

   

  i=i+1


path="/Users/24756376/data/Flamingo/L1000N1800/"
f=h5py.File(path+'compare_brightnest_central_3000kpc.hdf5','w')


s=f.create_group("PartType0")
s.create_dataset("centralid", data=mainhalo_id)
s.create_dataset("Brightid", data=Brightid)
s.create_dataset("Bright2ndid", data=Bright2nd_id)

f.close()

title="D of mbp BCG and central halo, r_star<100kpc"  
import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlabel("Offset/r200")
ax.set_ylabel("Counts")
ax.set_title(title)
#for i in range(len(mainhalo_id)):
#    ax.plot(bin,rho_g[i],color='b')
ax.hist(BCGoffset,bins=50)
ax.set_yscale("log")
#ax.plot(np.logspace(0,4,100),np.logspace(0,4,100)*0.1,color='r')
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/BCG_offset_exr100kpc.png")
'''