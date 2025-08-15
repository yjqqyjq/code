import numpy as np
import unyt
import swiftsimio as sw
import functions as fn
import h5py
from tqdm import tqdm
from PIL import Image#To plot image directly from pixle
#dir="/mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
#dir="../../../mnt/su3-pro/colibre/L0200N1504/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"




path="/Users/24756376/data/Flamingo/L1000N0900/"
#path="/home/jyang/data/Colibre/L0200N1504/"
f=h5py.File(path+'halos_ranked.hdf5','r')

halo_id=np.array(f["id"])

radius=(np.array(f["r200"]))[halo_id<=0]

mass=np.array(f["mass"])[halo_id<=0]

mbp=(np.array([f["centers_x"],f["centers_y"],f["centers_z"]]).T)[halo_id<=0]
com_star=np.array([f["com_star_100_x"],f["com_star_100_y"],f["com_star_100_z"]]).T

ms=np.array(f["ms_100"])
#msg=np.array(f["halos"]["mass_gas_bound"])
f.close()

#ignore all the halos without models of star lum and com

mainhalo_id=halo_id[halo_id<=0]

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

for id in tqdm(mainhalo_id):
  id=int(id)
  sub_id=halo_id[(halo_id<=-id+1)*(halo_id>=-id)]

  ms_sub=ms[(halo_id<=-id+1)*(halo_id>=-id)]#find BCG in mass
  com_sub=com_star[(halo_id<=-id+1)*(halo_id>=-id)]
 
  

  if len(ms_sub)==0:
#    BCGoffset[i]=fn.radial_distance(com_star[id][0]-mbp[id][0],
#                                     com_star[id][1]-mbp[id][1],com_star[id][2]-mbp[id][2])/radius[id]#the halo doesn't have any subhalo

    BCGid[i]=-1
    i=i+1
    continue
  
#  hosthalo is not included in the subhalo
  if ms[halo_id==id][0]>=np.max(ms_sub):
    d=com_star[halo_id==id][0]-mbp[-id]
   
    BCGoffset[i]=fn.radial_distance(d[0],d[1],d[2])/radius[-id]

    BCGid[i]=id
  else:
    
    BCG_id=sub_id[np.argmax(ms_sub)]
    BCGid[i]=BCG_id
    d=com_sub[np.argmax(ms_sub)]-mbp[-id]

    BCGoffset[i]=fn.radial_distance(d[0],d[1],d[2])/radius[-id]
#    Rbri[i]=np.max(star_lum_sub)/star_lumz[id]
    
   

  i=i+1

#mass=np.log10(mass)
bin_edge=10**np.linspace(4,5.5,21)
bins=np.digitize(mass,bin_edge)

fBCG=np.zeros(len(bin_edge)-1)
massbin=np.zeros(len(bin_edge)-1)
'''
for j in range(0,len(fBCG)):
  suboff=BCGoffset[bins==(j+1)]
  submain=mainhalo_id[bins==(j+1)]
  subbcg=BCGid[bins==(j+1)]
  massbin[j]=0.5*(bin_edge[j]+bin_edge[j+1])
  if len(suboff[(subbcg!=-1)])!=0:
    fBCG[j]=len(suboff[(submain!=subbcg)*(subbcg!=-1)])/len(suboff[(subbcg!=-1)])#the number of BCG miscenter/ number of clusters(exclude halos without satellites)
print()

path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'BCG_gasmass.hdf5','w')
#s=f["PartType0"]
#del s["fBCG_mscut"]
#s.create_dataset("fBCG_mscut", data=fBCG)


s=f.create_group("PartType0")
s.create_dataset("fBCG", data=fBCG)
s.create_dataset("massbin", data=massbin)
s.create_dataset("mainid", data=mainhalo_id)
s.create_dataset("BCGid", data=BCGid)
s.create_dataset("BCGoffset", data=BCGoffset)
s.create_dataset("Rbri", data=Rbri)
f.close()
'''
title=r"M>10^{14}M_\odot"  

import matplotlib.pyplot as plt
fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlabel("Distance/r200")
ax.set_ylabel("Frequency")
ax.set_title(r"$M>10^{14}M_\odot$"  )
#for i in range(len(mainhalo_id)):
#    ax.plot(bin,rho_g[i],color='b')
h1=np.histogram(BCGoffset[(BCGid<=0)*(BCGid!=-1)],bins=10**np.linspace(-4,1,51))
h2=np.histogram(BCGoffset[BCGid>0],bins=10**np.linspace(-4,1,51))

bin=10**np.linspace(-4,1,51)
ax.stairs(h1[0]/len(mainhalo_id),h1[1],color='b',label="BCG is center")
ax.stairs(h2[0]/len(mainhalo_id),h2[1],color='orange',label="BCG is not center")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()

#ax.plot(np.logspace(0,4,100),np.logspace(0,4,100)*0.1,color='r')
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/BCG_100kpc.png")
