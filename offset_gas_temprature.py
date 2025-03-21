import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio import load
import swiftgalaxy as sg
import functions as fn
#load the data
soap_dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/"

data_h=sw.load(soap_dir+"/halo_properties_0127.hdf5")
xc=data_h.bound_subhalo.centre_of_mass[:,0]
yc=data_h.bound_subhalo.centre_of_mass[:,1]
zc=data_h.bound_subhalo.centre_of_mass[:,2]
radius=data_h.bound_subhalo.enclose_radius
boxsize=data_h.metadata.boxsize[0]


host_id=data_h.soap.host_halo_index#central halo=-1\
halo_id=np.arange(0,len(host_id),1)
mass=data_h.spherical_overdensity_200_crit.total_mass
data_h=[]

mainhalo_id=halo_id[(host_id==-1)*(mass>1)]
x_dm=[[]]
y_dm=[[]]
z_dm=[[]]
x_g=[[]]
y_g=[[]]
z_g=[[]]
T=[[]]
for i in range(0, len(mainhalo_id)-1):
  
    x_dm.append([])
    x_g.append([])
    z_dm.append([])
    y_dm.append([])
    y_g.append([])
    z_g.append([])
    T.append([])


#get the coordinates of the dark matter and gas particles
for i in range (0,len(mainhalo_id)):
    
    id=mainhalo_id[i]
   
    mask=sw.mask(soap_dir+"colibre_with_SOAP_membership_0127.hdf5",spatial_only=False)
    mask.constrain_spatial([np.array([xc[id]-radius[id],xc[id]+radius[id]])*unyt.Mpc,
                           np.array([yc[id]-radius[id],yc[id]+radius[id]])*unyt.Mpc,
                           np.array([zc[id]-radius[id],zc[id]+radius[id]])*unyt.Mpc])
#    mask.constrain_mask("gas","temperatures",35000,1e6)
    data=sw.load(soap_dir+"colibre_with_SOAP_membership_0127.hdf5",mask=mask)
    x_dm[i]=data.dark_matter.coordinates[:,0]
    y_dm[i]=data.dark_matter.coordinates[:,1]
    z_dm[i]=data.dark_matter.coordinates[:,2]
    x_g[i]=data.stars.coordinates[:,0]
    y_g[i]=data.stars.coordinates[:,1]
    z_g[i]=data.stars.coordinates[:,2]
#    T[i]=data.gas.temperatures
#some particles are outside the box
    if xc[id]-radius[id]<0:
        mask=x_dm[i]>boxsize/2
        x_dm[i][mask]=x_dm[i][mask]-boxsize
        mask=x_g[i]>boxsize/2
        x_g[i][mask]=x_g[i][mask]-boxsize
    elif xc[id]+radius[id]>boxsize:
        mask=x_dm[i]<boxsize/2
        x_dm[i][mask]=x_dm[i][mask]+boxsize
        mask=x_g[i]<boxsize/2
        x_g[i][mask]=x_g[i][mask]+boxsize
    if yc[id]-radius[id]<0:
        mask=y_dm[i]>boxsize/2
        y_dm[i][mask]=y_dm[i][mask]-boxsize
        mask=y_g[i]>boxsize/2
        y_g[i][mask]=y_g[i][mask]-boxsize
    elif yc[id]+radius[id]>boxsize:
        mask=y_dm[i]<boxsize/2
        y_dm[i][mask]=y_dm[i][mask]+boxsize
        mask=y_g[i]<boxsize/2
        y_g[i][mask]=y_g[i][mask]+boxsize
    if zc[id]-radius[id]<0:
        mask=z_dm[i]>boxsize/2
        z_dm[i][mask]=z_dm[i][mask]-boxsize
        mask=z_g[i]>boxsize/2
        z_g[i][mask]=z_g[i][mask]-boxsize
    elif zc[id]+radius[id]>boxsize:
        mask=z_dm[i]<boxsize/2
        z_dm[i][mask]=z_dm[i][mask]+boxsize
        mask=z_g[i]<boxsize/2
        z_g[i][mask]=z_g[i][mask]+boxsize
#calculate the offset and make plots
Offsets=np.zeros(len(mainhalo_id))
for i in range(0,len(mainhalo_id)):
#    mask=(T[i]<35000)
    xg=x_g[i]#[mask]
    yg=y_g[i]#[mask]
    zg=z_g[i]#[mask]
    Offsets[i]=fn.offset(x_dm[i],y_dm[i],z_dm[i],xg,yg,zg)

import matplotlib.pyplot as plt


fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Offsets*1000, bins=100)
ax.set_xlabel("Offsets(kpc)")
ax.set_ylabel("Counts")
ax.set_title("M>10^10")
fig.savefig("/home/jyang/plot/Colibre/L0012N0094/Offset_stars.png")
plt.close()
'''
Offsets=np.zeros(len(mainhalo_id))
for i in range(0,len(mainhalo_id)):
    mask=(T[i]>35000)
    xg=x_g[i][mask]
    yg=y_g[i][mask]
    zg=z_g[i][mask]
    Offsets[i]=fn.offset(x_dm[i],y_dm[i],z_dm[i],xg,yg,zg)


fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Offsets*1000, bins=100)
ax.set_xlabel("Offsets(kpc)")
ax.set_ylabel("Counts")
ax.set_title("M>10^10")
fig.savefig("/home/jyang/plot/Colibre/L0012N0094/Offset_hot.png")
plt.close()
'''