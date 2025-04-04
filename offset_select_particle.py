import numpy as np
import unyt
import swiftsimio as sw
from velociraptor import load
import swiftgalaxy as sg
import functions as fn
from swiftgalaxy import SWIFTGalaxy, MaskCollection
from swiftsimio import cosmo_array
#load the data Flamingo


dir="/mnt/su3-pro/flamingo/L0200N0360/"

data_h=load(dir+"VR/halos_0008.properties.0")
radius=data_h.radii.rvir

is_mainhalo=data_h.centrals
halo_id=data_h.ids.id-1
mass=data_h.masses.mass_tot
xc=data_h.positions.xcmbp
yc=data_h.positions.ycmbp
zc=data_h.positions.zcmbp

data_h=[]
mainhalo_id=halo_id[is_mainhalo*(mass>10000)]
mainhalo_id=mainhalo_id.astype(int)

print((len(mainhalo_id)))


m_pdm=0.565
m_pg=0.1108
m_ps=0.075

'''

#load the data Colibre
dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/"

data_h=sw.load(dir+"/halo_properties_0127.hdf5")
host_id=data_h.soap.host_halo_index#central halo=-1\
halo_id=np.arange(0,len(host_id),1)
mass=data_h.spherical_overdensity_200_crit.total_mass
radius=data_h.bound_subhalo.enclose_radius
xc=data_h.bound_subhalo.centre_of_mass[:,0]
yc=data_h.bound_subhalo.centre_of_mass[:,1]
zc=data_h.bound_subhalo.centre_of_mass[:,2]

mainhalo_id=halo_id[(host_id==-1)*(mass>1)]
m_pdm=1.937#10^7
m_pg=1.493#10^7
m_ps=1.124#10^7

x_dm=[[]]
y_dm=[[]]
z_dm=[[]]
x_g=[[]]
y_g=[[]]
z_g=[[]]
x_stars=[[]]
y_stars=[[]]
z_stars=[[]]


for i in range(0, len(mainhalo_id)-1):
  
    x_dm.append([])
    x_g.append([])
    z_dm.append([])
    y_dm.append([])
    y_g.append([])
    z_g.append([])
    x_stars.append([])
    y_stars.append([])
    z_stars.append([])
   
'''   
    
Offsets_out=np.zeros(len(mainhalo_id))

#get the coordinates of the dark matter and gas particles
for i in range (0,len(mainhalo_id)):
    
    id=mainhalo_id[i]
    
    centre=np.array([xc[id],yc[id],zc[id]])*unyt.Mpc
    r=radius[id]

    sgi=sg.SWIFTGalaxy(dir+"snapshots/flamingo_0008.hdf5",#"colibre_with_SOAP_membership_0127.hdf5",
                   sg.Standalone(centre=centre,velocity_centre=np.array([0,0,0])*
                                unyt.km/unyt.s,spatial_offsets=np.array([[-r,r],[-r,r],[-r,r]])*unyt.Mpc,extra_mask=None))
    mask=sg.MaskCollection(dark_matter=(sgi.dark_matter.spherical_coordinates.r<0.1*r),
                        gas=(sgi.gas.spherical_coordinates.r<0.1*r)*(sgi.gas.temperatures>0),
                        stars=(sgi.stars.spherical_coordinates.r<r))
    sgi.mask_particles(mask)
    

   
    x_dm=np.array(sgi.dark_matter.cartesian_coordinates.x)
    y_dm=np.array(sgi.dark_matter.cartesian_coordinates.y)
    z_dm=np.array(sgi.dark_matter.cartesian_coordinates.z)
    x_g=np.array(sgi.gas.cartesian_coordinates.x)
    y_g=np.array(sgi.gas.cartesian_coordinates.y)
    z_g=np.array(sgi.gas.cartesian_coordinates.z)
    x_stars=np.array(sgi.stars.cartesian_coordinates.x)
    y_stars=np.array(sgi.stars.cartesian_coordinates.y)
    z_stars=np.array(sgi.stars.cartesian_coordinates.z)
    Offsets_out[i]=fn.offset(x_dm,y_dm,z_dm,x_g,y_g,z_g)/float(r)

  




  
#calculate the offset and make plots
S=np.zeros(len(mainhalo_id))

Offsets_b=np.zeros(len(mainhalo_id))
Center_diff=np.zeros(len(mainhalo_id))
#for i in range(0,len(mainhalo_id)):
    
#    id=mainhalo_id[i]
#    S[i]=fn.dissociation(x_dm[i],y_dm[i],z_dm[i],x_g[i],y_g[i],z_g[i])
#    Offsets[i]=fn.offset(x_dm[i],y_dm[i],z_dm[i],x_g[i],y_g[i],z_g[i])/float(radius[id])
#    Offsets[i]=fn.offsetb(x_dm[i],y_dm[i],z_dm[i],x_g[i],y_g[i],z_g[i],x_stars[i],y_stars[i],z_stars[i],m_pg,m_ps)/float(radius[id])   
#    Offsets_s[i]=fn.offset(x_dm[i],y_dm[i],z_dm[i],x_stars[i],y_stars[i],z_stars[i])/float(radius[id])
'''    

    x_cdm,y_cdm,z_cdm=fn.center_of_mass(x_dm[i],y_dm[i],z_dm[i])
    
    x_cg,y_cg,z_cg=fn.center_of_mass(x_g[i],y_g[i],z_g[i])
    x_cs,y_cs,z_cs=fn.center_of_mass(x_stars[i],y_stars[i],z_stars[i])
    
    x_c=(x_cdm*m_pdm*len(x_dm[i])+x_cg*m_pg*len(x_g[i])+x_cs*m_ps*len(x_stars[i]))/(m_pdm*len(x_dm[i])+m_pg*len(x_g[i])+m_ps*len(x_stars[i]))
    y_c=(y_cdm*m_pdm*len(y_dm[i])+y_cg*m_pg*len(y_g[i])+y_cs*m_ps*len(y_stars[i]))/(m_pdm*len(y_dm[i])+m_pg*len(y_g[i])+m_ps*len(y_stars[i]))
    z_c=(z_cdm*m_pdm*len(z_dm[i])+z_cg*m_pg*len(z_g[i])+z_cs*m_ps*len(z_stars[i]))/(m_pdm*len(z_dm[i])+m_pg*len(z_g[i])+m_ps*len(z_stars[i]))
    
    Center_diff[i]=fn.radial_distance(x_c,y_c,z_c)/float(radius[id])
'''

import matplotlib.pyplot as plt
title="M>1e14"
boxused="/Flamingo/L0200N0360/"
np.savetxt("/home/jyang/data/"+boxused+title+"_offset_in.txt",np.array([Offsets_out]),comments
="Select all the particles inside 0.1rvir")

#np.savetxt("/home/jyang/data/Flamingo/L0200N0360/M=5_to_6e12.txt",np.array([S,Offsets,Center_diff]))
'''
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(S, bins=10)
ax.set_xlabel("S")
ax.set_ylabel("Counts")
ax.set_title(title)
#fig.savefig("/home/jyang/plot/Flamingo/L0200N0360/S_small.png")
fig.savefig("/home/jyang/plot/"+boxused+"/S_hot.png")
plt.close()

'''


fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Offsets_out, bins=10)
ax.set_xlabel("Offsets_in")
ax.set_ylabel("Counts")
ax.set_title(title)
#fig.savefig("/home/jyang/plot/Flamingo/L0200N0360/Offset_small.png")
fig.savefig("/home/jyang/plot/"+boxused+"/Offset_in.png")
plt.close()
'''
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Center_diff, bins=10)
ax.set_xlabel("Center_diff")
ax.set_ylabel("Counts")
ax.set_title(title)
fig.savefig("/home/jyang/plot/"+boxused+"/Centerdiff.png")
plt.close()

fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(S,Offsets,s=1)
ax.set_xlabel("S")
ax.set_ylabel("Offset")
ax.set_title(title)
fig.savefig("/home/jyang/plot/"+boxused+"/Co_S_Offset.png")
plt.close()
'''