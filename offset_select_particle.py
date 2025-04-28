import numpy as np
import unyt
import swiftsimio as sw
from velociraptor import load
import swiftgalaxy as sg
import functions as fn
from swiftgalaxy import SWIFTGalaxy, MaskCollection
from swiftsimio import cosmo_array
import tqdm
#load the data Flamingo

/mnt/su3ctm/ludlow/Flamingo/L1000N1800/HYDRO_FIDUCIAL_HiResDM/SOAP-HBT/
/mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5
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
    
Offsets=np.zeros(len(mainhalo_id))
CoMr=np.zeros(len(mainhalo_id))
#get the coordinates of the dark matter and gas particles
for i in tqdm.tqdm(range(0,len(mainhalo_id))):
    
    id=mainhalo_id[i]
    
    centre=np.array([xc[id],yc[id],zc[id]])*unyt.Mpc
    r=radius[id]
#"snapshots/flamingo_0008.hdf5",#"colibre_with_SOAP_membership_0127.hdf5",
    sgi=sg.SWIFTGalaxy(dir+"halo_properties_0077.hdf5",
                    sg.SOAP(dir+"halo_properties_0077.hdf5",soap_index=0,extra_mask=None))
#                   sg.Standalone(centre=centre,velocity_centre=np.array([0,0,0])*
#                                unyt.km/unyt.s,spatial_offsets=np.array([[-r,r],[-r,r],[-r,r]])*unyt.Mpc,extra_mask=None))
    mask=sg.MaskCollection(dark_matter=(sgi.dark_matter.spherical_coordinates.r<r),
                        gas=(sgi.gas.spherical_coordinates.r<0.1*r)*(sgi.gas.temperatures>0),
                        stars=(sgi.stars.spherical_coordinates.r<r))
    sgi.mask_particles(mask)
    

   
    x_dm=np.array(sgi.dark_matter.cartesian_coordinates.x)
    y_dm=np.array(sgi.dark_matter.cartesian_coordinates.y)
    z_dm=np.array(sgi.dark_matter.cartesian_coordinates.z)
    x_g=np.array(sgi.gas.cartesian_coordinates.x)
    y_g=np.array(sgi.gas.cartesian_coordinates.y)
    z_g=np.array(sgi.gas.cartesian_coordinates.z)
#    x_stars=np.array(sgi.stars.cartesian_coordinates.x)
#    y_stars=np.array(sgi.stars.cartesian_coordinates.y)
#    z_stars=np.array(sgi.stars.cartesian_coordinates.z)
#    xlum=np.array(sgi.gas.xray_luminosities.erosita_high+sgi.gas.xray_luminosities.erosita_low+sgi.gas.xray_luminosities.ROSAT)
#    print(len(xlum[xlum==0])/len(xlum))
#    xlc=x_g[xlum==np.max(xlum)]
#    ylc=y_g[xlum==np.max(xlum)]/float(r)
#    zlc=z_g[xlum==np.max(xlum)]
#    Rho_g=np.array(sgi.gas.densities)
    
##    h=np.histogram2d(x_g,y_g,bins=1000,weights=xlum)
#    print(h[1])
#    lum=np.reshape(h[0],shape=1000*1000)
#    index=np.argmax(lum)
#    xlc=h[1][index//1000]+0.0005
#    ylc=h[2][index%1000]+0.0005
#    print(xlc,ylc)
    
    Offsets[i]=fn.offset(x_dm,y_dm,z_dm,x_g,y_g,z_g)/float(r)/10
#    print(Offsets[i])
#    Offsets[i]=fn.offset_lum(x_g,y_g,z_g,Rho_g)/float(r)
#    print(Offsets[i])

  



'''
  
#calculate the offset and make plots
S=np.zeros(len(mainhalo_id))

Offsets_b=np.zeros(len(mainhalo_id))
Center_diff=np.zeros(len(mainhalo_id))
#for i in range(0,len(mainhalo_id)):
    
#    id=mainhalo_id[i]
   

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
np.savetxt("/home/jyang/data/"+boxused+title+"_offset_in.txt",np.array([Offsets]), header="Select all the gas particles inside 0.1 rvir,calculate the offset between these gas and all the DM normorized by 0.1rvir"
            , comments='#' )

#np.savetxt("/home/jyang/data/Flamingo/L0200N0360/M=5_to_6e12.txt",np.array([S,Offsets,Center_diff]))




fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Offsets, bins=20)
ax.set_xlabel("Offsets")
ax.set_ylabel("Counts")
ax.set_title(title)
#fig.savefig("/home/jyang/plot/Flamingo/L0200N0360/Offset_small.png")
fig.savefig("/home/jyang/plot/"+boxused+"/Offset_in.png")
plt.close()

'''

fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(S,Offsets,s=1)
ax.set_xlabel("S")
ax.set_ylabel("Offset")
ax.set_title(title)
fig.savefig("/home/jyang/plot/"+boxused+"/Co_S_Offset.png")
plt.close()
'''