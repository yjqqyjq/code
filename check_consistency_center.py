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
radius=np.array(data_h.radii.rvir)

is_mainhalo=data_h.centrals
halo_id=data_h.ids.id-1
mass=data_h.masses.mass_tot
xc=data_h.positions.xcmbp
yc=data_h.positions.ycmbp
zc=data_h.positions.zcmbp
#CoM in the catalogue
xcm=data_h.positions.xc-xc
ycm=data_h.positions.yc-xc
zcm=data_h.positions.zc-xc
xcg=data_h.positions.xc_gas+xcm#Center to mbp
ycg=data_h.positions.yc_gas+ycm
zcg=data_h.positions.zc_gas+zcm
xcs=data_h.positions.xc_star+xcm
ycs=data_h.positions.yc_star+ycm
zcs=data_h.positions.zc_star+zcm
data_h=[]
mask=(is_mainhalo)*(mass>100)*(mass<101)
mainhalo_id=halo_id[mask]
mainhalo_id=mainhalo_id.astype(int)
xcm=np.array(xcm[mask])
ycm=np.array(ycm[mask])
zcm=np.array(zcm[mask])
xcg=np.array(xcg[mask])
ycg=np.array(ycg[mask])
zcg=np.array(zcg[mask])
xcs=np.array(xcs[mask])
ycs=np.array(ycs[mask])
zcs=np.array(zcs[mask])
print(len(mainhalo_id))

R_center=np.zeros(len(mainhalo_id))
R_center_g=np.zeros(len(mainhalo_id))
R_center_s=np.zeros(len(mainhalo_id))

m_pdm=0.565
m_pg=0.1108
m_ps=0.075


for i in range (0,len(mainhalo_id)):
    
    id=mainhalo_id[i]
    
    centre=np.array([xc[id],yc[id],zc[id]])*unyt.Mpc
    r=radius[id]

    sgi=sg.SWIFTGalaxy(dir+"snapshots/flamingo_0008.hdf5",#"colibre_with_SOAP_membership_0127.hdf5",
                   sg.Standalone(centre=centre,velocity_centre=np.array([0,0,0])*
                                unyt.km/unyt.s,spatial_offsets=np.array([[-r,r],[-r,r],[-r,r]])*unyt.Mpc,extra_mask=None))
    mask=sg.MaskCollection(dark_matter=(sgi.dark_matter.spherical_coordinates.r<r),
                        gas=(sgi.gas.spherical_coordinates.r<r)*(sgi.gas.temperatures>0),
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

    xcg_comp,ycg_comp,zcg_comp=fn.center_of_mass(x_g,y_g,z_g)
    xcs_comp,ycs_comp,zcs_comp=fn.center_of_mass(x_stars,y_stars,z_stars)
    R_center_g[i]=fn.radial_distance(xcg[i]-xcg_comp,ycg[i]-ycg_comp,zcg[i]-zcg_comp)/fn.radial_distance(xcg[i],ycg[i],zcg[i])
    R_center_s[i]=fn.radial_distance(xcs[i]-xcs_comp,ycs[i]-ycs_comp,zcs[i]-zcs_comp)/fn.radial_distance(xcs[i],ycs[i],zcs[i])
    mhalo=(m_pdm*len(x_dm)+m_pg*len(x_g)+m_ps*len(x_stars))
    ycm_comp=(np.sum(m_pdm*y_dm)+np.sum(m_pg*y_g)+np.sum(m_ps*y_stars))/mhalo
    xcm_comp=(np.sum(m_pdm*x_dm)+np.sum(m_pg*x_g)+np.sum(m_ps*x_stars))/mhalo
    zcm_comp=(np.sum(m_pdm*z_dm)+np.sum(m_pg*z_g)+np.sum(m_ps*z_stars))/mhalo
    R_center[i]=fn.radial_distance(xcm[i]-xcm_comp,ycm[i]-ycm_comp,zcm[i]-zcm_comp)/fn.radial_distance(xcm[i],ycm[i],zcm[i])

boxused="/Flamingo/L0200N0360/"
import matplotlib.pyplot as plt
fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.plot(mainhalo_id,R_center,'r',label="CoM")
ax.plot(mainhalo_id,R_center_g,'g',label="CoM_g")
ax.plot(mainhalo_id,R_center_s,'b',label="CoM_s")
ax.set_xlabel("Halo")
ax.set_ylabel("Ratio")
ax.set_title("M=1to1.01e12")
ax.legend()
#fig.savefig("/home/jyang/plot/Flamingo/L0200N0360/Offset_small.png")
fig.savefig("/home/jyang/plot/"+boxused+"/Check_CoM_verysmall.png")
plt.close()