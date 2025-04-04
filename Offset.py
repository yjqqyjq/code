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
xcm=data_h.positions.xc
ycm=data_h.positions.yc
zcm=data_h.positions.zc-xc
xcg=data_h.positions.xc_gas#Center to CoM
ycg=data_h.positions.yc_gas
zcg=data_h.positions.zc_gas
xcs=data_h.positions.xc_star
ycs=data_h.positions.yc_star
zcs=data_h.positions.zc_star
data_h=[]
mask=(is_mainhalo)*(mass>100)
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