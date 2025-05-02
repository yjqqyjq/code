import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw
import h5py
import functions as fn
from pathlib import Path
from tqdm import tqdm
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles/"
f=h5py.File(path+'66.hdf5','r')
Coord_g=np.array(f['PartType1']["Coordinates"])
x_g=Coord_g[:,0]
y_g=Coord_g[:,1]
z_g=Coord_g[:,2]

xray_lum=np.array(f['PartType1']['xray_lum_erosita_low'])+np.array(f['PartType1']['xray_lum_erosita_high'])
x_dm=np.array(f['PartType2']["Coordinates"]).T[0]
f.close()
ri=np.max(x_dm)
x_g/=ri
y_g/=ri
z_g/=ri
centeri=np.array([0.5,0.5,0.5])
x_g=x_g-centeri[0]
y_g=y_g-centeri[1]
z_g=z_g-centeri[2]
r_g=np.sqrt((x_g)**2+(y_g)**2+(z_g)**2)
rvir=1

for i in range(200):
    mask=(r_g<rvir)
    centeri=fn.center_of_mass(x_g[mask],y_g[mask],z_g[mask],xray_lum[mask])
    x_g=x_g-centeri[0]
    y_g=y_g-centeri[1]
    z_g=z_g-centeri[2]
    r_g=np.sqrt(x_g**2+y_g**2+z_g**2)
    rvir=rvir-1/201
    print(rvir,centeri)