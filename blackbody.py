import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw
import h5py
#xray_lum offset>0.2:2.,  13.,  21.,  50.,  64.,  66.,  82., 105
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles/"


f=h5py.File(path+'50.hdf5','r')

Coord_g=f['PartType1']["Coordinates"]
Coord_dm=f['PartType1']["Coordinates"]
xyz_dm=np.array(Coord_dm).T 
xyz_g=np.array(Coord_g).T
xray_lum=np.array(f['PartType1']['xray_lum_rosat'])
#xray_lum=np.array(f['PartType1']['xray_lum_erosita_low'])+np.array(f['PartType1']['xray_lum_erosita_high'])
T=np.array(f['PartType1']['temperatures'])
f.close()
x_dm=xyz_dm[0]
y_dm=xyz_dm[1]
x_g=xyz_g[0]
y_g=xyz_g[1]
z_g=xyz_g[2]
r_g=np.sqrt(x_g**2+y_g**2+z_g**2)/np.max(x_g)
import matplotlib.pyplot as plt
from matplotlib import colors
title="M>1e14"
boxused="/Flamingo/L0200N0720/"
#  x=np.append(x_dm[i],x_g[i])
#  y=np.append(y_dm[i],y_g[i])
fig = plt.figure()
ax=plt.subplot(1,1,1)
B=xray_lum
ax.scatter(r_g[xray_lum==0],xray_lum[xray_lum==0],color='r',s=0.5)
fig.savefig("/home/jyang/plot/"+boxused+"/blackbody50.png")
