from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw
import h5py
from scipy.spatial.transform import Rotation as Ro
from pathlib import Path
from tqdm import tqdm
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles/"
f=h5py.File(path+'66.hdf5','r')
Coord_g=np.array(f['PartType1']["Coordinates"])
xray_lum=np.array(f['PartType1']['xray_lum_erosita_low'])+np.array(f['PartType1']['xray_lum_erosita_high'])
x_dm=np.array(f['PartType2']["Coordinates"]).T[0]
f.close()
r=np.max(x_dm)
#histgram and smoothing
file=h5py.File("/home/jyang/data/Flamingo/L02000720/halo_2d_lum_rotate/2_kde1.hdf5",'r')
axis=np.array(file['PartType1']["rot_axis"]).T
file.close()
offset_3d=np.zeros(len(axis))
offset_2d=np.zeros(len(axis))
for i in tqdm(range(0,len(axis))):
  rotation = Ro.from_rotvec(-axis[i])
  rotated_points = rotation.apply(Coord_g)/r
  n=200  
#3d histogram
  hist=np.histogramdd(rotated_points,bins=n,range=[[-1,1],[-1,1],[-1,1]],weights=xray_lum)
  density=hist[0]
  edge=hist[1]
  density_smooth = gaussian_filter(density,sigma=10,mode='constant',truncate=3.0)
  max_position = np.unravel_index(np.argmax(density_smooth), density_smooth.shape)
  
  r3=np.sqrt((max_position[0]-n/2)**2+(max_position[1]-n/2)**2)/n*2
#2d histogram
  hist2d=np.histogram2d(rotated_points[:,0],rotated_points[:,1],bins=n,range=[[-1,1],[-1,1]],weights=xray_lum)
  density2d=hist2d[0]
  edge2d=hist2d[1]
  density2d_smooth = gaussian_filter(density2d,sigma=10,mode='constant',truncate=3.0)
  max_position_2d = np.unravel_index(np.argmax(density2d_smooth), density2d_smooth.shape)
  r2=np.sqrt((max_position_2d[0]-n/2)**2+(max_position_2d[1]-n/2)**2)/n*2
  offset_3d[i]=r3
  offset_2d[i]=r2
#  print(max_position,max_position_2d)
fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.set_title("3D and 2D offset,npixels=200,sigma=10")
ax.scatter(offset_3d,offset_2d,color='r',s=0.1)
ax.plot(np.arange(0,0.1,0.005),np.arange(0,0.1,0.005),color='k',linestyle='--')
ax.set_xlabel("3D offset")
ax.set_ylabel("2D offset")
fig.savefig("/home/jyang/plot/Flamingo/L0200N0720/halo66_2d_3d_offset_smoothed.png")