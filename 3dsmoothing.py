from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
import unyt
import functions as fn
import h5py
from scipy.spatial.transform import Rotation as Ro
from pathlib import Path
from tqdm import tqdm
from matplotlib import colors
from PIL import Image
path="/Users/24756376/data/Flamingo/L1000N0900/"
id=0

main_id=fn.halo_ids[fn.halo_ids<=0]
mainarg=np.argwhere((fn.halo_ids==-id))
center=fn.centers[fn.halo_ids<=0][id]
#histgram and smoothing

particle=fn.load_particles(path,fn.halo_ids[mainarg[0]],dm=0,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":["xray_lum_erosita_low"],"stars":[]},mode="halo")
Coord_g=particle[0][0]-center
#Coord_g=Coord_g/np.max(Coord_g[:,0])
xlum=particle[0][1]
#offset_3d=np.zeros(len(axis))
#offset_2d=np.zeros(len(axis))
#rotate
'''
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
'''
hist=np.histogramdd(Coord_g,bins=200,range=[[-8,8],[-8,8],[-8,8]],weights=xlum)
#print(np.min(Coord_g[:,0]),np.max(Coord_g[:,0]))
density=hist[0]
edge=hist[1]
bins=10**np.linspace(np.log10(np.min(density[density!=0])),np.log10(np.max(density)),10)
binnum=np.digitize(density,bins)
img=np.sum(density,axis=2)
img_smooth=gaussian_filter(img,sigma=2.4375,mode='constant',truncate=4.0)

#sigma=0.13Mpc
#for i in range(0,len(bins)-1):
#  density_part=np.where(binnum==i+1,density,density)
#  density_smooth += gaussian_filter(density_part,sigma=1.625,mode='constant',truncate=3.0)
#max_position = np.unravel_index(np.argmax(density_smooth), density_smooth.shape)
print(np.sum(img),np.sum(img_smooth))
img=np.where(img<=0,0,np.log10(img))
img_smooth=np.where(img_smooth<=0,0,np.log10(img_smooth))

#img=np.log10(img)
print(np.sum(img),np.sum(img_smooth))
print(np.max(img),np.max(img_smooth))
print(np.min(img),np.min(img_smooth))
fig = plt.figure()

'''
ax.set_title("3D and 2D offset,npixels=200,sigma=10")
ax.scatter(offset_3d,offset_2d,color='r',s=0.1)
ax.plot(np.arange(0,0.1,0.005),np.arange(0,0.1,0.005),color='k',linestyle='--')
ax.set_xlabel("3D offset")
ax.set_ylabel("2D offset")
fig.savefig("/home/jyang/plot/Flamingo/L0200N0720/halo66_2d_3d_offset_smoothed.png")
'''
i=plt.imshow(img,cmap="Reds", extent=[-8,8,-8,8],vmax=10,vmin=0)

#i = Image.new("RGB", (200,200))
plt.xlabel("x(Mpc)")
plt.ylabel("y(Mpc)")

fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/smoothing.png")
#i.putdata(img)
#i.save("/Users/24756376/plot/Flamingo/L1000N0900/smoothing.png")
fig = plt.figure()

#print(np.min(img2[img2!=0]))

#img2=np.log10(img2)


i2=plt.imshow(img_smooth,cmap="Reds",extent=[-8,8,-8,8],vmax=10,vmin=0)
plt.xlabel("x(Mpc)")
plt.ylabel("y(Mpc)")

fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/smoothing_b.png")