import numpy as np
from scipy.spatial.transform import Rotation as Ro
from PIL import Image
from pathlib import Path
import h5py
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles/"
f=h5py.File(path+'66.hdf5','r')
Coord_g=np.array(f['PartType1']["Coordinates"]).T
xray_lum=np.array(f['PartType1']['xray_lum_rosat'])


xg=Coord_g[0][xray_lum==np.max(xray_lum)][0]/3.048
yg=Coord_g[1][xray_lum==np.max(xray_lum)][0]/3.048
zg=Coord_g[2][xray_lum==np.max(xray_lum)][0]/3.048
print(xg,yg,zg)
xyz_maxlum=np.array([xg,yg,zg])

f.close()

file=h5py.File("/home/jyang/data/Flamingo/L02000720/halo_2d_lum_rotate/66_kde0.hdf5",'r')
axis=np.array(file['PartType1']["rot_axis"]).T
xlim=np.array(file['PartType1']['xlim'])
#rot=np.array(file['PartType1']['rot_matrix'])
file.close()
#print(axis)
#iterate over png file

rotation = Ro.from_rotvec(-axis)
rotated_points = rotation.apply(xyz_maxlum)

offsets=np.zeros(100)
offsets_3dproj=np.zeros(100)
i=0
folder = Path("/home/jyang/plot/Flamingo/L0200N0720/halo_2d_lum_rotate/66_kde0")
for png_file in folder.glob("*.png"):
   
#    rotated_points = np.matmul(rot[:,:,i].T,xyz_maxlum)
    img = Image.open(png_file)#.convert("L")
    img_array = np.array(img)
    pixel_x=len(img_array[0])
    pixel_y=len(img_array[:,0])
    max_position = np.unravel_index(np.argmax(img_array), img_array.shape)
    
    x=(max_position[1]/pixel_y-0.5)*2
    y=((1-max_position[0]/pixel_x)-0.5)*2
    r=np.sqrt(x**2+y**2)
    offsets[i]=r
    if r>0.5:
        print(png_file.name,axis[i])
    x3=rotated_points[i][0]
    y3=rotated_points[i][1]
    r3=np.sqrt(x3**2+y3**2)
    offsets_3dproj[i]=r3
    i=i+1
    
#    print(r/r3,png_file.name)
#        print(x,y)
import matplotlib.pyplot as plt
from matplotlib import colors
title="M>1e14"
boxused="/Flamingo/L0200N0720/"
'''
fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.hist(offsets,bins=10)
ax.set_xlabel("2doffsets/rvir")
ax.set_ylabel("Counts")
ax.set_title("Offsets of 2d xray lum max(from gagetry) and mbp")
fig.savefig("/home/jyang/plot/"+boxused+"/halo_2_2doffset_kde1.png")

fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(offsets,offsets_3dproj,s=1,color='r')
#ax.plot(np.linspace(0.1,0.2,100),np.linspace(0.1,0.2,100))

ax.set_xlabel("2doffsets/rvir")
ax.set_ylabel("3doffsets_proj/rvirs")
ax.set_title("Offsets of 2d and 3d xray lum max(from gagetry) and mbp")
fig.savefig("/home/jyang/plot/"+boxused+"/halo_66_2d_3d_offset_kde0.png")
'''