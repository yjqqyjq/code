#visualize the idealized simulations
import numpy as np
import matplotlib.pyplot as plt
import h5py
import functions as fn
import matplotlib.pyplot as plt
from matplotlib import colors
from tqdm import tqdm
import os
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles_recenter"



center3=np.zeros(103)
center2=np.zeros(103)
centers=np.zeros(103)
i=0
for filename in tqdm(os.listdir(path)):
    
    f=h5py.File(path+"/"+filename,'r')
    Coord_s=np.array(f['PartType2']["Coordinates"])
    star_lum=np.array(f['PartType2']["lum_gamaz"])
    Coord_dm=np.array(f['PartType1']["Coordinates"])
    x_dm=np.array(Coord_dm).T[0]
  
 
    
  
    f.close()

   
    r=np.max(x_dm)
    Coord_s=Coord_s/r
    xyz_s=Coord_s.T
    ra=np.sqrt(xyz_s[0]**2+xyz_s[1]**2+xyz_s[2]**2)
    r_s=fn.radial_distance(xyz_s[0],xyz_s[1],xyz_s[2])
    x,y,z=fn.center_of_mass(xyz_s[0][ra<0.1],xyz_s[1][(ra<0.1)],xyz_s[2][ra<0.1],weight=star_lum[ra<0.1])
    centers[i]=fn.radial_distance(x,y,z)
    if centers[i]<0.005:
        print(filename,centers[i])
    i=i+1
#    Coord_s=Coord_s[ra<0.1]
#    star_lum=star_lum[ra<0.1*r]
#    print(Coord_s)
'''
    h3=np.histogramdd(Coord_s[ra<0.1],bins=100,range=[[-0.1,0.1],[-0.1,0.1],[-0.1,0.1]],weights=star_lum[ra<0.1])
    h2=np.histogram2d(xyz_s[0][ra<0.1],xyz_s[1][ra<0.1],bins=100,range=[[-0.1,0.1],[-0.1,0.1]],weights=star_lum[ra<0.1])
    max_position3d = np.unravel_index(np.argmax(h3[0]), h3[0].shape)
    n=100
    r3=(np.sqrt((max_position3d[0]-n/2)**2+(max_position3d[1]-n/2)**2)/n*2)/10
    center3[i]=r3
    max_position2d = np.unravel_index(np.argmax(h2[0]), h2[0].shape)
    r2=(np.sqrt((max_position2d[0]-n/2)**2+(max_position2d[1]-n/2)**2)/n*2)/10
    center2[i]=r2
    if (r3-r2)/r2>1:

        print(filename,r3,r2)


    


fig = plt.figure()
ax=plt.subplot(1,1,1) 
ax.hist(centers,bins=10)
#ax.scatter(center2,center3, s=0.5)
#ax.plot([0, 0.1], [0, 0.1], 'r--')
ax.set_yscale('log')
ax.set_xlabel("Offsets")
ax.set_ylabel("COunts")
#ax.set_title("LOS accumulated(2d) Brightest pixel vs 3d, stars r<0.1rvir, npixels=100")
ax.set_title("CoM_lum vs deepest potential, lum_z>10^6.5")

fig.savefig("/home/jyang/plot/Flamingo/L0200N0720/central_galaxy_offset_lumicut.png")
'''