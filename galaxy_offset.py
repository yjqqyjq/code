#visualize the idealized simulations
import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw
import h5py
import functions as fn
import matplotlib.pyplot as plt
from matplotlib import colors
import os
path="/home/jyang/data/Flamingo/L0200N0360/halo_particles_recenter"



centers=np.zeros(103)
i=0
for filename in os.listdir(path):
    f=h5py.File(path+"/"+filename,'r')
    Coord_s=np.array(f['PartType2']["Coordinates"])
    star_lum=np.array(f['PartType2']["lum_gamaz"])
    Coord_dm=np.array(f['PartType1']["Coordinates"])
    x_dm=np.array(Coord_dm).T[0]
  
 
    
  
    f.close()

   
    r=np.max(x_dm)
    xyz_s=Coord_s.T/r
    r_s=fn.radial_distance(xyz_s[0],xyz_s[1],np.zeros(len(xyz_s[0])))
    x,y,z=fn.center_of_mass(xyz_s[0][r_s<0.1],xyz_s[1][r_s<0.1],np.zeros(len(xyz_s[0][r_s<0.1])),weight=star_lum[r_s<0.1])
    centers[i]=fn.radial_distance(x,y,z)
    if centers[i]>0.01:
        print(filename,centers[i])
    i=i+1
    


fig = plt.figure()
ax=plt.subplot(1,1,1) 
ax.hist(centers,bins=10)
ax.set_xlabel("Offsets")
ax.set_ylabel("Counts")
ax.set_title("CoM_lum, stars x^2+y^2<0.1rvir, vs the deepest potential")
fig.savefig("/home/jyang/plot/Flamingo/L0200N0720/central_galaxy_offset_proj.png")
