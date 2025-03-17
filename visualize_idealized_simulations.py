#visualize the idealized simulations
import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw
import h5py
import os
print(os.getcwd())
path="../../mnt/su3ctm/wmcdonald/merger_output_1_p2/0.00000_0.87034_15.000000/"
x_dm=[]
y_dm=[]
x_g=[]
y_g=[]
for i in range(0,16):
    if i<10:
      f=h5py.File(path+'snapshot_0'+str(i)+'0.hdf5','r')
    else:
      f=h5py.File(path+'snapshot_'+str(i)+'0.hdf5','r')
    Coord_g=f['PartType0']["Coordinates"]
    Coord_dm=f['PartType1']["Coordinates"]
    xyz_dm=np.array(Coord_dm).T 
    xyz_g=np.array(Coord_g).T
    f.close()
    x_dm.append(xyz_dm[0])
    y_dm.append(xyz_dm[1])
    x_g.append(xyz_g[0])
    y_g.append(xyz_g[1])
def radial_distance(x, y):
  return np.sqrt(x**2 + y**2)

def quadrupole(x, y, num_particles):
    r = radial_distance(x, y)
    
    return np.sqrt(15/(16*np.pi) * ( 4*np.sum(x*y/r)**2 + 
    np.sum(x**2/r - y**2/r)**2 )) / num_particles
   
def dissociation(x_dm, y_dm, x_g, y_g):
    n_dm = len(x_dm)
    n_g = len(x_g)
    r_mean_dm = np.average(radial_distance(x_dm, y_dm))
    r_mean_g = np.average(radial_distance(x_g, y_g))
    r_max = max(r_mean_dm, r_mean_g)
    q_dm = quadrupole(x_dm, y_dm, n_dm)
    q_g = quadrupole(x_g, y_g, n_g)
    return np.sqrt(4*np.pi/5) * (q_dm - q_g) / r_max

S=[]
for i in range(0,16):
  S.append(dissociation(x_dm[i], y_dm[i], x_g[i], y_g[i]))
import matplotlib.pyplot as plt
from matplotlib import colors
fig = plt.figure()
ax=plt.subplot(1,1,1) 
ax.plot(np.arange(0,16,1),S)
fig.savefig("/home/jyang/plot/idealizedsim/0.00000_0.87034_15.000000/S.png")

for i in range(0,16):
#  x=np.append(x_dm[i],x_g[i])
#  y=np.append(y_dm[i],y_g[i])
  
  fig = plt.figure()
  ax=plt.subplot(1,1,1)
  h=ax.hist2d(x_g[i],y_g[i],bins=1000,norm=colors.LogNorm())
  plt.colorbar(h[3], ax=ax)
  if i<10:
    fig.savefig('/home/jyang/plot/idealizedsim/0.00000_0.87034_15.000000/xyg0'+str(i)+'0.png')
  else:

    fig.savefig('/home/jyang/plot/idealizedsim/0.00000_0.87034_15.000000/xyg'+str(i)+'0.png')
  plt.close(fig)
