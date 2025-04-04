#visualize the idealized simulations
import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw
import h5py
import functions as fn
import matplotlib.pyplot as plt
from matplotlib import colors
para="1.00000_0.32404_15.000000/"
path="../../mnt/su3ctm/wmcdonald/merger_output_1_p2/"+para+"/"
x_dm=[]
y_dm=[]
z_dm=[]
x_g=[]
y_g=[]
z_g=[]
n_files=100
'''
f0=h5py.File(path+'snapshot_000.hdf5','r')
Coord_g0=f0['PartType0']["Coordinates"]
Coord_dm0=f0['PartType1']["Coordinates"]
xyz_dm0=np.array(Coord_dm0).T
xyz_g0=np.array(Coord_g0).T
f0.close()
x_dm0=xyz_dm0[0]
x_g0=xyz_g0[0]
mask_dm=np.ones(len(x_dm0))#(x_dm0<(np.max(x_dm0)+np.min(x_dm0))/2)
mask_g=np.ones(len(x_g0))#(x_g0<(np.max(x_g0)+np.min(x_g0))/2)
'''
for i in range(0,n_files):
    if i<10:
      f=h5py.File(path+'snapshot_00'+str(i)+'.hdf5','r')
    elif i<100:
      f=h5py.File(path+'snapshot_0'+str(i)+'.hdf5','r')
    else:
      f=h5py.File(path+'snapshot_'+str(i)+'.hdf5','r')
    Coord_g=f['PartType0']["Coordinates"]
    Coord_dm=f['PartType1']["Coordinates"]
    xyz_dm=np.array(Coord_dm).T 
    xyz_g=np.array(Coord_g).T
    f.close()
    x_dm.append(xyz_dm[0])
    y_dm.append(xyz_dm[1])
    z_dm.append(xyz_dm[2])
    x_g.append(xyz_g[0])
    y_g.append(xyz_g[1])
    z_g.append(xyz_g[2])
print(len(x_g))

S=np.zeros(n_files)
Offset=np.zeros(n_files)
for i in range(0,n_files):
#  S[i]=fn.dissociation(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
  Offset[i]=fn.offset(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])


fig = plt.figure()
ax=plt.subplot(1,1,1) 
ax.plot(np.arange(0,n_files,1)*0.1,Offset)
ax.set_xlabel("Time(Gyr)")
ax.set_ylabel("Offset")
fig.savefig("/home/jyang/plot/idealizedsim/"+para+"/Offset.png")
'''
for i in range(0,10):
#  x=np.append(x_dm[i],x_g[i])
#  y=np.append(y_dm[i],y_g[i])
  
  fig = plt.figure()
  ax=plt.subplot(1,1,1)
  h=ax.hist2d(x_dm[i*10],y_dm[i*10],bins=1000,norm=colors.LogNorm())
  plt.colorbar(h[3], ax=ax)
  if i*10<10:
    fig.savefig('/home/jyang/plot/idealizedsim/'+para+'/xydm0'+str(i)+'01h.png')
  else:

    fig.savefig('/home/jyang/plot/idealizedsim/'+para+'/xydm'+str(i)+'01h.png')
  plt.close(fig)

for i in range(0,10):
#  x=np.append(x_dm[i],x_g[i])
#  y=np.append(y_dm[i],y_g[i])
  
  fig = plt.figure()
  ax=plt.subplot(1,1,1)
  h=ax.hist2d(x_g[i*10],y_g[i*10],bins=1000,norm=colors.LogNorm())
  plt.colorbar(h[3], ax=ax)
  if i*10<10:
    fig.savefig('/home/jyang/plot/idealizedsim/'+para+'/xyg0'+str(i)+'01h.png')
  else:

    fig.savefig('/home/jyang/plot/idealizedsim/'+para+'/xyg'+str(i)+'01h.png')
  plt.close(fig)'
  ='''