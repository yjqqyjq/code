#LSS
import matplotlib.pyplot as plt
import numpy as np
import h5py
path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'halos_ranked.hdf5','r')
centers_14=np.array([f["centers_x"],f["centers_y"],f["centers_z"]]).T
f.close()
f=h5py.File(path+'halos_13_ranked.hdf5','r')
centers_13=np.array([f["centers_x"],f["centers_y"],f["centers_z"]]).T
f.close()
centers_14=centers_14[(centers_14[:,2]>995)*(centers_14[:,0]>000)*(centers_14[:,1]>000)]
centers_13=centers_13[(centers_13[:,2]>995)*(centers_13[:,0]>000)*(centers_13[:,1]>000)]
f=h5py.File('/Users/24756376/data/Flamingo/L1000N0900/cluster_particles0.8333333333333333_unbound.hdf5','r')
Coord_dm=np.array(f["PartType1"]["Coordinates"], dtype=np.float32)
Coord_dm=Coord_dm[(Coord_dm[:,2]>995)*(Coord_dm[:,0]>000)*(Coord_dm[:,1]>000)]
print(np.min(Coord_dm[:,2]),np.max(Coord_dm[:,2]))
f.close()
fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(centers_14[:,0],centers_14[:,1],s=1,alpha=1,label="14")
ax.scatter(centers_13[:,0],centers_13[:,1],s=1,alpha=1,label="13")
ax.scatter(Coord_dm[:,0],Coord_dm[:,1],s=0.001,alpha=0.1,color="grey")
ax.legend()
ax.set_xlabel("X")
ax.set_ylabel("Y")

fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/halo_positions.png")