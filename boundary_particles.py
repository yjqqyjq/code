#boundary particles
import swiftsimio as sw
import numpy as np
import h5py
import tracemalloc
from tqdm import tqdm
import sys
path="/Users/24756376/data/Flamingo/L1000N0900/"
#path="/Users/24756376/data/Colibre/L0012N0094/"
f=h5py.File(path+'halos_ranked.hdf5','r')
id=np.array(f['id'])
N_dm=np.array(f['N_dm'])
c=np.array(f['cross_bound'])
f.close()
print(id[(c==1)])
'''
main_id=id[id<0]
cross=np.zeros(len(N_dm))
f=h5py.File(path+'particles_ranked.hdf5','r')
Coord_dm=np.array(f["PartType1"]["Coordinates"],dtype=np.float32)
#print(Coord_dm)
f.close()
for i in tqdm(range(0,len(id))):
   
    sub_particles=Coord_dm[int(np.sum(N_dm[0:i])):int((np.sum(N_dm[0:i])+N_dm[i]))]
#    print(int(np.sum(N_dm[0:i])),int((np.sum(N_dm[0:i])+N_dm[i])))
    x_s=sub_particles[:,0]
    y_s=sub_particles[:,1]
    z_s=sub_particles[:,2]
    if np.max(x_s)-np.min(x_s)>500:
        cross[i]+=1
    if np.max(y_s)-np.min(y_s)>500:
        cross[i]+=0.1
    if np.max(z_s)-np.min(z_s)>500:
        cross[i]+=0.01
f=h5py.File(path+'halos_ranked.hdf5','a')
#del f["cross_bound"]
f.create_dataset("cross_bound",data=cross)
f.close()
'''