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
N_dm_c=np.array(f['N_dm_c'])
N_g_c=np.array(f['N_g_c'])
N_s_c=np.array(f['N_s_c'])
#c=np.array(f['cross_bound'])
f.close()
#print(id[(c==1)])

main_id=id[id<0]

f=h5py.File(path+'particles_ranked.hdf5','r')
Coord_dm=np.array(f["PartType1"]["Coordinates"])
#Coord_g=np.array(f["PartType0"]["Coordinates"])
#Coord_s=np.array(f["PartType2"]["Coordinates"])
#print(Coord_dm)
f.close()
for i in tqdm(range(0,100)):
   
    sub_particles=Coord_dm[int(np.sum(N_dm_c[0:i])):int((np.sum(N_dm_c[0:i])+N_dm_c[i]))]
#    dm_particles=Coord_dm[int(np.sum(N_dm[0:i])):int((np.sum(N_dm[0:i])+N_dm[i]))]
#    dm_particles=Coord_dm[int(np.sum(N_dm[0:i])):int((np.sum(N_dm[0:i])+N_dm[i]))]

    x_s=sub_particles[:,0]
    y_s=sub_particles[:,1]
    z_s=sub_particles[:,2]
    if np.max(x_s)-np.min(x_s)>500:
        print(i,np.max(x_s),np.min(x_s))
