#load membership
import unyt
import swiftsimio as sw
import numpy as np
import h5py
import tracemalloc
import sys
from tqdm import tqdm
tracemalloc.start()
path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'halos.hdf5','r')
halo_id=np.array(f["halos"]["id"])
host_id=np.array(f["halos"]["hostid"])
input_id=np.array(f["halos"]["input_ids"])
mass=np.array(f["halos"]["mass"])
print(len(host_id))
f.close()
ids=np.arange(0,len(halo_id),1) 
mrank=np.argsort(mass)
#print(mrank)
print(tracemalloc.get_traced_memory())
f=h5py.File(path+'cluster_particles.hdf5','r')
dm_group=np.array(f["PartType1"]["member"],dtype=np.int32)#align with input id
dm_group_new=np.array(f["PartType1"]["member"],dtype=np.float32)
f.close()


print(tracemalloc.get_traced_memory())
tracemalloc.stop()
#np.where(host_id==-1,dm_group_new,)
for i in tqdm(range(0,len(input_id))):
    if host_id[mrank[i]]==-1:
        dm_group_new[dm_group==input_id[mrank[i]]]=i
    else:
        n_p=len(dm_group_new[dm_group==input_id[mrank[i]]])
        dm_group_new[dm_group==input_id[mrank[i]]]=ids[halo_id==host_id[i]]+n_p/10**len(str(n_p))
f=h5py.File(path+'cluster_particles.hdf5','a')
dm=f["PartType1"]
dm.create_dataset("member_new",dm_group_new)
f.close()