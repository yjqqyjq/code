import numpy as np
import h5py
import tracemalloc
from tqdm import tqdm
import sys
import itertools
import multiprocessing as mp 
from functools import partial
import datetime
#path="/Users/24756376/data/Colibre/L0012N0094/"
path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'halos_ranked.hdf5','r')
id=np.array(f['id']).astype(np.float64)

input_id=np.array(f['input_ids']).astype(np.int32)
f.close()

f=h5py.File(path+'cluster_particles.hdf5','r')
member_dm=np.array(f['PartType1']['member'],dtype=np.int32)
member_g=np.array(f['PartType0']['member'],dtype=np.int32)
X=np.array(f['PartType0']['xray_lum_erosita_low'],dtype=np.float32).T
f.close()

order = np.empty(np.max(member_dm)+1, dtype=np.int32)
#print(max(input_id))
order[input_id] = np.arange(len(input_id), dtype=np.int32)
mask=np.argsort(order[member_g])

sorted_B = X[mask]
member_new=member_g[mask]
print(member_new[-1])
print(np.argwhere(input_id==member_new[-1]))