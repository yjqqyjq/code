#find neighbor, deal with space
import faiss
import h5py
import glob
from tqdm import tqdm
import numpy as np
import datetime
import tracemalloc
dir="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(dir+'halos_ranked.hdf5','r')
tracemalloc.start()
centers=np.array([f["center_x"],f["center_y"],f["center_z"]],dtype=np.float32).T
r200=np.array(f["r200"])
ids=np.array(f['id'])
f.close()
r200=r200[ids<=0]

centers=centers[ids<=0]

files_dm=glob.glob(dir+'faiss_index/*dm.faiss')
files_g=glob.glob(dir+'faiss_index/*g.faiss')
#files_s=glob.glob(dir+'faiss_index/*s.faiss')
i=0

for file in files_dm:
    if i==0:
      
        index_dm=faiss.read_index(file)
 
        
    else:
        index_dm.merge_from(faiss.read_index(file))
        
      
    i+=1
i=0
for file in files_g:
    if i==0:
        index_g=faiss.read_index(file)
    else:
        index_g.merge_from(faiss.read_index(file))
        
    i+=1
centers=np.array_split(centers, 20)
print(centers[0])
#quary
#print(len(centers),len(index_dm))
print(tracemalloc.get_traced_memory())
print("start")
print(datetime.datetime.now())
dm=[]
g=[]
faiss.omp_set_num_threads(10)
for center in centers:
  dm_c=np.array(index_dm.range_search(center, 2**2)[0])
  g_c=np.array(index_g.range_search(center, 2**2)[0])
  dm.append(dm_c[:1]-dm_c[:-1])
  g.append(g_c[:1]-g_c[:-1])
  print(tracemalloc.get_traced_memory())
print(datetime.datetime.now())
dm=np.concatenate(dm)
g=np.concatenate(g)
tracemalloc.stop()
f=h5py.File(dir+'neighbor/2Mpc.hdf5','w')
f.create_dataset('dm', data=dm)
f.create_dataset('g', data=g)
f.close()