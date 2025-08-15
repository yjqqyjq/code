#rearange_particles
import numpy as np
import h5py
import tracemalloc
from tqdm import tqdm
import sys
import itertools
import multiprocessing as mp 
from functools import partial
import datetime
from collections import Counter

#     member_new_dm.append(index)
#     return len(index)



print(datetime.datetime.now())
path="/Users/24756376/data/Flamingo/L1000N0900/"
#path="/home/jyang/Colibre/L0012N0094/"
f=h5py.File(path+'halos_ranked.hdf5','r')
id=np.array(f['id']).astype(np.float64)
input_id=np.array(f['input_ids']).astype(np.int64)
f.close()

tracemalloc.start()
keys=[]
#read the membership
print("#read the membership")
f=h5py.File(path+'cluster_particles.hdf5','r')
member_dm=np.array(f['PartType1']['member'],dtype=np.int32)
member_g=np.array(f['PartType0']['member'],dtype=np.int32)
member_s=np.array(f['PartType2']['member'],dtype=np.int32)
f.close()
  

print(tracemalloc.get_traced_memory())




#path="/Users/24756376/data/Colibre/L0012N0094/"
print("#count")

# np.array(list(chain(*b)))
order = np.empty(np.max(member_dm)+1, dtype=np.int32)
order[input_id] = np.arange(len(input_id), dtype=np.int32)
mask_dm=np.argsort(order[member_dm]).astype(np.int32)

member_dm_new=order[member_dm].astype(np.int32)
N_dm = np.bincount(member_dm_new, minlength=len(input_id)).astype(np.int32)
print(tracemalloc.get_traced_memory())
member_dm_new=[]

mask_g=np.argsort(order[member_g]).astype(np.int32)
member_g_new=order[member_g].astype(np.int32)
N_g= np.bincount(member_g_new, minlength=len(input_id)).astype(np.int32)
member_g_new=[]

mask_s=np.argsort(order[member_s]).astype(np.int32)
member_s_new=order[member_s].astype(np.int32)
N_s= np.bincount(member_s_new, minlength=len(input_id)).astype(np.int32)
member_s_new=[]
order=[]


#clusters
N_dm_c=np.zeros(len(id[id<=0]),dtype=np.int32)
N_g_c=np.zeros(len(id[id<=0]),dtype=np.int32)
N_s_c=np.zeros(len(id[id<=0]),dtype=np.int32)
for i in tqdm(range(0,len(id[id<=0]))):
    i_s=int(np.argwhere(id==-i))
    if len(np.argwhere(id==-i-1))==0:
        i_e=len(id)
    else:
        i_e=int(np.argwhere(id==-i-1))
    
    N_dm_c[i]=np.sum(N_dm[i_s:i_e])
    N_g_c[i]=np.sum(N_g[i_s:i_e])
    N_s_c[i]=np.sum(N_s[i_s:i_e])

print("finish sorting")

member_dm=[]
member_g=[]
member_s=[] 




print(tracemalloc.get_traced_memory())
tracemalloc.stop()
print(datetime.datetime.now())

keys_dm=[]
keys_g=[]
keys_s=[]
i_dm=0
i_g=0
i_s=0
cokey_dm=0
cokey_g=0
cokey_s=0
#read the table and sort by the new membership sequence
print("save the data")
f=h5py.File(path+'cluster_particles.hdf5','r')
PartType0=[]
PartType1=[]
PartType2=[]

for key in f['PartType1'].keys():
      if key=='member':
         continue   
      if key=='Coordinates':
         keys_dm.append('Coordinates')
         keys_dm.append('Coordinates')
         keys_dm.append('Coordinates')
         Coord=np.array(f['PartType1']['Coordinates'],dtype=np.float32)
         PartType1.append(Coord[:,0])
         PartType1.append(Coord[:,1])
         PartType1.append(Coord[:,2])
         Coord=[]
         cokey=i_dm
      else:
        keys_dm.append((str(key)))
        PartType1.append(np.array(f['PartType1'][str(key)],dtype=np.float32))
        i_dm=i_dm+1
PartType1=np.array(PartType1,dtype=np.float32)
for key in f['PartType0'].keys():
      if key=='member':
         continue   
      if key=='Coordinates':
         keys_g.append('Coordinates')
         keys_g.append('Coordinates')
         keys_g.append('Coordinates')
         Coord=np.array(f['PartType0']['Coordinates'],dtype=np.float32)
         PartType0.append(Coord[:,0])
         PartType0.append(Coord[:,1])
         PartType0.append(Coord[:,2])
         Coord=[]
         cokey=i_g
      else:
        keys_g.append((str(key)))
        PartType0.append(np.array(f['PartType0'][str(key)],dtype=np.float32))
        i_g=i_g+1
PartType0=np.array(PartType0,dtype=np.float32)
for key in f['PartType2'].keys():
      if key=='member':
         continue   
      if key=='Coordinates':
         keys_s.append('Coordinates')
         keys_s.append('Coordinates')
         keys_s.append('Coordinates')
         Coord=np.array(f['PartType2']['Coordinates'],dtype=np.float32)
         PartType2.append(Coord[:,0])
         PartType2.append(Coord[:,1])
         PartType2.append(Coord[:,2])
         Coord=[]
         cokey=i_s
      else:
        keys_s.append((str(key)))
        PartType2.append(np.array(f['PartType2'][str(key)],dtype=np.float32))
        i_s=i_s+1
PartType2=np.array(PartType2,dtype=np.float32)  
f.close()

dm_new=PartType1[:,mask_dm]
PartType1=[]
gas_new=PartType0[:,mask_g]
PartType0=[]
star_new=PartType2[:,mask_s]
PartType2=[]


#save the data to new file

f=h5py.File(path+'halos_ranked.hdf5','a')
f.create_dataset("N_g",data=N_g)
f.create_dataset("N_dm",data=N_dm)
f.create_dataset("N_s",data=N_s)
f.create_dataset("N_g_c",data=N_g_c)
f.create_dataset("N_dm_c",data=N_dm_c)
f.create_dataset("N_s_c",data=N_s_c)
f.close()
'''
f=h5py.File(path+'particles_ranked_13.hdf5','w')
dm=f.create_group("PartType1")
for j in range(0,len(dm_new)):
    if keys_dm[j]=="Coordinates":
        continue
    dm.create_dataset(keys_dm[j],data=dm_new[j])
dm.create_dataset("Coordinates",data=np.array([dm_new[cokey_dm],dm_new[cokey_dm+1],dm_new[cokey_dm+2]],dtype=np.float32).T)
g=f.create_group("PartType0")
for j in range(0,len(gas_new)):
    if keys_g[j]=="Coordinates":
        continue
    g.create_dataset(keys_g[j],data=gas_new[j])
g.create_dataset("Coordinates",data=np.array([gas_new[cokey_g],gas_new[cokey_g+1],gas_new[cokey_g+2]],dtype=np.float32).T)
s=f.create_group("PartType2")
for j in range(0,len(star_new)):
    if keys_s[j]=="Coordinates":
        continue
    s.create_dataset(keys_s[j],data=star_new[j])
s.create_dataset("Coordinates",data=np.array([star_new[cokey_s],star_new[cokey_s+1],star_new[cokey_s+2]],dtype=np.float32).T)
f.close()

 '''
 
