#rearange_particle_region
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
path="/Users/24756376/data/Flamingo/L1000N1800/"
#path="/home/jyang/Colibre/L0012N0094/"

tracemalloc.start()
keys=[]
#read the membership
print("#read the membership")
f=h5py.File(path+'particles_radius.hdf5','r')
member_dm=np.array(f['PartType1']['members'],dtype=np.int32)
member_g=np.array(f['PartType0']['members'],dtype=np.int32)
member_s=np.array(f['PartType2']['members'],dtype=np.int32)
f.close()
  

print(tracemalloc.get_traced_memory())




#path="/Users/24756376/data/Colibre/L0012N0094/"
print("#count")

# np.array(list(chain(*b)))

mask_dm = np.argsort(member_dm,kind="mergesort").astype(int)
mask_g = np.argsort(member_g,kind="mergesort").astype(int)
mask_s = np.argsort(member_s,kind="mergesort").astype(int)

N_dm = np.bincount(member_dm[mask_dm]).astype(np.int32)   
N_g = np.bincount(member_g[mask_g]).astype(np.int32)  
N_s = np.bincount(member_s[mask_s]).astype(np.int32)  

print("finish sorting")
print(member_dm[mask_dm])
print(len(member_dm))
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
vkey_dm=0
vkey_g=0
vkey_s=0
#read the table and sort by the new membership sequence
print("save the data")
f=h5py.File(path+'particles_radius.hdf5','r')
PartType0=[]
PartType1=[]
PartType2=[]

for key in f['PartType1'].keys():
      if key=='members':
         continue   
      if key=='Coordinates':
         keys_dm.append('Coordinates')
         keys_dm.append('Coordinates')
         keys_dm.append('Coordinates')
         Coord=np.array(f['PartType1']['Coordinates'])
         PartType1.append(Coord[:,0])
         PartType1.append(Coord[:,1])
         PartType1.append(Coord[:,2])
         Coord=[]
         cokey_dm=i_dm
         i_dm=i_dm+3
      elif key=='Velocities':
         keys_dm.append('Velocities')
         keys_dm.append('Velocities')
         keys_dm.append('Velocities')
         Velo=np.array(f['PartType1']['Velocities'])
         PartType1.append(Velo[:,0])
         PartType1.append(Velo[:,1])
         PartType1.append(Velo[:,2])
         Velo=[]
         vkey_dm=i_dm
         i_dm=i_dm+3
      else:
        keys_dm.append((str(key)))
        PartType1.append(np.array(f['PartType1'][str(key)],dtype=np.float32))
        i_dm=i_dm+1
dm_new=np.array(PartType1)[:,mask_dm]
for key in f['PartType0'].keys():
      
      if key=='members':
         continue   
      if key=='Coordinates':
         keys_g.append('Coordinates')
         keys_g.append('Coordinates')
         keys_g.append('Coordinates')
         Coord=np.array(f['PartType0']['Coordinates'])
         PartType0.append(Coord[:,0])
         PartType0.append(Coord[:,1])
         PartType0.append(Coord[:,2])
         Coord=[]
         cokey_g=i_g
         i_g=i_g+3
      elif key=='Velocities':
         keys_g.append('Velocities')
         keys_g.append('Velocities')
         keys_g.append('Velocities')
         Velo=np.array(f['PartType0']['Velocities'])
         PartType0.append(Velo[:,0])
         PartType0.append(Velo[:,1])
         PartType0.append(Velo[:,2])
         Velo=[]
         vkey_g=i_g
         i_g=i_g+3
         

      else:
        
        keys_g.append((str(key)))
        PartType0.append(np.array(f['PartType0'][str(key)],dtype=np.float32))
        i_g=i_g+1

gas_new=np.array(PartType0)[:,mask_g]
for key in f['PartType2'].keys():
      if key=='members':
         continue   
      if key=='Coordinates':
         keys_s.append('Coordinates')
         keys_s.append('Coordinates')
         keys_s.append('Coordinates')
         Coord=np.array(f['PartType2']['Coordinates'])
         PartType2.append(Coord[:,0])
         PartType2.append(Coord[:,1])
         PartType2.append(Coord[:,2])
         Coord=[]
         cokey_s=i_s
         i_s=i_s+3
      elif key=='Velocities':
         keys_s.append('Velocities')
         keys_s.append('Velocities')
         keys_s.append('Velocities')
         Velo=np.array(f['PartType2']['Velocities'])
         PartType2.append(Velo[:,0])
         PartType2.append(Velo[:,1])
         PartType2.append(Velo[:,2])
         Velo=[]
         vkey_s=i_s
         i_s=i_s+3
      else:
        keys_s.append((str(key)))
        PartType2.append(np.array(f['PartType2'][str(key)],dtype=np.float32))
        i_s=i_s+1
star_new=np.array(PartType2)[:,mask_s]
f.close()



#print(table_new)
#print(N_g)




#print(id_comb)
#print(id_comb)
#save the data to new file

f=h5py.File(path+'halos_ranked.hdf5','a')
del f["N_g_region"]
del f["N_dm_region"]
del f["N_s_region"]
f.create_dataset("N_g_region",data=N_g)
f.create_dataset("N_dm_region",data=N_dm)
f.create_dataset("N_s_region",data=N_s)
f.close()
'''
f=h5py.File(path+'particles_region_ranked.hdf5','w')
dm=f.create_group("PartType1")
for j in range(0,len(dm_new)):
    if keys_dm[j]=="Coordinates" or keys_dm[j]=="Velocities":
        continue
    else:
      dm.create_dataset(keys_dm[j],data=dm_new[j])
dm.create_dataset("Coordinates",data=np.array([dm_new[cokey_dm],dm_new[cokey_dm+1],dm_new[cokey_dm+2]]).T)
dm.create_dataset("Velocities",data=np.array([dm_new[vkey_dm],dm_new[vkey_dm+1],dm_new[vkey_dm+2]]).T)
g=f.create_group("PartType0")
for j in range(0,len(gas_new)):
    print(j,keys_g[j],len(gas_new))
    if keys_g[j]=="Coordinates" or keys_g[j]=="Velocities":
        continue
    else:
      print(keys_g[j],gas_new[j])
      g.create_dataset(keys_g[j],data=gas_new[j])
g.create_dataset("Coordinates",data=np.array([gas_new[cokey_g],gas_new[cokey_g+1],gas_new[cokey_g+2]]).T)
g.create_dataset("Velocities",data=np.array([gas_new[vkey_g],gas_new[vkey_g+1],gas_new[vkey_g+2]]).T)
s=f.create_group("PartType2")
for j in range(0,len(star_new)):
    if keys_s[j]=="Coordinates" or keys_s[j]=="Velocities":
        continue
    else:
      s.create_dataset(keys_s[j],data=star_new[j])
s.create_dataset("Coordinates",data=np.array([star_new[cokey_s],star_new[cokey_s+1],star_new[cokey_s+2]]).T)
s.create_dataset("Velocities",data=np.array([star_new[vkey_s],star_new[vkey_s+1],star_new[vkey_s+2]]).T)
f.close()
'''