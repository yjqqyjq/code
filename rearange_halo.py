#rearange_halo
import numpy as np
import h5py
import tracemalloc
from tqdm import tqdm
import sys
#Qpath="/Users/24756376//data/Colibre/L0012N0094/"
path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'halos_cen13.hdf5','r')
table=[]
masskey=0
keys=[]
i=0
for key in f['halos'].keys():
    length=np.array(np.array(f['halos'][str(key)]).shape)
    
    if len(length)==2:
       keys.append(str(key)+'_x')
       keys.append(str(key)+'_y')
       keys.append(str(key)+'_z')
       table.append(np.array(f['halos'][str(key)])[:,0])
       table.append(np.array(f['halos'][str(key)])[:,1])
       table.append(np.array(f['halos'][str(key)])[:,2])
       i=i+3
    else:
      keys.append((str(key)))
      table.append(np.array(f['halos'][str(key)]))
      if key=="mass":
          masskey=i
      i=i+1
f.close()


#table=np.array(table).T
#print(table[0])
print(keys)
#rank by reverse mass first
table=np.array(table)
maskm = np.argsort(table[masskey],kind="mergesort")[::-1]
table=table[:,maskm]

id_comb=np.zeros(len(table[0]))
maincount=0

for i in tqdm(range(0,len(table[0]))):
    if table[keys.index("hostid")][i]==-1:
        id_comb[i]=-maincount
        mask=(table[keys.index("hostid")]==table[keys.index("id")][i])
        
        sub_id=table[keys.index("id")][mask]
        
        if len(sub_id>0):
            
            index=np.arange(1,len(sub_id)+1,1)
           
           
            id_comb[mask]=maincount+1/(len(sub_id)+2)*index
        maincount+=1
#print(table[3])
id_sort=np.where(id_comb<0,-id_comb,id_comb)
mask = np.argsort(id_sort,kind="mergesort")
table_new=table[:,mask]
#print(table_new[0])

#print(id_comb)
#print(id_comb)
f=h5py.File(path+'halos_cen13_ranked.hdf5','w')
f.create_dataset("id",data=id_comb[mask])
for i in range(0,len(table_new)):
  if keys[i]=="hostid" or keys[i]=="id":
      continue
  f.create_dataset(keys[i],data=table_new[i])
f.close()