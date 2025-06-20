#rearange_halo
import numpy as np
import h5py
import tracemalloc
from tqdm import tqdm
import sys
#Qpath="/Users/24756376//data/Colibre/L0012N0094/"
path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'halos.hdf5','r')
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

table=np.array(table)

#print(table[3])
mask = np.argsort(table[masskey],kind="mergesort")[::-1]
table_new=table[:,mask]
#print(table_new[0])
id_comb=np.zeros(len(table_new[0]))
maincount=0

for i in tqdm(range(0,len(table_new[0]))):
    if table_new[keys.index("hostid")][i]==-1:
        id_comb[i]=-maincount
        mask=(table_new[keys.index("hostid")]==table_new[keys.index("id")][i])
        
        sub_id=table_new[keys.index("id")][mask]
        
        if len(sub_id>0):
            
            index=np.arange(1,len(sub_id)+1,1)
           
           
            id_comb[mask]=maincount+1/(len(sub_id)+2)*index
        maincount+=1
#print(id_comb)
#print(id_comb)
f=h5py.File(path+'halos_ranked.hdf5','w')
f.create_dataset("id",data=id_comb)
for i in range(0,len(table_new)):
  if keys[i]=="hostid" or keys[i]=="id":
      continue
  f.create_dataset(keys[i],data=table_new[i])
f.close()