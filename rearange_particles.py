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

def sort_dm(dm,ids):
     return np.argwhere(dm==ids).astype(np.int32).T[0]
#     member_new_dm.append(index)
#     return len(index)
def sort_g(g,ids):
     return np.argwhere(g==ids).astype(np.int32).T[0]
def sort_s(s,ids):
     return np.argwhere(s==ids).astype(np.int32).T[0]

if __name__ == '__main__':
  path="/Users/24756376/data/Flamingo/L1000N0900/"
#  path="/Users/24756376/data/Colibre/L0012N0094/"
  f=h5py.File(path+'halos_ranked.hdf5','r')
  id=np.array(f['id']).astype(np.float64)
  input_id=np.array(f['input_ids']).astype(np.int64)
  f.close()


  keys=[]
#read the membership
  f=h5py.File(path+'cluster_particles.hdf5','r')
  member_dm=np.array(f['PartType1']['member'],dtype=np.int32)
  member_g=np.array(f['PartType0']['member'],dtype=np.int32)
  member_s=np.array(f['PartType2']['member'],dtype=np.int32)
  f.close()
  mp.freeze_support()
  tracemalloc.start()




#from tqdm_multiprocess.logger import setup_logger_tqdm
#



  member_new_dm=np.zeros(len(member_dm))
  member_new_g=np.zeros(len(member_g))
  member_new_s=np.zeros(len(member_g))

  N_g=np.zeros(len(input_id))
  N_dm=np.zeros(len(input_id))
  N_s=np.zeros(len(input_id))

#calculate the index of h enew membership, defined by the sequence of input_id
#def sort_part(j):

 
   

# np.array(list(chain(*b)))


 
  with mp.Pool(1) as p:
     fdm = partial(sort_dm, member_dm)
     fg=partial(sort_g, member_g)
     fs=partial(sort_s, member_s)
#     print(partial_function)
     print(tracemalloc.get_traced_memory())
     member_new_dm=p.map(fdm,input_id)
     member_new_g=p.map(fg,input_id)
     member_new_s=p.map(fs,input_id)
     
  tracemalloc.stop()
  for k in tqdm(range(0,len(input_id))):
      N_dm[k]=len(member_new_dm[k])
      N_g[k]=len(member_new_g[k])
      N_s[k]=len(member_new_s[k])
  member_new_dm=np.array(list(itertools.chain(*member_new_dm)))   
  member_new_g=np.array(list(itertools.chain(*member_new_g)))  
  member_new_s=np.array(list(itertools.chain(*member_new_s)))  
#for k in tqdm(range(0,len(N_dm))):
#   sort_part(k)

#member_new_dm=member_new_dm.astype(np.int32)
##member_new_g=member_new_g.astype(np.int32)
#member_new_s=member_new_s.astype(np.int32)
  member_dm=[]
  member_g=[]
  member_s=[]
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

  dm_new=PartType1[:,member_new_dm]
  PartType1=[]
  gas_new=PartType0[:,member_new_g]
  PartType0=[]
  star_new=PartType2[:,member_new_s]
  PartType2=[]


#print(table_new)
#print(N_g)




#print(id_comb)
#print(id_comb)
#save the data to new file

  f=h5py.File(path+'halos_ranked.hdf5','a')
  del f['N_g']
  del f['N_dm']
  del f['N_s']
  f.create_dataset("N_g",data=N_g)
  f.create_dataset("N_dm",data=N_dm)
  f.create_dataset("N_s",data=N_s)
  f.close()

  f=h5py.File(path+'particles_ranked.hdf5','w')
  dm=f.create_group("PartType1")
  for j in range(0,len(dm_new)):
    if keys_dm[j]=="Coordinates":
        continue
    dm.create_dataset(keys_dm[j],data=dm_new[j])
  dm.create_dataset("Coordinates",data=np.array([dm_new[cokey_dm],dm_new[cokey_dm+1],dm_new[cokey_dm+2]],dtype=np.float32))
  g=f.create_group("PartType0")
  for j in range(0,len(gas_new)):
    if keys_g[j]=="Coordinates":
        continue
    g.create_dataset(keys_g[j],data=gas_new[j])
  g.create_dataset("Coordinates",data=np.array([gas_new[cokey_g],gas_new[cokey_g+1],gas_new[cokey_g+2]],dtype=np.float32))
  s=f.create_group("PartType2")
  for j in range(0,len(star_new)):
    if keys_s[j]=="Coordinates":
        continue
    s.create_dataset(keys_s[j],data=star_new[j])
  s.create_dataset("Coordinates",data=np.array([star_new[cokey_s],star_new[cokey_s+1],star_new[cokey_s+2]],dtype=np.float32))
  f.close()


'''
for j in tqdm(range(0,len(N_dm))):
   index_dm=np.argwhere(member_dm==input_id[j]).astype(np.int32).T[0]
   
   index_g=np.argwhere(member_g==input_id[j]).astype(np.int32).T[0]
 
   index_s=np.argwhere(member_s==input_id[j]).astype(np.int32).T[0]
  

   if len(index_dm)!=0:

     member_new_dm[int(np.sum(N_dm)):int(np.sum(N_dm)+len(index_dm))]=index_dm
     N_dm[j]=len(index_dm)
   if len(index_g)!=0:

     member_new_g[int(np.sum(N_g)):int(np.sum(N_g)+len(index_g))]=index_g
     N_g[j]=len(index_g)
   if len(index_s)!=0:
   
     member_new_s[int(np.sum(N_s)):int(np.sum(N_s)+len(index_s))]=index_s
     N_s[j]=len(index_s)

   if j%1000==0:
      print(len(N_dm[N_dm!=0]))
'''





 
 
