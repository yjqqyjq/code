#Xray Lum
import numpy as np
import h5py
from tqdm import tqdm
path="/Users/24756376/data/Flamingo/L1000N0900/"
f=h5py.File(path+'halos_ranked.hdf5','r')
N_g=np.array(f['N_g'])
id=np.array(f['id'])
f.close()
f=h5py.File(path+'particles_ranked.hdf5','r')
Xray=np.array(f['PartType0']['xray_lum_erosita_low'],dtype=np.float32)
f.close()
main_id=id[id<=0]

Rlum=np.zeros(len(main_id))
Rmass=np.zeros(len(main_id))
for i in tqdm(range(0,len(main_id))):
   subarg=np.argwhere((id>i)*(id<i+1)).astype(int)
   subhalos=[]
   if len(subarg)==0:
      Rlum[i]=-1
      continue
   for arg in subarg:
      arg=int(arg)
      if arg!=0:
        id_s=int(np.sum(N_g[0:arg]))
        id_e=int(np.sum(N_g[0:arg+1]))
        lum=np.sum(Xray[id_s:id_e])
        subhalos.append(lum)
   if subhalos!=[]:
     
     lum_sub=np.max(subhalos)
     mainarg=int(np.argwhere(id==-i))
   
     lum_main=np.sum(Xray[int(np.sum(N_g[0:mainarg])):int(np.sum(N_g[0:mainarg+1]))])
    
     Rlum[i]=lum_sub/lum_main
print(np.histogram(Rlum[Rlum>0],bins=20))