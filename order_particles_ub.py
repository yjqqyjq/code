import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool
def create_order(k):
    print(k)
    path="/Users/24756376/data/Flamingo/L1000N1800/"
    Tree_dir="/Users/24756376/data/Flamingo/L1000N1800/Trees/"
    index=joblib.load(path+"Trees/"+str(k)+'.joblib')


    member_dm=index[0][1].astype(np.int32)
    member_g=index[1][1].astype(np.int32)
    index=[]
  
    return [member_dm,member_g]
#rank index according to membership here

if __name__ == '__main__':
  tracemalloc.start()
  path="/Users/24756376/data/Flamingo/L1000N1800/"
  Tree_dir="/Users/24756376/data/Flamingo/L1000N1800/Trees/"

  print("start") 

  
  m_dm=[]
  m_g=[]
  

  with Pool(3) as p: # Create a pool with 5 worker processes


       results=p.map(create_order, range(0,64))



  for i in range(0,64):
      m_dm.append(results[i][0])

      m_g.append(results[i][1])
     
  results=[]
  
  m_dm=np.concatenate(m_dm).astype(np.int32)
  m_g=np.concatenate(m_g).astype(np.int32)

 
  #rank unbound
  mask_dm_ub=np.argsort(m_dm,kind="mergesort").astype(np.int32)
  mask_g_ub=np.argsort(m_g,kind="mergesort").astype(np.int32)


  
 
  

 
 
#  N_dm = np.bincount(m_dm, minlength=np.max(m_dm)+1).astype(np.int32)
#  N_g = np.bincount(m_g, minlength=np.max(m_g)+1).astype(np.int32)


  

  tracemalloc.stop()
#  f=h5py.File(path+'/cluster1/mask.hdf5','r')
  
  
#  mask=np.array(f["mask_g_ub"])

##  print(mask,mask_g_ub)
#  f.close()

  f=h5py.File(path+'cluster14/mask.hdf5','w')
  
  f.create_dataset("mask_dm_ub",data=mask_dm_ub)

  f.create_dataset("mask_g_ub",data=mask_g_ub)


  f.close()
   
  
  '''
  fh=h5py.File(path+'halos_13_ranked.hdf5','a')
  del fh["N_g_region"]
  del fh["N_dm_region"]

  fh.create_dataset("N_g_region",data=N_g)

  

  fh.create_dataset("N_dm_region",data=N_dm)

  



  

  fh.close()

  '''