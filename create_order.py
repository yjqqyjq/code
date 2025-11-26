#Create_Order
import numpy as np
import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool
def create_order(k):
    Tree_dir="/Users/24756376/data/Flamingo/L1000N0900/Trees/"
    path="/Users/24756376/data/Flamingo/L1000N0900/"
    folder_path ="/Users/24756376/data/Flamingo/L1000N0900/test/flamingo_0077/flamingo_0077."
    bound_path="/Users/24756376/data/Flamingo/L1000N0900/test/member/membership_0077."
    index=joblib.load(Tree_dir+str(k)+'.joblib')
    member_dm=(index[0][1].astype(np.int32))
    member_g=(index[1][1].astype(np.int32))
    index_dm=index[0][0].astype(np.int32)
    index_g=index[1][0].astype(np.int32)
    data=h5py.File(bound_path+str(k)+".hdf5",'r')
    bound_dm=np.array(data["PartType1"]["GroupNr_bound"],dtype=np.int32)[index_dm]
    bound_g=np.array(data["PartType1"]["GroupNr_bound"],dtype=np.int32)[index_g]
    data.close()
    return [member_dm,bound_dm,member_g,bound_g]
#rank index according to membership here
    
if __name__ == '__main__':
  tracemalloc.start()
  Tree_dir="/Users/24756376/data/Flamingo/L1000N0900/Trees/"

  folder_path ="/Users/24756376/data/Flamingo/L1000N0900/test/flamingo_0077/flamingo_0077."
  bound_path="/Users/24756376/data/Flamingo/L1000N0900/test//member/membership_0077."
  path="/Users/24756376/data/Flamingo/L1000N0900/"
  
  print("start") 


  m_dm=[]
  m_g=[]

  b_dm=[]
  b_g=[]

  with Pool(1) as p: # Create a pool with 5 worker processes


       results=p.map(create_order, range(0,2))



  for i in range(0,2):
      m_dm.append(results[i][0])

      m_g.append(results[i][2])
      b_dm.append(results[i][1])
      b_g.append(results[i][3])
  results=[]
  f=h5py.File(path+'halos_ranked.hdf5','r')
  id=np.array(f['id']).astype(np.float64)
  input_id=np.array(f['input_ids']).astype(np.int32)
  f.close()
  m_dm=np.concatenate(m_dm).astype(np.int32)
  m_g=np.concatenate(m_g).astype(np.int32)

  b_dm=np.concatenate(b_dm).astype(np.int32)
  b_g=np.concatenate(b_g).astype(np.int32)
  b_dm[np.isin(b_dm,input_id,invert=True)]=-1#call everything not belongs to halos unbound
  b_g[np.isin(b_g,input_id,invert=True)]=-1
  mask_dm=np.argsort(b_dm,kind="mergesort").astype(np.int32)
  mask_g=np.argsort(b_g,kind="mergesort").astype(np.int32)
  b_dm=b_dm[mask_dm]
  b_g=b_g[mask_g]
  print(b_dm[0])
  b_dm_c=b_dm[b_dm!=-1]
  b_g_c=b_g[b_g!=-1]
  m_dm_ub=m_dm[mask_dm][b_dm==-1]
  m_g_ub=m_g[mask_g][b_g==-1]
  
  
  
  

 
  order = np.empty(np.max(input_id)+1, dtype=np.int32)
  order[input_id] = np.arange(len(input_id), dtype=np.int32)
  
  mask_dm_c=np.argsort(order[b_dm_c],kind="mergesort").astype(np.int32)
  mask_g_c=np.argsort(order[b_g_c],kind="mergesort").astype(np.int32)
  mask_dm_ub=np.argsort(m_dm_ub,kind="mergesort").astype(np.int32)
  mask_g_ub=np.argsort(m_g_ub,kind="mergesort").astype(np.int32)
  N_dm = np.bincount(order[b_dm_c], minlength=np.max(m_dm)).astype(np.int32)
  N_g = np.bincount(order[b_g_c], minlength=np.max(m_dm)).astype(np.int32)
  N_dm_ub = np.bincount(m_dm_ub, minlength=np.max(m_dm)).astype(np.int32)
  N_g_ub = np.bincount(m_g_ub, minlength=np.max(m_dm)).astype(np.int32)

  N_dm_r = np.bincount(m_dm).astype(np.int32)   
  N_g_r = np.bincount(m_g).astype(np.int32)
  N_dm_c=np.zeros(len(id[id<=0]),dtype=np.int32)
  N_g_c=np.zeros(len(id[id<=0]),dtype=np.int32)

  b_dm=[]
  b_g=[]
  m_dm=[]
  m_g=[]
  for i in tqdm(range(0,len(id[id<=0]))):
    i_s=int(np.argwhere(id==-i))
    if len(np.argwhere(id==-i-1))==0:
        i_e=len(id)
    else:
        i_e=int(np.argwhere(id==-i-1))
    
    N_dm_c[i]=np.sum(N_dm[i_s:i_e])
    N_g_c[i]=np.sum(N_g[i_s:i_e])
     
  mask_dm=mask_dm[np.concatenate([mask_dm_ub,mask_dm_c+len(mask_dm_ub)])]
  mask_d=mask_g[np.concatenate([mask_g_ub,mask_g_c+len(mask_g_ub)])]
  tracemalloc.stop()
  '''
  f=h5py.File(path+'particle_orders.hdf5','w')
  f.create_dataset("mask_dm",data=o_dm)

  f.create_dataset("mask_g",data=o_g)


  f.close()
  


  fh=h5py.File(path+'halos_ranked.hdf5','a')


  fh.create_dataset("N_g_region",data=N_g)

  

  fh.create_dataset("N_dm_region",data=N_dm)

  



  

  fh.close()


  '''