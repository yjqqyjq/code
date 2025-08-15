#Kdtree
from sklearn.neighbors import KDTree
import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool

dir="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/snapshots/flamingo_0077/flamingo_0077."
data_dir="/home/jyang/data/Flamingo/L1000N0900/Trees/"
def build_tree(k):
   print("starting", k)
   file=dir+str(k)+'.hdf5'
   f = h5py.File(file, 'r')
   dm = np.array(f['PartType1']['Coordinates'])
   Tree_dm=KDTree(dm)
   dm=[]
   print(tracemalloc.get_traced_memory())
   g=np.array(f['PartType0']['Coordinates'])
   Tree_g=KDTree(g)
   g=[]
   f.close()
   print(tracemalloc.get_traced_memory())
   joblib.dump(Tree_dm, data_dir+str(k)+'dm.pkl')
   joblib.dump(Tree_g, data_dir+str(k)+'g.pkl')
   print("Done")


def query_tree(k):
   print("starting", k)
   Tree_dm=joblib.load(data_dir+str(k)+'dm.pkl')
   
   index,D=Tree_dm.query_radius(
      centers,
      r=1,
      return_distance=True,
      sort_results=True,
      
   )
   joblib.dump([index,D], str(k)+'dm.joblib')

   
   Tree_dm=[]

   print(tracemalloc.get_traced_memory())
   Tree_g=joblib.load(data_dir+str(k)+'g.pkl')
   index,D=Tree_g.query_radius(
        centers,
        r=1,
        return_distance=True,
        sort_results=True,
        
    )
   joblib.dump([index,D], str(k)+'g.joblib')
   Tree_g=[]
   
   print(tracemalloc.get_traced_memory())
if __name__ == '__main__':
    tracemalloc.start()
    f=h5py.File("/home/jyang/data/Flamingo/L1000N0900/halos_ranked.hdf5", 'r')
    centers=np.array([f['center_x'],f['center_y'],f['center_z']]).T
    ids=np.array(f['id'])
    centers = centers[ids <= 0]
    print(len(centers))
    f.close()
'''
    with Pool() as p:
        p.map(query_tree, range(0,32))
    tracemalloc.stop()
'''
