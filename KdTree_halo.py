#KdTree_halo
from scipy.spatial import KDTree
import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool

dir="/Users/24756376/data/Flamingo/L1000N0900/halos_cen13_ranked.hdf5"
data_dir="/Users/24756376/data/Flamingo/L1000N0900/Trees/"
def build_tree():
   print("starting")
   file=dir#+str(k)+'.hdf5'
   f = h5py.File(file, 'r')
   halos = np.array([f['center_x'],f['center_y'],f['center_z']]).T
   f.close()
   Tree=KDTree(halos)
   
   print(tracemalloc.get_traced_memory())
   joblib.dump(Tree, data_dir+'halos.pkl')

   print("Done")


def query_tree():
   print("starting")
   Tree=joblib.load(data_dir+'halos.pkl')
   file=dir#+str(k)+'.hdf5'
   f = h5py.File(file, 'r')
#   r200=np.array(f['r200'])
   mass=np.array(f['mass'])
   halos = np.array([f['center_x'],f['center_y'],f['center_z']]).T
   f.close()
   indexes=Tree.query_ball_point(
      centers,
      r=5*r200,
      return_sorted=True,
      
      
   )
   index=np.zeros(len(indexes))
   D=np.zeros(len(indexes))
   for i in range(0,len(indexes)):
      sunmass=mass[indexes[i][1:]]
#      print(sunmass)
      subin=np.array(indexes[i])[1:][sunmass>0.1*mass14[i]]
      if len(subin)==0:
         index[i]=-1
         D[i]=0
      else:
         index[i]=subin[0]
         D[i]=np.sqrt((centers[i]-halos[subin[0]])@(centers[i]-halos[subin[0]]))
        
   joblib.dump([index,D], "/Users/24756376/data/Flamingo/L1000N0900/Trees/halo.joblib")
   return index,D
   



   
   print(tracemalloc.get_traced_memory())
tracemalloc.start()

f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/halos_ranked.hdf5", 'r')
centers=np.array([f['centers_x'],f['centers_y'],f['centers_z']]).T
ids=np.array(f['id'])
centers = centers[ids <= 0]
r200=np.array(f['r200'])[ids <= 0]
mass14=np.array(f['mass'])[ids <= 0]
f.close()

i,d=query_tree()

tracemalloc.stop()
d=d/r200
print(np.max(d))

f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/neighbour.hdf5", 'w')
f.create_dataset("index",data=i)
f.create_dataset("distance",data=d)
f.close()
