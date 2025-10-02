#Kdtree
from sklearn.neighbors import KDTree
import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool

dir="/Users/24756376/data/Flamingo/L1000N0900/particles_ranked.hdf5"
data_dir="/Users/24756376/data/Flamingo/L1000N0900/Trees/"
def build_tree(k):
   print("starting", k)
   file=dir#+str(k)+'.hdf5'
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
      r=6,
      return_distance=True,
      sort_results=True,
      
   )
   joblib.dump([index,D], "/Users/24756376/data/Flamingo/L1000N0900/Trees/"+str(k)+'dm.joblib')

   
   Tree_dm=[]

   print(tracemalloc.get_traced_memory())
   Tree_g=joblib.load(data_dir+str(k)+'g.pkl')
   index,D=Tree_g.query_radius(
        centers,
        r=6,
        return_distance=True,
        sort_results=True,
        
    )
   joblib.dump([index,D], "/Users/24756376/data/Flamingo/L1000N0900/Trees/"+str(k)+'g.joblib')
   Tree_g=[]
   
   print(tracemalloc.get_traced_memory())
'''
if __name__ == '__main__':
    tracemalloc.start()
    f=h5py.File("/home/jyang/data/Flamingo/L1000N0900/halos_ranked.hdf5", 'r')
    centers=np.array([f['centers_x'],f['centers_y'],f['centers_z']]).T
    ids=np.array(f['id'])
    centers = centers[ids <= 0]
    print(len(centers))
    f.close()

    with Pool() as p:
        p.map(query_tree, range(0,32))
    tracemalloc.stop()
'''
f=h5py.File("/Users/24756376/data/Flamingo/L1000N0900/halos_ranked.hdf5", 'r')
centers=np.array([f['centers_x'],f['centers_y'],f['centers_z']]).T
ids=np.array(f['id'])
centers = centers[ids <= 0]
r100=np.array(f['r100'])[ids <= 0]
print(len(centers))
f.close()
file="/Users/24756376/data/Flamingo/L1000N0900/Trees/halos.joblib"
index,D=joblib.load(file)
D=D[0:11965]/r100
print(np.histogram(D,bins=20))
print(np.argwhere(D<1))
print(D[np.argwhere(D<1)])
