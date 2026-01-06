#cosma mass profile/S
from scipy.spatial import KDTree
import h5py
import glob
from tqdm import tqdm
import numpy as np
import joblib
import tracemalloc
from multiprocessing import Pool
import sys
def radial_distance(x, y,z):
    return np.sqrt(x**2 + y**2+z**2)

def spherical_harmonic_0(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(5/16/np.pi)*(3*z**2-r**2)/r**2
def spherical_harmonic_1(x, y,z):
    r = radial_distance(x, y,z)  
    return np.sqrt(15/4/np.pi)*x*z/r**2
def spherical_harmonic_2(x, y,z): 
    r = radial_distance(x, y,z)
    return np.sqrt(15/16/np.pi)*(x**2-y**2)/r**2
def spherical_harmonic__1(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(15/4/np.pi)*y*z/r**2 
def spherical_harmonic__2(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(15/4/np.pi)*x*y/r**2 
def quadrupole(x, y,z):
    r = radial_distance(x, y,z)
    
    f0=np.sum(spherical_harmonic_0(x, y,z)*r) 
    f1=np.sum(spherical_harmonic_1(x, y,z)*r) 
    f2=np.sum(spherical_harmonic_2(x, y,z)*r)
    f_1=np.sum(spherical_harmonic__1(x, y,z)*r)
    f_2=np.sum(spherical_harmonic__2(x, y,z)*r)
    return np.array(f0,f1,f2,f_1,f_2)
def check_boundry(Coord):#do this only after recentering

    x=Coord[:,0]
    y=Coord[:,1]
    z=Coord[:,2]
    if len(x)>0:
      if np.max(x)-np.min(x)>500:
        if np.min(x)<-500:
            Coord[Coord[:,0]<-500,0]+=1000
        elif np.max(x)>500:
            Coord[Coord[:,0]>500,0]-=1000
      if np.max(y)-np.min(y)>500:
        if np.min(y)<-500:
            Coord[Coord[:,1]<-500,1]+=1000
        elif np.max(y)>500:
            Coord[Coord[:,1]>500,1]-=1000
      if np.max(z)-np.min(z)>500:
        if np.min(z)<-500:
            Coord[Coord[:,2]<-500,2]+=1000
        elif np.max(z)>500:
            Coord[Coord[:,2]>500,2]-=1000

    return Coord
def build_tree(k):

   print("starting", k)
   file=path+"snapshots/flamingo_00"+str(snap)+"/flamingo_00"+str(snap)+"."+str(k)+'.hdf5'
   f = h5py.File(file, 'r')
  

   print(tracemalloc.get_traced_memory())
  


   dm = np.asarray(f['PartType1']['Coordinates'],dtype=np.float32)
   mass_dm=np.asarray(f['PartType1']['Masses'])
   
   Tree_dm=KDTree(dm,boxsize=1000.00001,leafsize=256)
   dm=[]
   index_dm=Tree_dm.query_ball_point(
      centers,
      r=3*r50,
      
    

      
   )
   #mass profile
   #quadrupole
   #number of particles
   dm_index=[]
   
   
   for i in tqdm(range (0,len(index_dm))):
      if len(index_dm[i])>0:
        
       

         coord=Tree_dm.data[index_dm[i]].astype(np.float32)
         coord=check_boundry(coord-centers[i])
         r=np.sqrt(np.sum(coord**2,axis=1))
         profile_dm[i]+=np.histogram(r,bins=bins,weights=mass_dm[index_dm[i]])[0]
         coord=coord[r<=1]-com_dm
         f_dm[i]+=quadrupole(coord[:,0],coord[:,1],coord[:,2])
         
         dm_index.append([i,len(index_dm[i])])
         
       
     

   
   index_dm=np.concatenate(index_dm).astype(np.int32)
   dm_index=np.array(dm_index)
   Tree_dm=[]
   print("dm Done")
   mass_g=np.asarray(f['PartType0']['Masses'])
   Temperature=np.asarray(f['PartType1']['Temperatures'])
   g=np.asarray(f['PartType0']['Coordinates'],dtype=np.float32)
   Tree_g=KDTree(g,boxsize=1000.00001,leafsize=256)
   g=[]
   f.close()
   index_g=Tree_g.query_ball_point(
      centers,
      r=3*r50,
     
    
      
   )
   g_index=[]
   
   
   for i in tqdm(range (0,len(index_g))):
      if len(index_g[i])>0:
        
       

         coord=Tree_g.data[index_g[i]].astype(np.float32)
         coord=check_boundry(coord-centers[i])
         r=np.sqrt(np.sum(coord**2,axis=1))
         T=Temperature[index_g[i]]
         profile_g[i]+=np.histogram(r,bins=bins,weights=mass_g[index_g[i]])[0]
         profile_hg[i]+=np.histogram(r[T>50000],bins=bins,weights=mass_g[index_g[i]][T>50000])[0]
         coord=coord[r<=1]-com_gas
         f_g[i]+=quadrupole(coord[:,0],coord[:,1],coord[:,2])
         
         g_index.append([i,len(index_g[i])])
         
       
     

   
   index_g=np.concatenate(index_g).astype(np.int32)
   g_index=np.array(g_index)
   Tree_g=[]
   f=h5py.File(save+str(k)+'.hdf5', 'w')

   f.create_dataset("index_dm", data=index_dm)   
   f.create_dataset("halos_dm",data=dm_index)
   f.create_dataset("halos_g",data=g_index)
   f.create_dataset("index_g", data=index_g)
  
   f.close()
   

 



def query_s(k):
   file=path+"snapshots/flamingo_0077/flamingo_0077."+str(k)+'.hdf5'
   f = h5py.File(file, 'r')
  

   print(tracemalloc.get_traced_memory())
  
   s = np.asarray(f['PartType4']['Coordinates'])
   Tree_s=KDTree(s,boxsize=1000,leafsize=256)
   s=[]
   index_s=Tree_s.query_ball_point(

      centers,
      r=3*r50,
     

    
      
   )
   member_s=[]
   for i in tqdm(range (0,len(index_s))):
      if len(index_s[i])>0:
         member_s.append(np.ones(len(index_s[i]))*i)


   
   index_s=np.concatenate(index_s,dtype=np.int32).astype(np.int32)
   member_s=np.concatenate(member_s).astype(np.int32)
   Coord_s=Tree_s.data[index_s]
   Tree_s=[]
   print("star Done")

   

tracemalloc.start()
path="/Users/24756376/data/Flamingo/L1000N0900/test/flamingo_0077/"
data_dir="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"
snap=sys.argv[1]    
    
save=data_dir+"flamingo_00"+str(snap)+"/"
    
f=h5py.File(save+"halos_central_12.hdf5", 'r')
centers=np.asarray(f['centers'])

r50=np.asarray(f['r50'])
    
com_dm=0
com_gas=0

f.close()


profile_dm=np.zeros((len(centers),50))
bins=np.linspace(0,3,51)
f_dm=np.zeros(len(centers),5)

profile_g=np.zeros((len(centers),50))
profile_hg=p.zeros((len(centers),50))#hot gas
f_g=np.zeros(len(centers))


for k in range(0,1):
  build_tree(int(k))
  if k//5==0 or k==63:
     f=h5py.File(save+"halo_properties.hdf5", 'w')
     dm=f.create_group("PartType1")
     g=f.create_group("PartType0")
     dm.create_dataset("mass_profile",data=profile_dm)
     dm.create_dataset("harmonics",data=f_dm)
     g.create_dataset("mass_profile",data=profile_g)
     g.create_dataset("mass_profile_hot",data=profile_hg)
     g.create_dataset("harmonics",data=f_g)
     f.close()
      

   
 
    

tracemalloc.stop()
