#find com
from scipy.spatial import KDTree
import h5py
import glob
import datetime
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
def quadrupole(x, y,z,m):
    r = radial_distance(x, y,z)

    f0=np.sum(spherical_harmonic_0(x, y,z)*r*m)

    f1=np.sum(spherical_harmonic_1(x, y,z)*r*m)
    f2=np.sum(spherical_harmonic_2(x, y,z)*r*m)
    f_1=np.sum(spherical_harmonic__1(x, y,z)*r*m)
    f_2=np.sum(spherical_harmonic__2(x, y,z)*r*m)
    return np.array([f0,f1,f2,f_1,f_2])
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
def check_cross(Coord):#do this only after recentering

    x=Coord[:,0]
    y=Coord[:,1]
    z=Coord[:,2]
    if len(x)>0: 
      if np.max(x)-np.min(x)>500:
        return 1
      if np.max(y)-np.min(y)>500:
        return 1
      if np.max(z)-np.min(z)>500:
        return 1

    return 0
def load_tree(k):

   print("starting", k)
   file=path+"snapshots/flamingo_00"+str(snap)+"/flamingo_00"+str(snap)+"."+str(k)+'.hdf5'
   f=h5py.File(save+str(k)+'.hdf5', 'r')

   index_dm=np.asarray(f["index_dm"])
   halos_id_dm=np.asarray(f["halos_dm"])
   halos_id_g=np.asarray(f["halos_g"])
   index_g=np.asarray(f["index_g"])
   f.close()
   f = h5py.File(file, 'r')


   print(tracemalloc.get_traced_memory())



   dm = np.asarray(f['PartType1']['Coordinates'],dtype=np.float32)
   mass_dm=np.asarray(f['PartType1']['Masses'])
   
   for i in tqdm(range(len(halos_id_dm))):
      
      id=halos_id_dm[:,0][i]
      i_ini=np.sum((halos_id_dm[:,1][0:i]))
      i_end=np.sum((halos_id_dm[:,1][0:i+1]))
      index=index_dm[i_ini:i_end].astype(int)
      



      coord=dm[index]
      if check_cross(coord):
        coord=check_boundry(coord-centers[id])
        r=np.sqrt(np.sum(coord**2,axis=1))
      
        coord=coord[r<r50[id]]#particles within r50
        com_dm[id]=np.average(coord,axis=0,weights=mass_dm[index][r<r50[id]])
      
   mass_dm=[]
   dm=[]
   mass_g=np.asarray(f['PartType0']['Masses'])

   g=np.asarray(f['PartType0']['Coordinates'],dtype=np.float32)
   for i in tqdm(range(len(halos_id_g))):
      
      id=halos_id_g[:,0][i]
      i_ini=np.sum((halos_id_g[:,1][0:i]))
      i_end=np.sum((halos_id_g[:,1][0:i+1]))
      index=index_g[i_ini:i_end].astype(int)
      coord=g[index]
      
      if check_cross(coord):
        coord=check_boundry(coord-centers[id])
        r=np.sqrt(np.sum(coord**2,axis=1))
      
        coord=coord[r<r50[id]]#particles within r50
        com_gas[id]=np.average(coord,axis=0,weights=mass_g[index][r<r50[id]])
   f.close()
def load_star(k):
   print("starting", k)
   #read the index and halos id
   file=path+"snapshots/flamingo_00"+str(snap)+"/flamingo_00"+str(snap)+"."+str(k)+'.hdf5'
   f=h5py.File(save+str(k)+'stars.hdf5', 'r')

   index_s=np.asarray(f["index_s"])
   halos_id_s=np.asarray(f["halos_s"])

   f.close()
   #read coordinate and mass
   f = h5py.File(file, 'r')


   print(tracemalloc.get_traced_memory())



   s = np.asarray(f['PartType4']['Coordinates'],dtype=np.float32)
   mass_s=np.asarray(f['PartType4']['Masses'])
   #load coordiantes of single halos and calculate the quadrupole moment
   for i in tqdm(range(len(halos_id_s))):
      id=halos_id_s[:,0][i]
      i_ini=np.sum((halos_id_s[:,1][0:i]))
      i_end=np.sum((halos_id_s[:,1][0:i+1]))
      index=index_s[i_ini:i_end].astype(int)
      



      coord=s[index]
      coord=check_boundry(coord-centers[id])
      r=np.sqrt(np.sum(coord**2,axis=1))
      
      coord=coord[r<r50[id]]#particles within r50
      coord=coord-com[id]
      mass=mass_s[index][r<r50[id]]
      r=np.sqrt(np.sum(coord**2,axis=1))
      
      f_s[i]+=quadrupole(coord[:,0],coord[:,1],coord[:,2],mass)
      r_s[i]+=np.sum(np.sqrt(np.sum(coord**2,axis=1))*mass)
      
tracemalloc.start()
path="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/"
data_dir="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"
snap=sys.argv[1]

save=data_dir+"flamingo_00"+str(snap)+"/"
#print(save)    
f=h5py.File(save+"halos_central_12.hdf5", 'r')
centers=np.asarray(f['centers'])

r50=np.asarray(f['r50'])



com_dm=np.asarray(f['com_dm'])
com_gas=np.asarray(f['com_g'])
com_s=np.asarray(f['com_s'])

M_g=np.asarray(f['M_g'])
M_s=np.asarray(f['M_s'])
com_b=((M_s*com_s.T+M_g*com_gas.T)/(M_g+M_s)).T
#com_b=com_gas
f.close()
com=0.5*(com_dm+com_b)
#print(centers)

time=datetime.datetime.now()
profile_dm=np.zeros((len(centers),50))
bins=np.linspace(0,3,51)
f_dm=np.zeros((len(centers),5))

profile_g=np.zeros((len(centers),50))
profile_hg=np.zeros((len(centers),50))#hot gas
f_g=np.zeros((len(centers),5))
f_s=np.zeros((len(centers),5))
r_dm=np.zeros(len(centers))
r_g=np.zeros(len(centers))
r_s=np.zeros(len(centers))
for k in range(0,64):
  load_star(int(k))
  if k//10==0 or k==63:
    if k//10==0:
       continue
#      f=h5py.File(save+"halo_properties"+str(k)+".hdf5", 'a')
#      dm=f.create_group("PartType1")
#      g=f.create_group("PartType0")
#      s=f.create_group("PartType2")
    elif k==63:
      f=h5py.File(save+"halo_properties_stars.hdf5", 'a')
      s=f["PartType2"]
    s.create_dataset("f_s",data=f_s)
    s.create_dataset("r_s",data=r_s)
    f.close()
'''
      dm=f["PartType1"]
      g=f["PartType0"]
    dm.create_dataset("r_dm",data=r_dm)
    g.create_dataset("r_g",data=r_g)

    if "harmonics" in dm:
      dm["harmonics"]=f_dm
    else:
      dm.create_dataset("harmonics",data=f_dm)
    if "harmonics" in g:
      g["harmonics"]=f_g
    else:
      g.create_dataset("harmonics",data=f_g)
    f.close()
'''
tracemalloc.stop()