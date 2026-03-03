#Half mass z
import numpy as np
import unyt
import swiftsimio as sw
#import Falco_S as fd
import matplotlib.pyplot as plt

import h5py
from tqdm import tqdm
from multiprocessing import Pool
import subprocess
import sys
def load_mass(i):
     f=h5py.File(dir+str(i)+".hdf5")
     print(i)
     track=np.array(f["InputHalos"]["HBTplus"]['TrackId'],dtype=int)

     in_z0=np.isin(track,track_z0_ranked,assume_unique=True)
     track_z=track[in_z0]
     args_z=np.argsort(track_z)
     track_z_ranked=track_z[args_z]

     in_z=np.isin(track_z0_ranked,track_z_ranked,assume_unique=True)
     

     m200=np.array(f["SO"]['200_crit']['TotalMass'])[in_z0][args_z]
     m200_z=np.where(in_z,m200,0)

     f.close()
     
     return m200_z
def test(i):
     print(i)
     f=h5py.File(dir+str(i*2-1)+".hdf5")

     halo_ids=np.array(f["InputHalos"]["HBTplus"]['TrackId'])
     host_ids=np.array(f["BoundSubhalo"]["TotalMass"])
     print(host_ids[halo_ids==id])
     f.close()
     

    
     
def load_redshift(i):#To z=3
    data=sw.load(dir+str(i)+".hdf5")
    z=data.metadata.redshift
    print(z)

if __name__ == '__main__':



    dir="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_00"

    save="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"


    f=h5py.File(save+"flamingo_0077/halos_central_12.hdf5","r")


    track_z0=np.array(f['track_id'],dtype=int)
    m200_z0=np.array(f['m200'],dtype=int)


    f.close()
 
    args_z0=np.argsort(track_z0)
    args_rev_z0=np.argsort(args_z0)
    track_z0_ranked=track_z0[args_z0]
    
    with Pool(4) as p: # Create a pool with 4 worker processes





       halo_mass=p.map(load_mass, range(76,77))
    halo_mass=np.array(halo_mass)[:,args_rev_z0]


    print(halo_mass[0])



    f=h5py.File(save+"mass_evolution_12.hdf5","w")
    

#    del f["m50"]
    f.create_dataset("m200",data=halo_mass)

    f.close()
       

