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
     f=h5py.File(dir+str(i*2-1)+".hdf5")

     halo_ids=np.array(f["InputHalos"]["HBTplus"]['TrackId'])
     mask=np.isin(halo_ids,ids,assume_unique=True)
     m50=np.array(f["SO"]['50_crit']['TotalMass'])[mask]
     halo_ids=halo_ids[mask]
     mask_id=np.argsort(halo_ids)[args_rev]
     m50=m50[mask_id]
     return m50
def test(i):
     f=h5py.File(dir+str(i*2-1)+".hdf5")
     halo_ids=np.array(f["InputHalos"]["HBTplus"]['TrackId'])
     host_ids=np.array(f["InputHalos"]["HostHaloIndex"])
     print(host_ids[halo_ids==id])
     f.close()
     
    
     
def load_redshift(i):#To z=3
    data=sw.load(dir+str(i)+".hdf5")
    z=data.metadata.redshift
    print(z)

if __name__ == '__main__':
    dir="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_00"

    save="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"
    f=h5py.File(save+"central_halos.hdf5","r")

    ids=np.array(f['Track_id'])
    id=ids[2]

    args_id=np.argsort(ids)
    args_rev=np.argsort(args_id)
    f.close()
    with Pool(1) as p: # Create a pool with 5 worker processes


       halo_mass=p.map(test, range(25,40))
    halo_mass=np.array(halo_mass)

    print(halo_mass)
    f=h5py.File(save+"mass_evolution.hdf5","w")
    f.create_dataset("m50",data=halo_mass)
    f.close()
       
