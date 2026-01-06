#track_halo
import h5py
import numpy as np
from tqdm import tqdm
data_path="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"
snap_path="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/"
#z=0


f=h5py.File(data_path+'flamingo_0077/halos_satellites_12.hdf5','r')
track_id=np.asarray(f["track_id"])
host_track=np.asarray(f["parent_track"])
f.close()
halo_left=np.full(len(track_id),True)
snapshots=np.arange(37,76,1)[::-1]
i=0
snap_born=np.zeros(len(track_id))

for snap in tqdm(snapshots):
    dir="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_00"+str(snap).zfill(2)+".hdf5"
    f=h5py.File(dir,'r' )
    track=np.asarray(f["InputHalos"]["HBTplus"]["TrackId"])
    mask=np.isin(track,track_id,assume_unique=True)
    #Find what's left
   
    #order the track_id
    
    host_id=np.asarray(f["SOAP"]["HostHaloIndex"])[mask]
    mass=np.asarray(f['BoundSubhalo']['TotalMass'])[mask]
    f.close()
    parent_track=np.where(host_id==-1,track[mask],track[host_id])
    track=track[mask]
    
    halo_lost=np.isin(track_id[halo_left],track,assume_unique=True,invert=True)
    snap_born[halo_left][halo_left]=i
    halo_left[halo_left][halo_lost]=False
    i=i+1
    f=h5py.File(data_path+"flamigno_00"+str(snap).zfill(2)+'/halos_satellite_12_z0.hdf5','w')
    f.create_dataset("mass", data=mass)
    f.create_dataset("track_id", data=track)
    f.create_dataset("host_id", data=host_id)
   
    f.create_dataset("host_track",data=parent_track)
    f.close()
f=h5py.File(data_path+'flamingo_0077/halos_satellites_12.hdf5','a')
f.create_dataset("snap_born",data=snap_born)
f.close()