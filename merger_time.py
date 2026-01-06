#merger time
import h5py
import numpy as np
from tqdm import tqdm
data_path="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"
snap_path="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/"
#z=0
f=h5py.File(data_path+'flamingo_0077/halos_central_12_idorder.hdf5','r')
track_main=np.asarray(f["track_id"])
mass_main=np.asarray(f["mass"])
halo_id_main=np.asarray(f["id"])
z_born=np.asarray(f["zborn"])
f.close()
main_id=np.arange(0,len(mass_main),1)
f=h5py.File(data_path+'flamingo_0077/halos_satellites_12.hdf5','r')
track_id0=np.asarray(f["track_id"])
parent_track0=np.asarray(f["parent_track"])
z_born0=np.asarray(f["z_born"])
f.close()
halo_left=np.full(len(track_id0),True)
snapshots=np.arange(37,76,1)[::-1]
i=0
snap_born=np.zeros(len(track_id0))

for snap in tqdm(snapshots):
    halo_left=(snap_born>=i+1)
    f=h5py.File(data_path+"flamigno_00"+str(snap).zfill(2)+'/halos_satellite_12_z0.hdf5','r')
    mass=np.asarry(f["mass"])
    track_id=np.asarray(["track_id"])#auusme it's in the same order as z=0
    host_id=np.asarray(f["host_id"])
   
    parent_track=np.array(f["parent_track"])
    f.close()
    mask=np.isin(track,track_id,assume_unique=True)
    #Find what's left
   
    #order the track_id
    
    host_id=np.asarray(f["SOAP"]["HostHaloIndex"])[mask]
    mass=np.asarray(f['BoundSubhalo']['TotalMass'])[mask]
    f.close()
    parent_track=np.where(host_id==-1,track[mask],track[host_id])
    track=track[mask]
    
    halo_lost=np.isin(track_id[halo_left],track,assume_unique=True,invert=True)
    snap_born[halo_left][halo_left]=i+1
    halo_left[halo_left][halo_lost]=False
    i=i+1
    f=h5py.File(data_path+"flamigno_00"+str(snap).zfill(2)+'/halos_satellite_12_z0.hdf5','w')
    f.create_dataset("mass", data=mass)
    f.create_dataset("track_id", data=track)
    f.create_dataset("host_id", data=host_id)
   
    f.create_dataset("host_track",data=parent_track)
    f.close()

    
   
    if len(track)>0:
        
        host_done=host_track0[mask]
        mask_track*=np.isin(host_done,host_track0,invert=True,assume_unique=True)
        
        host_ids=track_host0[track]
        ids=host_halo0[host_ids]
        if np.max(z_m[ids])>-1:
            print("calculating the same halo")
        z_m[ids]=i*0.05
        #mass
        
        mass_ratios[ids]=mass_ratio
        i=i+1
        
    
