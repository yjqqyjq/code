#halo born time
import h5py
import numpy as np
from tqdm import tqdm
data_path="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/"
snap_path="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/"
#z=0
f=h5py.File(data_path+'flamingo_0077/halos_central_12_idorder.hdf5','r')
track_id=np.asarray(f["track_id"])
halo_left=np.full(len(track_id),True)
f.close()
z_born=np.ones(len(track_id))*2
snapshots=np.arange(37,76,1)[::-1]
i=0
for snap in tqdm(snapshots):
    f=h5py.File(data_path+'flamingo_00'+str(snap)+'/halos_central_12_idorder.hdf5','r')
    track_main=np.asarray(f["track_id"])
    
  
    f.close()
    mask=np.isin(track_id[halo_left],track_main,invert=True,assume_unique=True)
    z_born[halo_left][mask]=i*0.05
    print(len(z_born[halo_left][mask]))
    halo_left[halo_left][mask]=False
    
    if np.min(z_born[mask])<2:
        print("error")
    i=i+1

f=h5py.File(data_path+'flamingo_0077/halos_central_12_idorder.hdf5','a')
f.create_dataset("z_born",data=z_born)

f.close()