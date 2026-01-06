#load_sub
import h5py
import numpy as np
import sys
k=sys.argv[1]
path="/cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/flamingo_00"+str(k).zfill(2)+"/"
snap_path="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/"
#z=0
'''
f=h5py.File(data_path+'flamingo_0077/halos_central_12.hdf5','r')
track_id=np.asarray(f["track_id"])
f.close()
'''
dir="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_00"+str(k).zfill(2)+".hdf5"
f=h5py.File(dir,'r' )

host_id=np.asarray(f["SOAP"]["HostHaloIndex"])#central halo=-1\
track_id=np.asarray(f["InputHalos"]["HBTplus"]["TrackId"])

mass=np.asarray(f['BoundSubhalo']['TotalMass'])
mask=(mass>10**2.5)*(host_id==-1)
halo_id=np.arange(0,len(host_id),1)
main_id=halo_id[mask]
main_mass=mass[mask]


mask=np.isin(host_id,main_id)
sub_mass=mass[mask]
sub_host=host_id[mask]
sub_track=track_id[mask]
sub_id=halo_id[mask]
mask=np.argsort(sub_host)
sub_mass=sub_mass[mask]
sub_host=sub_host[mask]
sub_track=sub_track[mask]
sub_id=sub_id[mask]
f=h5py.File(path+'halos_satellites_12.hdf5','w')
f.create_dataset("mass", data=sub_mass)
f.create_dataset("track_id", data=sub_track)
f.create_dataset("host_id", data=sub_host)
f.create_dataset("id",data=sub_id)
f.close()