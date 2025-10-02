#load cell
import swiftsimio as sw
import numpy as np
import h5py


def load_cell(i):
    
  return 0

if __name__ == "__main__":
  data_dir="/Users/24756376/data/Flamingo/L1000N0900/test/flamingo_0077.0.hdf5"
  
  dir="/Users/24756376/data/Flamingo/L1000N0900/halo_catalogue/halo_properties_0077.hdf5"
  f=h5py.File(data_dir,'r')
  cell_length=float(f["Cells"]["Centres"][0][0]*2)
  cell_center=np.array(f["Cells"]["Centres"])
  cell_num=32
 

  f.close()
  data_h=sw.load(dir)

  host_id=np.array(data_h.soap.host_halo_index)#central halo=-1\

  mass=np.array(data_h.bound_subhalo.total_mass)

  mask=(host_id==-1)*(mass>10000)
  centers=np.array(data_h.input_halos.halo_centre)[mask]
  n_cell=centers[:,2]//cell_length+(centers[:,1]//cell_length)*cell_num+(centers[:,0]//cell_length)*cell_num*cell_num
  n_cell=np.unique(n_cell).astype(int)#cell list
  print(len(np.unique(n_cell)))