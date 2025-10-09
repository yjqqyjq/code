#load cell
import swiftsimio as sw
import numpy as np
import h5py
from multiprocessing import Pool
import tracemalloc

def load_cell(i):
  center=centers[i]
  length=cell_length
  mask = sw.mask(data_dir)
  load_region=[[center[0]-length/2,center[0]+length/2],[center[1]-length/2,center[1]+length/2],[center[2]-length/2,center[2]+length/2]] 
  mask.constrain_spatial(load_region)
  sw.subset_writer.write_subset("/data/L1000N1800/Nocool/cell"+str(int(i))+".hdf5", mask)
  print(tracemalloc.get_traced_memory())
  
if __name__ == "__main__":
  tracemalloc.start()
  data_dir="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_ADIABATIC/SOAP-HBT/flamingo_0077.hdf5"
  
  dir="/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_ADIABATIC/SOAP-HBT/halo_properties_0077.hdf5"
  f=h5py.File(data_dir,'r')
  cell_length=float(f["Cells"]["Centres"][0][0]*2)
  cell_center=np.array(f["Cells"]["Centres"])
  cell_counts=np.array(f["Cells"]["Counts"])
  cell_num=32
 

  f.close()
  data_h=sw.load(dir)

  host_id=np.array(data_h.soap.host_halo_index)#central halo=-1\

  mass=np.array(data_h.bound_subhalo.total_mass)

  mask=(host_id==-1)*(mass>10000)
  centers=np.array(data_h.input_halos.halo_centre)[mask]
  
  n_cell=centers[:,2]//cell_length+(centers[:,1]//cell_length)*cell_num+(centers[:,0]//cell_length)*cell_num*cell_num
  n_cell=np.unique(n_cell).astype(int)#cell list
  centers=centers[n_cell]
  '''
  with Pool(32) as p: # Create a pool with 5 worker processes
    p.map(load_cell, range(0,len(n_cell)))
  '''
  tracemalloc.stop()