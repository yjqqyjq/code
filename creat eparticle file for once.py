#creat eparticle file for once
import numpy as np
import unyt


import matplotlib.pyplot as plt

import h5py
from tqdm import tqdm

import subprocess
import sys
f=h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos_ranked.hdf5','r')

ids=np.array(f['id'])
print(f["r100"][ids==-196])
print(f["r50"][ids==-196])
f.close()
f=h5py.File('/Users/24756376/data/Flamingo/L1000N0900/cluster_particles.hdf5','r')
xray=np.array(f['PartType0']["xray_lum_erosita_high"])
coord=np.array(f['PartType0']['Coordinates'])
member=np.array(f['PartType0']['member'])
f.close()
coord=coord[member==member[0]]
xray=xray[member==member[0]]
f=h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halo_xray.hdf5','w')
g=f.create_group('PartType0')
g.create_dataset('Coordinates',data=coord)
g.create_dataset('xray_lum_erosita_high',data=xray)
f.close()