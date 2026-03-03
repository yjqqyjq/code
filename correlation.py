#correlation
import numpy as np
from Corrfunc.theory.xi import xi
import h5py
import swiftsimio as sw
import matplotlib.pyplot as plt
file='/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_ADIABATIC/SOAP-HBT/flamingo_0077.hdf5'
mask=sw.mask(file)
boxsize = mask.metadata.boxsize
load_region = [[0.0 * b, 0.2 * b] for b in boxsize]
mask.constrain_spatial(load_region)
data=sw.load(file,mask=mask)
x_dm=data.dark_matter.coordinates[:,0]
y_dm=data.dark_matter.coordinates[:,1]
z_dm=data.dark_matter.coordinates[:,2]
x_gas=data.gas.coordinates[:,0]
y_gas=data.gas.coordinates[:,1]
z_gas=data.gas.coordinates[:,2]

h=0.6774
boxsize=0.2*boxsize[0]*h #Mpc/h
nthreads = 2
rmin = 0.1
rmax = 10.0
nbins = 10
rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins)
xi_dm = xi(boxsize, nthreads, rbins, x_dm,y_dm,z_dm)
xi_gas = xi(boxsize, nthreads, rbins, x_gas,y_gas,z_gas)
f=h5py.File('/cosma8/data/do012/dc-yang9/Flamingo/L1000N1800_NoCool/correlation_functions.hdf5','w')
f.create_dataset('dm',data=xi_dm)
f.create_dataset('gas',data=xi_gas)
f.create_dataset('rbins',data=rbins)
f.close()