#reduce snapshots

import numpy as np
import h5py
import tracemalloc
from pathlib import Path
from tqdm import tqdm
import datetime

from multiprocessing import Pool
from functools import partial
import os
print("start")
print(datetime.datetime.now())
dir="../../../mnt/su3-pro/L2800N5040/snapshots/"
#save_dir="../../../mnt/su3-pro/L2800N5040/snapshots_reduced/"
#soap_dir="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
tracemalloc.start()
#dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/colibre_with_SOAP_membership_0127.hdf5"
#soap_dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"

         
def reduce_snapshots(k):
   
    
    
    data=h5py.File(dir+"flamingo_0078."+str(k)+".hdf5",'a')

    #dm particles
    #KeysViewHDF5 ['Coordinates', 'FOFGroupIDs', 'Masses', 'ParticleIDs', 'Softenings', 'Velocities']>
    dm=data["PartType1"]
    for key in list(dm.keys()):
      #Entries I delete, below is the same
      if np.isin(key,['Softenings','Masses','FOFGroupIDs']):
         
         del dm[key]
      else:
         print(key)
    print(tracemalloc.get_traced_memory())
    #gas particles
    #['ComptonYParameters', 'Coordinates', 'Densities', 'ElectronNumberDensities', 'FOFGroupIDs', 'LastAGNFeedbackScaleFactors', 'Masses', 'MaximalTemperatureScaleFactors', 'MaximalTemperatures', 'MetalMassFractions', 'ParticleIDs', 'Pressures', 'SmoothedElementMassFractions', 'SmoothedMetalMassFractions', 'SmoothingLengths', 'StarFormationRates', 'Temperatures', 'Velocities', 'VelocityDivergences', 'XrayLuminosities', 'XrayPhotonLuminosities']
    g=data["PartType0"]
    for key in list(g.keys()):
       if np.isin(key,['ElectronNumberDensities','SmoothedElementMassFractions', 'SmoothingLengths', 'ComptonYParameters', 
                            'VelocityDivergences', 'XrayPhotonLuminosities','Masses','MaximalTemperatures','FOFGroupIDs']):
         
         del g[key]
       else:
         print(key)
    print(tracemalloc.get_traced_memory())
    #star particles
    #['BirthDensities', 'BirthScaleFactors', 'Coordinates', 'FOFGroupIDs', 'InitialMasses', 'Luminosities', 'Masses', 'MetalMassFractions', 'ParticleIDs', 'SmoothedElementMassFractions', 'SmoothedMetalMassFractions', 'SmoothingLengths', 'Velocities']
    s=data["PartType4"]
    for key in list(s.keys()):
       if np.isin(key,['SmoothingLengths', 'BirthDensities', 'SmoothedElementMassFractions','Luminosities','FOFGroupIDs']):
         
         del s[key]
       else:
         print(key)
    print(tracemalloc.get_traced_memory())
    #BH particles
    #['AccretionRates', 'Coordinates', 'DynamicalMasses', 'FOFGroupIDs', 'FormationScaleFactors', 'LastAGNFeedbackScaleFactors', 'LastHighEddingtonFractionScaleFactors', 'LastMajorMergerScaleFactors', 'LastMinorMergerScaleFactors', 'NumberOfAGNEvents', 'NumberOfHeatingEvents', 'NumberOfMergers', 'ParticleIDs', 'SmoothingLengths', 'SubgridMasses', 'TotalAccretedMasses', 'Velocities']
    bh=data["PartType5"]
    for key in list(bh.keys()):
        if np.isin(key, ["Coordinates", 'DynamicalMasses', 'AccretionRates' , 'ParticleIDs', 'LastAGNFeedbackScaleFactors', 'NumberOfAGNEvents', 'NumberOfMergers']):
          print(key)
        else:
          del bh[key]
    print(tracemalloc.get_traced_memory())
    del data['PartType6']
    data.close()
     

#    potential_gas=potential_gas[pmask2]
    
    

    print(datetime.datetime.now())
if __name__ == '__main__':

  tracemalloc.start()
  print("start")  
  with Pool() as p: # Create a pool 
        results = p.map(reduce_snapshots, range(700,800))
  tracemalloc.stop()