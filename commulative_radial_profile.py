#commulative radius
import matplotlib.pyplot as plt
import numpy as np
import h5py
f=h5py.File("/Users/24756376/data/Flamingo/L1000N1800/2d_fgas_H.hdf5",'r')
hdm=np.array(f['Ndm'])
hg=np.array(f['Ng'])
hs=np.array(f['Ns'])
bin=np.array(f['bin'])
binv=np.array(f['binv'])

f.close()
f=h5py.File("/Users/24756376/data/Flamingo/L1000N1800_NoCOol/2d_fgas_H.hdf5",'r')
hdm_ad=np.array(f['Ndm'])
hg_ad=np.array(f['Ng'])

f.close()
hga=np.cumsum(np.sum(hg,axis=2),axis=1)
hdma=np.cumsum(np.sum(hdm,axis=2),axis=1)
hsa=np.cumsum(np.sum(hs,axis=2),axis=1)
hfa=hga+hsa*6.54/8.56
hga_ad=np.cumsum(np.sum(hg_ad,axis=2),axis=1)
hdma_ad=np.cumsum(np.sum(hdm_ad,axis=2),axis=1)
hfa_ad=hga_ad[:,0:9]
hfa=hfa[:,0:9]
hdma=hdma[:,0:9]
hdma_ad=hdma_ad[:,0:9]
'''
for i in range(len(hga)):
    
    hdma[i]=hdma[i]/hdma[:,8][i]
    hfa[i]=hfa[i]/hfa[:,8][i]
   
    hdma_ad[i]=hdma_ad[i]/hdma_ad[:,8][i]
    hfa_ad[i]=hfa_ad[i]/hfa_ad[:,8][i]
'''
print(hdma[0])
bin=bin[0:9]
print(bin)
f=h5py.File("/Users/24756376/data/Flamingo/L1000N1800/cum_profiles.hdf5",'w')
f.create_dataset("r",data=bin)
f.create_dataset("baryons",data=hfa)
f.create_dataset("dm",data=hdma)
f.close()
f=h5py.File("/Users/24756376/data/Flamingo/L1000N1800_NoCool/cum_profiles.hdf5",'w')
f.create_dataset("r",data=bin)
f.create_dataset("baryons",data=hfa_ad)
f.create_dataset("dm",data=hdma_ad)
f.close()