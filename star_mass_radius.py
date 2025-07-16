import numpy as np
import unyt
import functions as fn
import h5py
from tqdm import tqdm
import matplotlib.pyplot as plt
plt.close()
id=fn.halo_ids
fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.scatter(fn.ms100[id<=0],fn.ms3000[id<=0],s=0.5,alpha=0.3,label="centrals")
ax.scatter(fn.ms100[id>0],fn.ms3000[id>0],s=0.5,alpha=0.3,label="satellites")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
fig.savefig("/Users/24756376/plot/Flamingo/L1000N0900/ms_radius.png")
print(id[(id>0)*(fn.ms3000>(4.5*fn.ms100))])