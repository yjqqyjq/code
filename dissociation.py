# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio import load
import swiftgalaxy as sg
import functions as fn
from swiftgalaxy import SWIFTGalaxy, MaskCollection


#work with SOAP
#load the data
soap_dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/"

data_h=sw.load(soap_dir+"/halo_properties_0127.hdf5")
#select the main halo
host_id=data_h.soap.host_halo_index#central halo=-1\
halo_id=np.arange(0,len(host_id),1)
mass=data_h.spherical_overdensity_200_crit.total_mass
mainhalo_id=halo_id[(host_id==-1)*(mass>100)]

data_h=[]

x_dm=[[]]
y_dm=[[]]
z_dm=[[]]
x_g=[[]]
y_g=[[]]
z_g=[[]]
for i in range(0, len(mainhalo_id)-1):
  
    x_dm.append([])
    x_g.append([])
    z_dm.append([])
    y_dm.append([])
    y_g.append([])
    z_g.append([])


#get the coordinates of the dark matter and gas particles

#analyse the main halo
sgs=sg.SWIFTGalaxies(soap_dir+"colibre_with_SOAP_membership_0127.hdf5",
    sg.SOAP(soap_dir+"/halo_properties_0127.hdf5",soap_index=mainhalo_id,extra_mask=None),
    preload={"dark_matter.cartesian_coordinates","gas.cartesian_coordinates"})

i=0
for sgi in sgs:
    fn.analyse(sgi,i,x_dm, y_dm, z_dm,x_g, y_g,z_g)
    i+=1  
print(x_dm[0][0])  
#calculate the dissociation
'''
S=np.zeros(len(x_dm))
for i in  range(0,len(x_dm)):
      S[i]=fn.dissociation(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
print(S)
'''
Offset=np.zeros(len(x_dm))
for i in  range(0,len(x_dm)):
      Offset[i]=fn.offset(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
print(Offset)
#plot
import matplotlib.pyplot as plt
plt.close()
'''
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(S, bins=100)
ax.set_xlabel("S")
ax.set_ylabel("Counts")
ax.set_title("M>10^10")
fig.savefig("/home/jyang/plot/Colibre/L0012N0094/S.png")


plt.close()
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Offset*1000, bins=100)
ax.set_xlabel("Offset(kpc)")
ax.set_ylabel("Counts")
ax.set_title("M>10^10")
fig.savefig("/home/jyang/plot/Colibre/L0012N0094/Offset_central.png")
'''