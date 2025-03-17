# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio import load
import swiftgalaxy as sg

def radial_distance(x, y):
    return np.sqrt(x**2 + y**2)

def quadrupole(x, y, num_particles):
    r = radial_distance(x, y)
    
    return np.sqrt(15/(16*np.pi) * ( 4*np.sum(x*y/r)**2 + 
    np.sum(x**2/r - y**2/r)**2 )) / num_particles
   
def dissociation(x_dm, y_dm, x_g, y_g):
    n_dm = len(x_dm)
    n_g = len(x_g)
    r_mean_dm = np.average(radial_distance(x_dm, y_dm))
    r_mean_g = np.average(radial_distance(x_g, y_g))
    r_max = max(r_mean_dm, r_mean_g)
    q_dm = quadrupole(x_dm, y_dm, n_dm)
    q_g = quadrupole(x_g, y_g, n_g)
    return np.sqrt(4*np.pi/5) * (q_dm - q_g) / r_max


#work with SOAP
#load the data
soap_dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/"

data_h=sw.load(soap_dir+"/halo_properties_0127.hdf5")




#select the main halo
host_id=data_h.soap.host_halo_index#central halo=-1\
halo_id=np.arange(0,len(host_id),1)
mass=data_h.spherical_overdensity_200_crit.total_mass
mainhalo_id=halo_id[(host_id==-1)*(mass>100)]

x_dm=[[]]
y_dm=[[]]
x_g=[[]]
y_g=[[]]
for i in range(0, len(mainhalo_id)-1):
  
    x_dm.append([])
    x_g.append([])
    y_dm.append([])
    y_g.append([])


#get the coordinates of the dark matter and gas particles
def analyse(sgi,i):
      x_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates.x)
      y_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates.y)
      x_g[i]=np.array(sgi.gas.cartesian_coordinates.x)
      y_g[i]=np.array(sgi.gas.cartesian_coordinates.y)
#analyse the main halo
for i in range(0,len(mainhalo_id)):
  print(i)
  sga=sg.SWIFTGalaxy(soap_dir+"colibre_with_SOAP_membership_0127.hdf5",
    sg.SOAP(soap_dir+"/halo_properties_0127.hdf5",soap_index=mainhalo_id[i],extra_mask=None))
  analyse(sga,i)
#calculate the dissociation
S=[]
for i in  range(0,len(x_dm)):
      x_dmc=np.sum(x_dm[i])/len(x_dm[i])
      x_gc=np.sum(x_g[i])/len(x_g[i])
      x_dm[i]-=x_dmc
      x_g[i]-=x_gc
      y_dmc=np.sum(y_dm[i])/len(y_dm[i])
      y_dm[i]-=y_dmc
      y_gc=np.sum(y_g[i])/len(y_g[i])
      y_g[i]-=y_gc
      S.append(dissociation(x_dm[i], y_dm[i], x_g[i], y_g[i]))
#S max is 0.866?
#plot
import matplotlib.pyplot as plt
plt.close()
fig = plt.figure()
ax=plt.subplot(1,1,1)
ax.plot()
ax.set_xlabel()
ax.set_ylabel()
fig.savefig("./mass_function.png")

