# The code snippet `import numpy as np` imports the NumPy library and assigns it an alias `np` for
# easier reference in the code. NumPy is a popular library in Python used for numerical computations.
import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio import load
import swiftgalaxy as sg
import functions as fn
from swiftgalaxy import SWIFTGalaxy, MaskCollection
import h5py
from tqdm import tqdm


#work with SOAP
#load the data
path="/Users/24756376/data/Flamingo/L1000N0900/"

f=h5py.File(path+'halos_ranked.hdf5','r')
N_g=np.array(f['N_g'])
N_dm=np.array(f['N_dm'])
N_s=np.array(f['N_s'])
N_dm_c=np.array(f['N_dm_c'])
N_g_c=np.array(f['N_g_c'])
N_s_c=np.array(f['N_s_c'])
input_id=np.array(f['input_ids'])
halo_ids=np.array(f['id'])
mass=np.array(f['mass'])
cross=np.array(f['cross_bound'])
f.close()

f=h5py.File(path+'particles_ranked.hdf5','r')
Coord_dm=np.array(f['PartType1']['Coordinates'],dtype=np.float32)
Coord_g=np.array(f['PartType0']['Coordinates'],dtype=np.float32)
Coord_s=np.array(f['PartType2']['Coordinates'],dtype=np.float32)
f.close()

def load_halo(id,dm=0,g=0,s=0):
      """
      Load halo data from the specified path and ID.
      """
   
      arg=np.nonzero(halo_ids==id)[0]
      
    
    
      arg=int(arg)
      slide=[]
      if dm==1:
   
        dm_s=int(np.sum(N_dm[0:arg]))
        dm_e=int(np.sum(N_dm[0:arg+1]))
        slide.append(slice(dm_s,dm_e))
      if g==1:
        g_s=int(np.sum(N_g[0:arg]))
        g_e=int(np.sum(N_g[0:arg+1]))
        slide.append(slice(g_s,g_e))
      if s==1:
        s_s=int(np.sum(N_s[0:arg]))
        s_e=int(np.sum(N_s[0:arg+1]))
        slide.append(slice(s_s,s_e))
      
      slide=np.array(slide).T
      key=[]
      if dm==1:
        key.append('dm')
      if g==1:
        key.append('gas')
      if s==1:
        key.append('stars')
      return key,slide#in the form of [dm,g,s,[[dm_s,dm_e],[g_s,g_e],[s_s,s_e]]]
def load_cluster(id,dm=0,g=0,s=0):
   if id>0:
      id=-int(id) 
      print("warning: this is a satellite, loading the cluster it belongs to")
   arg=np.nonzero(halo_ids[halo_ids<=0]==id)[0]
   arg=int(arg)
   slide=[]
   if dm==1:
     
   
     dm_s=int(np.sum(N_dm_c[0:arg]))
     dm_e=int(np.sum(N_dm_c[0:arg+1]))
     slide.append(slice(dm_s,dm_e))
   #  print(arg)
   if g==1:
      g_s=int(np.sum(N_g_c[0:arg]))
      g_e=int(np.sum(N_g_c[0:arg+1]))
      slide.append(slice(g_s,g_e))
   if s==1:
      s_s=int(np.sum(N_s_c[0:arg]))
      s_e=int(np.sum(N_s_c[0:arg+1]))
      slide.append(slice(s_s,s_e))

   slide=np.array(slide)
   key=[]
   if dm==1:
       key.append('dm')
   if g==1:
       key.append('gas')
   if s==1:
       key.append('stars')
   return key,slide
# create the slice and then find al the particles in the slice, that save the spae by avoinding loading evertything  
def load_particles(path,id,dm=0,g=0,s=0,coordinate=1,extra_entry=[],mode="halo"):
   if mode=="halo":
      keys,slides=load_halo(id,dm,g,s)
   elif mode=="cluster":
      keys,slides=load_cluster(id,dm,g,s)
   else:
      raise ValueError("What on earth do you want to do?")
#   print(keys)
   # create the slice and then find all the particles in the slice, that save the space by avoiding loading everything
   f=h5py.File(path+'particles_ranked.hdf5','r')
   dataset=[]
   if dm==1:
      
      dataset.append(f['PartType1'])
   if g==1:
      
      dataset.append(f['PartType0'])
   if s==1:
       
      dataset.append(f['PartType2'])
   particles=[]

   for i in range(0,len(dataset)):
      comp=[]
      data=dataset[i]
#      if slides[i].start==slides[i].stop:
##         print("no particles in this halo")
#         continue
      
      if coordinate==1:
             
             Coord=np.array(data['Coordinates'][slides[i]],dtype=np.float32)
           
             comp.append(Coord)
             
      if extra_entry[keys[i]]!=[]:
             for entry in extra_entry[keys[i]]:
            
               for slide in slides:
                 
                 entry_data=np.array(data[entry][slides[i]],dtype=np.float32)
                 
                 comp.append(entry_data)
      comp=np.array(comp,dtype=np.float32)#in shape [Coord,entry1, entry2...]
      
      particles.append(comp)#in shape dm, g, s
   f.close()      
      
   return particles
S=np.zeros(len(halo_ids[halo_ids<=0]))
main_id=halo_ids[halo_ids<=0]
for i in tqdm(range(len(S))):
  if np.sum(cross[(halo_ids>i)*(halo_ids<i+1)+(halo_ids==-i)])!=0:
    particle=load_particles(path,main_id[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="cluster")
  
  
    S[i]=fn.dissociation(particle[0][0],particle[1][0])
  else:
    S[i]=-10
print(np.histogram(S[S>=-1],bins=10))
print(len(S[S==-10]),len(cross[cross>0]))
import subprocess
import sys

# Run the first script
#subprocess.run([sys.executable, "script_to_run_first.py", "argument_value"])      
#get the coordinates of the dark matter and gas particles

#analyse the main halo

#calculate the dissociation
'''
S=np.zeros(len(x_dm))exit
for i in  range(0,len(x_dm)):
      S[i]=fn.dissociation(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
      Offset[i]=fn.offset(x_dm[i], y_dm[i], z_dm[i],x_g[i], y_g[i],z_g[i])
'''   


 
#plot
import matplotlib.pyplot as plt
plt.close()
'''
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(S, bins=100)
ax.set_xlabel("S-S_rotate")
ax.set_ylabel("Counts")
ax.set_title("M>10^11")
fig.savefig("/home/jyang/plot/Colibre/L0012N0094/Offset_central_norm.png")


plt.close()
fig = plt.figure()
ax=plt.subplot(1,1,1)
h=ax.hist(Offset, bins=100)
ax.set_xlabel("Offset")
ax.set_ylabel("Counts")
ax.set_title("M>10^11")
fig.savefig("/home/jyang/plot/Colibre/L0012N0094/Offset_norm.png")
'''