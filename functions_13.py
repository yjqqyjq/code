#math functions
#float 32 ,7 digit, 0.1pc
#Bullet clusters Can; 8
import numpy as np
def radial_distance(x, y,z):
    return np.sqrt(x**2 + y**2+z**2)

def spherical_harmonic_0(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(5/16/np.pi)*(3*z**2-r**2)/r**2
def spherical_harmonic_1(x, y,z):
    r = radial_distance(x, y,z)  
    return np.sqrt(15/4/np.pi)*x*z/r**2
def spherical_harmonic_2(x, y,z): 
    r = radial_distance(x, y,z)
    return np.sqrt(15/16/np.pi)*(x**2-y**2)/r**2
def spherical_harmonic__1(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(15/4/np.pi)*y*z/r**2 
def spherical_harmonic__2(x, y,z):
    r = radial_distance(x, y,z)
    return np.sqrt(15/4/np.pi)*x*y/r**2 
def quadrupole(x, y,z):
    r = radial_distance(x, y,z)
    
    f0=np.sum(spherical_harmonic_0(x, y,z)*r) / len(r)
    f1=np.sum(spherical_harmonic_1(x, y,z)*r) / len(r)
    f2=np.sum(spherical_harmonic_2(x, y,z)*r) / len(r)
    f_1=np.sum(spherical_harmonic__1(x, y,z)*r) / len(r)
    f_2=np.sum(spherical_harmonic__2(x, y,z)*r) / len(r)
    return np.sqrt(f0**2+f1**2+f2**2+f_1**2+f_2**2)
#dissociation functions of an array of coordinates  
def dissociation(Coord_dm ,Coord_g,rdmg):
    x_dm, y_dm, z_dm = Coord_dm[:,0], Coord_dm[:,1], Coord_dm[:,2]
    x_g, y_g, z_g = Coord_g[:,0], Coord_g[:,1], Coord_g[:,2]
    #center of mass
    if len(x_dm)==0 or len(x_g)==0:
        return 0
    else:
       
        xdmc,ydmc,zdmc= center_of_mass(x_dm, y_dm,z_dm)
    
        xgc,ygc,zgc= center_of_mass(x_g, y_g,z_g) 
        xc=(xdmc*rdmg+xgc)/(rdmg+1)
        yc=(ydmc*rdmg+ygc)/(rdmg+1)
        zc=(zdmc*rdmg+zgc)/(rdmg+1)
        xdm=np.array(x_dm-xc)   
        ydm=np.array(y_dm-yc)
        zdm=np.array(z_dm-zc)
        xg=np.array(x_g-xc)
        yg=np.array(y_g-yc)
        zg=np.array(z_g-zc)
        n_dm = len(xdm)
        n_g = len(xg)
    
        r_mean_dm = np.average(radial_distance(xdm, ydm,zdm))
        q_dm = quadrupole(xdm, ydm,zdm)
   
        r_mean_g = np.average(radial_distance(xg, yg,zg))
        q_g = quadrupole(xg, yg, zg)
   
        r_max = max(r_mean_dm, r_mean_g)
   
    
        return np.sqrt(4*np.pi/5) * (q_dm - q_g) / r_max
def dissociation_old(Coord_dm ,Coord_g,rdmg=1):
    x_dm, y_dm, z_dm = Coord_dm[:,0], Coord_dm[:,1], Coord_dm[:,2]
    x_g, y_g, z_g = Coord_g[:,0], Coord_g[:,1], Coord_g[:,2]
    #center of mass
    xdmc,ydmc,zdmc= center_of_mass(x_dm, y_dm,z_dm)
    xgc,ygc,zgc= center_of_mass(x_g, y_g,z_g) 
   
    xdm=np.array(x_dm-xdmc)   
    ydm=np.array(y_dm-ydmc)
    zdm=np.array(z_dm-zdmc)
    xg=np.array(x_g-xgc)
    yg=np.array(y_g-ygc)
    zg=np.array(z_g-zgc)
    n_dm = len(xdm)
    n_g = len(xg)
    r_mean_dm = np.average(radial_distance(xdm, ydm,zdm))
    r_mean_g = np.average(radial_distance(xg, yg,zg))
    r_max = max(r_mean_dm, r_mean_g)
    q_dm = quadrupole(xdm, ydm,zdm)
    q_g = quadrupole(xg, yg, zg)
    return np.sqrt(4*np.pi/5) * (q_dm - q_g) / r_max
#offset functions of an array of coordinates, CoM
def offset(Coord_dm ,Coord_g):
    x_dm, y_dm, z_dm = Coord_dm[:,0], Coord_dm[:,1], Coord_dm[:,2]
    x_g, y_g, z_g = Coord_g[:,0], Coord_g[:,1], Coord_g[:,2]
    x_dmc=np.sum(x_dm)/len(x_dm)
    x_gc=np.sum(x_g)/len(x_g)
    y_dmc=np.sum(y_dm)/len(y_dm)
    y_gc=np.sum(y_g)/len(y_g)
    z_dmc=np.sum(z_dm)/len(z_dm)
    z_gc=np.sum(z_g)/len(z_g)
    return radial_distance(x_dmc-x_gc, y_dmc-y_gc, z_dmc-z_gc)
def offset_lum(x_g, y_g, z_g,rho_g):
    
    x_gc=np.sum(x_g*rho_g)/np.sum(rho_g)
    y_gc=np.sum(y_g*rho_g)/np.sum(rho_g)
    z_gc=np.sum(z_g*rho_g)/np.sum(rho_g)
    return radial_distance(x_gc,y_gc,z_gc)
#analyse the swift galaxy object
def offsetb(x_dm, y_dm,z_dm ,x_g, y_g, z_g,x_stars,y_stars,z_stars,Mg,Ms):
    xdmc,ydmc,zdmc= center_of_mass(x_dm, y_dm,z_dm)
    xbc=(np.sum(x_g)*Mg+np.sum(x_stars)*Ms)/(Mg*len(x_g)+Ms*len(x_stars))
    ybc=(np.sum(y_g)*Mg+np.sum(y_stars)*Ms)/(Mg*len(x_g)+Ms*len(x_stars))
    zbc=(np.sum(z_g)*Mg+np.sum(z_stars)*Ms)/(Mg*len(x_g)+Ms*len(x_stars))
  
    return radial_distance(xdmc-xbc,ydmc-ybc,zdmc-zbc)
def analyse(sgi,i,x_dm,y_dm,z_dm,x_g,y_g,z_g):#swift galaxy object and index
      x_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates.x)
     
      y_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates.y)
      z_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates
      .z)
      x_g[i]=np.array(sgi.gas.cartesian_coordinates.x)
      y_g[i]=np.array(sgi.gas.cartesian_coordinates.y)
      z_g[i]=np.array(sgi.gas.cartesian_coordinates.z)

#Calculate the center of mass of the dark matter and gas particles
def center_of_mass(x,y,z,weight=None):#R_m=Mp_dm/Mp_g
    if weight is None:
 
      xc=np.average(np.array(x))
      yc=np.average(np.array(y))
      zc=np.average(np.array(z))
    else:
      xc=np.sum(np.array(x)*weight)/np.sum(weight)
      yc=np.sum(np.array(y)*weight)/np.sum(weight)
      zc=np.sum(np.array(z)*weight)/np.sum(weight)
 
    return xc,yc,zc


def massfunction(mass):#in 10^10Msun
   
    mass=np.log10(mass)
    bins= np.linspace(0,4.5, 101)
    h= np.histogram(mass, bins=bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    return bin_centers, h[0]  # Normalize the histogram to get the mass function




import unyt
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


import h5py
from tqdm import tqdm


#work with SOAP
#load the data
path="/Users/24756376/data/Flamingo/L1000N1800/"

f=h5py.File(path+'halos_13_ranked.hdf5','r')
N_g=np.array(f['N_g'])

N_dm=np.array(f['N_dm'])
#N_s=np.array(f['N_s'])
N_dm_c=np.array(f['N_dm_c'])
N_g_c=np.array(f['N_g_c'])
#N_s_c=np.array(f['N_s_c'])
N_dm_region=np.array(f['N_dm_region'])
N_g_region=np.array(f['N_g_region'])
#N_s_region=np.array(f['N_s_region'])


#N_s_region_unbound=np.array(f['N_s_region_unbound'])
halo_ids=np.array(f['id'])
mass=np.array(f['mass'])
N_s=np.zeros(len(N_dm))
N_s_c=np.zeros(len(N_dm_c))
N_s_region=np.zeros(len(N_dm_c))
N_s_region_unbound=np.zeros(len(N_dm_c))
#103 halos with m>10**15, 659 for m>5*10**14
m50=np.array(f['m50'])
centers=np.array([f["centers_x"],f["centers_y"],f["centers_z"]]).T
r50=np.array(f["r50"])
#r100=np.array(f["r100"])
#ms100=np.array(f['mass_star_100kpc'])
#ms3000=np.array(f['mass_star_1000kpc'])
r200=np.array(f["r200"])

f.close()

def check_boundry(Coord):#do this only after recentering
   
    x=Coord[:,0]
    y=Coord[:,1]
    z=Coord[:,2]
    if len(x)>0:
      if np.max(x)-np.min(x)>500:
        if np.min(x)<-500:
            Coord[Coord[:,0]<-500,0]+=1000
        elif np.max(x)>500:
            Coord[Coord[:,0]>500,0]-=1000
      if np.max(y)-np.min(y)>500:
        if np.min(y)<-500:
            Coord[Coord[:,1]<-500,1]+=1000
        elif np.max(y)>500:
            Coord[Coord[:,1]>500,1]-=1000
      if np.max(z)-np.min(z)>500:
        if np.min(z)<-500:
            Coord[Coord[:,2]<-500,2]+=1000
        elif np.max(z)>500:
            Coord[Coord[:,2]>500,2]-=1000
   
    return Coord
      
        
        

def load_halo(id,dm=0,g=0,s=0):
      """
      Load halo data from the specified path and ID.
      """
   
      arg=np.nonzero(halo_ids==float(id))[0]
      
    
    
    
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

def load_satellite(id,dm=0,g=0,s=0):
  
   if id>0:
      id=-int(id) 
      print("warning: this is a satellite, loading the cluster it belongs to")
   arg=np.nonzero(halo_ids[halo_ids<=0]==id)[0]
   
   arg=int(arg)
   slide=[]
   if dm==1:
     
   
     dm_s=int(np.sum(N_dm_c[0:arg])+N_dm[halo_ids==id])
     dm_e=int(np.sum(N_dm_c[0:arg+1]))
     slide.append(slice(dm_s,dm_e))
    
   #  print(arg)
   if g==1:
      g_s=int(np.sum(N_g_c[0:arg])+N_g[halo_ids==id])
      g_e=int(np.sum(N_g_c[0:arg+1]))
      slide.append(slice(g_s,g_e))
   if s==1:
      s_s=int(np.sum(N_s_c[0:arg])+N_s[halo_ids==id])
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
def load_particles(path,id,dm=0,g=0,s=0,coordinate=1,extra_entry={},mode="halo"):
   
   if mode=="halo":
      keys,slides=load_halo(id,dm,g,s)
      center=centers[int(np.nonzero(halo_ids==id)[0])]
   elif mode=="cluster":
      keys,slides=load_cluster(id,dm,g,s)
      center=centers[halo_ids<=0][-int(id)]
   elif mode=="satellite":
      keys,slides=load_satellite(id,dm,g,s)
      center=centers[halo_ids<=0][-int(id)]
   else:
      raise ValueError("What on earth do you want to do?")
#   print(keys)
   # create the slice and then find all the particles in the slice, that save the space by avoiding loading everything
   f=h5py.File(path+'cluster_particles_13.hdf5','r')
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

      
      if coordinate==1:
             
            
             Coord=np.array(data['Coordinates'][slides[i]])-center
             Coord=check_boundry(Coord)
             comp.append(Coord)
             
             
      if extra_entry[keys[i]]!=[]:
             for entry in extra_entry[keys[i]]:
               
                 entry_data=np.array(data[entry][slides[i]],dtype=np.float32)
             
                 comp.append(entry_data)
#      comp=np.array(comp,dtype=np.float32)#in shape [Coord,entry1, entry2...]
      
      particles.append(comp)#particle in shape dm, g, s
   f.close()      
      
   return particles



def load_regions_unbound(id,dm,g,s):#id is negative
   arg=int(id)
   
   
   slide=[]
   if dm==1:
     
   
     dm_s=int(np.sum(N_dm_region[0:-arg]))
     dm_e=int(np.sum(N_dm_region[0:-arg+1]))
     slide.append(slice(dm_s,dm_e))
    
   #  print(arg)
   if g==1:
      g_s=int(np.sum(N_g_region[0:-arg]))
      g_e=int(np.sum(N_g_region[0:-arg+1]))
      slide.append(slice(g_s,g_e))
   if s==1:
      s_s=int(np.sum(N_s_region[0:-arg]))
      s_e=int(np.sum(N_s_region[0:-arg+1]))
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
def load_regions(path,id,radius,dm=0,g=0,s=0,coordinate=1,extra_entry={}):
   if id>0:#correct id input to 
      id=-int(id) 
   else:
      id=int(id)
   center=centers[halo_ids<=0][-int(id)]
   

   keys,slides=load_regions_unbound(id,dm,g,s)
   f=h5py.File(path+'particles_region_13.hdf5','r')
      
  
   # create the slice and then find all the particles in the slice, that save the space by avoiding loading everything
   dataset=[]
   
   if dm==1:
      
      dataset.append(f['PartType1'])
   

      
   if g==1:
      
      dataset.append(f['PartType0'])
      

   if s==1:
       
      dataset.append(f['PartType2'])
  
  
    #load entry keys

   particles=[]

   for i in range(0,len(dataset)):
      comp=[]
      data=dataset[i]
      Coord=np.array(data['Coordinates'][slides[i]])-centers[halo_ids<=0][-id]
#      print(Coord,centers[halo_ids<=0][-id]) 
      Coord=check_boundry(Coord)
      r2=Coord[:,0]**2+Coord[:,1]**2+Coord[:,2]**2
      
      Coord=Coord[r2<radius**2]
      
      
      
      
      if coordinate==1:
             
             
             comp.append(Coord)
             
             
      if extra_entry[keys[i]]!=[]:
             for entry in extra_entry[keys[i]]:
                 
                 entry_data=np.array(data[entry][slides[i]],dtype=data[entry][0].dtype)[r2<radius**2]
                 
                 comp.append(entry_data)
#      comp=np.array(comp,dtype=np.float32)#in shape [Coord,entry1, entry2...]
      
      particles.append(comp)#particle in shape dm, g, s
   f.close()      
   
   return particles
