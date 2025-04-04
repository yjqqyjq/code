#math functions
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
def quadrupole(x, y,z,num_particles):
    r = radial_distance(x, y,z)
    f0=np.sum(spherical_harmonic_0(x, y,z)*r)
    f1=np.sum(spherical_harmonic_1(x, y,z)*r)
    f2=np.sum(spherical_harmonic_2(x, y,z)*r)
    f_1=np.sum(spherical_harmonic__1(x, y,z)*r)
    f_2=np.sum(spherical_harmonic__2(x, y,z)*r)
    return np.sqrt(f0**2+f1**2+f2**2+f_1**2+f_2**2) / num_particles
#dissociation functions of an array of coordinates  
def dissociation(x_dm, y_dm,z_dm ,x_g, y_g, z_g):
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
    q_dm = quadrupole(xdm, ydm,zdm, n_dm)
    q_g = quadrupole(xg, yg, zg,n_g)
    return np.sqrt(4*np.pi/5) * (q_dm - q_g) / r_max
#offset functions of an array of coordinates
def offset(x_dm, y_dm,z_dm ,x_g, y_g, z_g):
    x_dmc=np.sum(x_dm)/len(x_dm)
    x_gc=np.sum(x_g)/len(x_g)
    y_dmc=np.sum(y_dm)/len(y_dm)
    y_gc=np.sum(y_g)/len(y_g)
    z_dmc=np.sum(z_dm)/len(z_dm)
    z_gc=np.sum(z_g)/len(z_g)
    return radial_distance(x_dmc-x_gc, y_dmc-y_gc, z_dmc-z_gc)
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
def center_of_mass(x,y,z):#R_m=Mp_dm/Mp_g
    xc=np.average(np.array(x))
    yc=np.average(np.array(y))
    zc=np.average(np.array(z))
 
    return xc,yc,zc