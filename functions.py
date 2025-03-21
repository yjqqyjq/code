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
    x_dmc=np.sum(x_dm)/len(x_dm)
    x_gc=np.sum(x_g)/len(x_g)
    x_dm-=x_dmc
    x_g-=x_gc
    y_dmc=np.sum(y_dm)/len(y_dm)
    y_dm-=y_dmc
    y_gc=np.sum(y_g)/len(y_g)
    y_g-=y_gc
    z_dmc=np.sum(z_dm)/len(z_dm)
    z_dm-=z_dmc
    z_gc=np.sum(z_g)/len(z_g) 
    z_g-=z_gc
    n_dm = len(x_dm)
    n_g = len(x_g)
    r_mean_dm = np.average(radial_distance(x_dm, y_dm,z_dm))
    r_mean_g = np.average(radial_distance(x_g, y_g,z_g))
    r_max = max(r_mean_dm, r_mean_g)
    q_dm = quadrupole(x_dm, y_dm,z_dm, n_dm)
    q_g = quadrupole(x_g, y_g, z_g,n_g)
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
def analyse(sgi,i,x_dm,y_dm,z_dm,x_g,y_g,z_g):#swift galaxy object and index
      x_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates.x)
      y_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates.y)
      z_dm[i]=np.array(sgi.dark_matter.cartesian_coordinates
      .z)
      x_g[i]=np.array(sgi.gas.cartesian_coordinates.x)
      y_g[i]=np.array(sgi.gas.cartesian_coordinates.y)
      z_g[i]=np.array(sgi.gas.cartesian_coordinates.z)