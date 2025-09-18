import numpy as np
# import matplotlib.pyplot as plt
# import scipy.spatial as spatial
import math
from scipy.special import sph_harm

# some useful numbers
dim = 3 # number of dimensions
N = 100 # number of particles
mass_max = 10 # in 1e10 Msun
size_max = 10 # in Mpc

# first generate random array of coordinates and masses
# these are the basic data points that we have for each particle
coords_gas = np.random.rand(N, dim) * size_max # in Mpc
coords_DM = np.random.rand(N, dim) * size_max # in Mpc
masses_gas = np.random.rand(N) * mass_max # in 1e10 Msun
masses_DM = np.random.rand(N) * mass_max # in 1e10 Msun

# we begin with defining all the functions that we'll use
# the total mass
def calculate_total_mass(masses):
    """
    Calculates the total mass from an array of masses:
        
        M = sum_(i=1)^(N) m_i.
        
     - 'masses' are the individual particle masses (in 1e10 Msun): m_i.
    """
    total_mass = np.sum(masses)
    return total_mass

# the COM coordinates
def calculate_COM(coords, masses, dim):
    """
    Calculates the COM of a set of coordinates with the corresponding masses,
    in a certain number of dimensions:
        
        X_COM = sum_(i=1)^(N) m_i X_i.
    
     - 'coords' are the coordinates of the particles (in Mpc): X_i;
     - 'masses' are the corresponding masses (in 1e10 Msun): m_i;
     - 'dim' is the number of dimensions [usually 3].
    """
    total_mass = calculate_total_mass(masses)
   
    COM = np.zeros(dim)
    for i in range(dim):
        COM[i] = np.sum(coords[:,i] * masses) / total_mass
    
    return COM

# the new origin
def calculate_origin(coords_gas, coords_DM, masses_gas, masses_DM, dim):
    """
    Calculates the new origin of the system of gas and DM particles, based on
    the average of the COM of the gas and DM particles, in a certain number of
    dimensions. Based on Equation (5) of McDonald et al. (2022):
        
        X_0 = 1/2M_DM sum_(i=1)^(N_DM) m_i^DM X_i^DM
            + 1/2M_gas sum_(i=1)^(N_gas) m_i^gas X_i^gas.
    
     - 'coords_gas' are the coordinates of the gas particles (in Mpc): X_i^gas;
     - 'coords_DM' are the coordinates of the DM particles (in Mpc): X_i^DM;
     - 'masses_gas' are the gas particle masses (in 1e10 Msun): m_i^gas;
     - 'masses_DM' are the DM particle masses (in 1e10 Msun): m_i^DM;
     - 'dim' is the number of dimensions [usually 3].
    """
    COM_gas = calculate_COM(coords_gas, masses_gas, dim)
    COM_DM = calculate_COM(coords_DM, masses_DM, dim)
    
    origin = [(COM_gas[i] + COM_DM[i]) / 2 for i in range(dim)]
    
    return origin

# the spherical coordinates
def calculate_sph_coords(coords, origin, dim, N):
    """
    Calculates the spherical coordinates of a set of cartesian coordinates, 
    relative to a certain origin, in a certain number of dimensions, and for
    a certain number of particles.
    
     - 'coords' are the cartesian coordinates of the particles;
     - 'origin' are the three cartesian coordinates of the origin;
     - 'dim' is the number of dimensions [usually 3];
     - 'N' is the number of particles.
    """
    # first we calculate the relative difference to the origin
    diff = [coords[:,i] - origin[i] for i in range(dim)]
    diff = np.array(diff)
    
    x, y, z = diff[0], diff[1], diff[2]
    
    # and then we calculate the spherical coordinates in the usual way
    radii = np.linalg.norm(diff, axis = 0)

    '''
    #This is the orginal code
    thetas = np.zeros(N)
    phis = np.zeros(N)
    for i in range(N):
#        thetas[i] = math.atan2(z[i], math.sqrt(x[i]**2 + y[i]**2))#This is the orginal one
        thetas[i] = math.atan2( math.sqrt(x[i]**2 + y[i]**2),z[i])
        phis[i] = math.atan2(y[i], x[i])
    '''
#    thetas = np.arctan2(np.sqrt(x**2+y**2),z)
    thetas = np.arctan2(z,np.sqrt(x**2+y**2))
    phis = np.arctan2(y, x)
    
    return radii, thetas, phis

# the mean radius
def calculate_mean_radius(mass, masses, radii, N):
    """
    Calculates the mean radius r_ of a component (gas or DM), through Equation
    (7) of McDonald et al. (2022):
        
        r_ = 1/M * sum_(i=1)^(N) r_i m_i.
    
     - 'mass' is the total mass of the component (gas or DM, in 1e10 Msun): M;
     - 'masses' are the individual particle masses (in 1e10 Msun): m_i;
     - 'radii' are the corresponding (spherical coordinate) radii of each
               particle with respect to the newly set origin (in Mpc): r_i;
     - 'N' is the number of particles.
    """    
    '''
    ####This is the orginal code 
    tot = 0
    for i in range(N):
        tot += radii[i] * masses[i]
    
    mean_radius = (1/mass) * tot # in Mpc
    '''
    mean_radius=np.sum(radii*masses)/mass
    return mean_radius

# the quadrupole coefficients
def calculate_quad_coeff(mass, masses, radii, thetas, phis, m, N):
    """
    Calculates the quadrupole coefficients f_m of a component (gas or DM),
    through Equation (6) of McDonald et al. (2022):
        
        f_m = 1/M * sum_(i=1)^(N) r_i m_i Y_2^m(x_i).
    
     - 'mass' is the total mass of the component (gas or DM, in 1e10 Msun): M;
     - 'masses' are the individual particle masses (in 1e10 Msun): m_i;
     - 'radii' are the corresponding (spherical coordinate) radii of each
               particle with respect to the newly set origin (in Mpc): r_i;
     - 'thetas' are the polar coordinates (between 0 and pi);
     - 'phis' are the azimuthal coordinates (between 0 and 2pi);
     - 'm' is the order of the spherical harmonic (m = -2, -1, 0, 1, 2);
     - 'N' is the number of particles.
     
    Note: theta and phi are switched around in the definition of sph_harm
    of the scipy package.
    """
    '''
    ###This is the orginal code
    tot = 0
    for i in range(N):
        tot += radii[i] * masses[i] * sph_harm(m, 2, phis[i], thetas[i])
    
    f_m = (1/mass) * tot # idk what the units are here, not important
    '''
    f_m=np.sum(radii*masses*sph_harm(m,2,phis,thetas))/mass
    return f_m

# the quadrupole amplitude
def calculate_quadrupole(mass, masses, radii, thetas, phis, N):
    """
    Calculates the quadrupole amplitude q from the quadrupole coefficients f_m,
    through Equation (2) of McDonald et al. (2022):
        
        q = (sum_(m=-2)^(2) |f_m|^2)^(1/2).
    
     - 'mass' is the total mass of the component (gas or DM, in 1e10 Msun): M;
     - 'masses' are the individual particle masses (in 1e10 Msun): m_i;
     - 'radii' are the corresponding (spherical coordinate) radii of each
               particle with respect to the newly set origin (in Mpc): r_i;
     - 'thetas' are the polar coordinates (between 0 and pi);
     - 'phis' are the azimuthal coordinates (between 0 and 2pi);
     - 'N' is the number of particles.
    """
    tot = 0
    for m in [-2, -1, 0, 1, 2]:
        f_m = calculate_quad_coeff(mass, masses, radii, thetas, phis, m, N)
        sq = (abs(f_m))**2
      
        tot += sq
    
    q = np.sqrt(tot)
    
    return q

# the dissociation index
def calculate_S(coords_gas, coords_DM, masses_gas, masses_DM, dim = 3):
    """
    Calculates the dissociation index S from the quadrupole amplitudes of gas
    and DM, and the maximum mean radius r_max, through Equation (8) of
    McDonald et al. (2022):
        
        S = sqrt(4pi / 5) * (q_DM - q_gas)/(max(r__DM, r__gas)).
    
    All of these quantities are calculated in their respective functions.
    
     - 'coords_gas' are the coordinates of the gas particles (in Mpc);
     - 'coords_DM' are the coordinates of the DM particles (in Mpc);
     - 'masses_gas' are the individual gas particle masses (in 1e10 Msun);
     - 'masses_DM' are the individual DM particle masses (in 1e10 Msun);
     - 'dim' is the number of dimensions [usually 3].
    """
    # we start from scratch, so by defining the new origin
    origin = calculate_origin(coords_gas, coords_DM, masses_gas, masses_DM, dim)
    
    # a number of particles
    N_gas = len(masses_gas)
    N_DM = len(masses_DM)
    
    # then we need the spherical coordinates
    radii_gas, thetas_gas, phis_gas = calculate_sph_coords(coords_gas, origin, dim, N_gas)
    radii_DM, thetas_DM, phis_DM = calculate_sph_coords(coords_DM, origin, dim, N_DM)
    
    # we can calculate the total masses of the gas and DM
    total_mass_gas = calculate_total_mass(masses_gas)
    total_mass_DM = calculate_total_mass(masses_DM)
    
    # next we can calculate the mean radii of the gas and DM
    mean_radius_gas = calculate_mean_radius(total_mass_gas, masses_gas, radii_gas, N_gas)
    mean_radius_DM = calculate_mean_radius(total_mass_DM, masses_DM, radii_DM, N_DM)
    r_max = max(mean_radius_gas, mean_radius_DM) # the max value of the two
    
    # then we can calculate the quadrupole amplitudes
    q_gas = calculate_quadrupole(total_mass_gas, masses_gas, radii_gas,
                                 thetas_gas, phis_gas, N_gas)
    q_DM = calculate_quadrupole(total_mass_DM, masses_DM, radii_DM,
                                thetas_DM, phis_DM, N_DM)
    
    # and finally from all of this we calculate S
    S = np.sqrt(4*np.pi/5) * (q_DM - q_gas)/r_max
    
    return S

# after all of these definitions, we can calculate the dissociation index from
# the very basic data (coordinates, masses, dimensions, and number of particles)
# S = calculate_S(coords_gas, coords_DM, masses_gas, masses_DM)
# print(S)




































