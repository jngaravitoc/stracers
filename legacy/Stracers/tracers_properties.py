"""
Script that compute the velocity dispersion and densities of DM particles and
stellar particles using the weights scheme explained in Laporte 13.


author: Nico Garavito-Camargo
university of arizona
May 23/2018.

"""

import numpy as np
import h5py
#import reading_snapshots

from sklearn.neighbors import NearestNeighbors
from astropy.coordinates import SkyCoord
from astropy import units as u



def load_weights(file_name, variable):
        """
        """
        print('loading weights : {}'.format(file_name))
        f = h5py.File(file_name, 'r')
        return f[variable]


def den_profile(r, mass, rbins, rcut):
    """

    Computes the density profiles of a halo.

    Paramters:
    -----------

    r : numpy.array
        Array with the distances to the particles

    mass : numpy.array
        Masses of the particles.
    nbins : int
        Number of radial bins to use.
    rcut : Float
        Radius at which truncate the halo.




    To-do:
    ------

    """



    r_bins = np.linspace(0, rcut, rbins)
    dr = (r_bins[1] - r_bins[0])/2.0
    rho_bins = np.zeros(rbins-1)
    for i in range(1, len(r_bins)):
       index = np.where((r<r_bins[i]) & (r>=r_bins[i-1]))[0]
       V = 4/3. * np.pi * (r_bins[i]**3-r_bins[i-1]**3)
       rho_bins[i-1] = np.sum(mass[index]) / V

    return r_bins+dr, rho_bins


def pos_cartesian_to_galactic(pos, vel):
    """
    Make mock observations in galactic coordinates.
    uses Astropy SkyCoord module.
    l and b coordinates are computed from cartesian coordinates.
    l : [-180, 180]
    b : [-90, 90]

    Parameters:
    -----------
    pos : 3d-numpy array.
    vel : 3d-numpy array.
    lmin : float.
           Minimum latitute of the observation in degrees.
    lmax : float.
           Maximum latitute of the observation in degrees.
    bmin : float.
           Minimum longitude of the observation in degrees.
    bmax : float.
           Maximum longitude of the observation in degrees.
    rmax : float.
           Maximum radius to compute the anisotropy and dispersion
           profiles.
    n_bins_r : int
           Number of bins to do the radial measurements.

    """
    ## transforming to galactic coordinates.

    c_gal = SkyCoord(pos, representation='cartesian',frame='galactic')
    c_gal.representation = 'spherical'

    ## to degrees and range of l.

    l_degrees = c_gal.l.wrap_at(180 * u.deg).radian
    b_degrees = c_gal.b.radian

    ## Selecting the region of observation.


    return l_degrees, b_degrees
    

def vel_cartesian_to_galactic(pos, vel, err_p=0):
    """
    Computes velocities in spherical coordinates from Cartesian.
    Assuming the ISO convention. Where theta is the inclination angle
    and phi the azimuth.

    It also assign a *constant* error to the positions assuming a Gaussian
    distribution.

    Parameters:
    -----------

    pos : numpy.ndarray
        3-D Cartesian positions of the particles.
    vel : numpy.ndarray
        3-D Cartesian velocities of the particles.

    err_p : float
        Error in the distance estimates.

    Returns:
    --------

    vr : numpy.ndarray
        Radial component of the velocity.
    vthehta : numpy.ndarray
        Theta component of the velocity.
    vphi : numpy.ndarray
        Phi component of the velocity.

    """

    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5

    r = np.random.normal(r, r*err_p, size=len(r))

    theta = np.arcsin(pos[:,2]/r)
    phi = np.arctan2(pos[:,1], pos[:,0])

    vr = np.cos(theta)*np.cos(phi)*vel[:,0] \
         +  np.cos(theta)*np.sin(phi)*vel[:,1] \
         +  np.sin(theta)*vel[:,2]

    v_theta = -np.sin(theta)*np.cos(phi)*vel[:,0]\
              - np.sin(theta)*np.sin(phi)*vel[:,1]\
              + np.cos(theta)*vel[:,2]

    v_phi = -np.sin(phi)*vel[:,0] + np.cos(phi)*vel[:,1]

    return vr, v_theta, v_phi



def vel_cartesian_to_spherical(pos, vel, err_p=0):
    """
    Computes velocities in spherical coordinates from Cartesian.
    Assuming the ISO convention. Where theta is the inclination angle
    and phi the azimuth.

    It also assign a *constant* error to the positions assuming a Gaussian
    distribution.

    Parameters:
    -----------

    pos : numpy.ndarray
        3-D Cartesian positions of the particles.
    vel : numpy.ndarray
        3-D Cartesian velocities of the particles.

    err_p : float
        Error in the distance estimates.

    Returns:
    --------

    vr : numpy.ndarray
        Radial component of the velocity.
    vthehta : numpy.ndarray
        Theta component of the velocity.
    vphi : numpy.ndarray
        Phi component of the velocity.

    """

    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5

    r = np.random.normal(r, r*err_p, size=len(r))

    theta = np.arccos(pos[:,2]/r)
    phi = np.arctan2(pos[:,1], pos[:,0])

    vr = np.sin(theta)*np.cos(phi)*vel[:,0] \
         +  np.sin(theta)*np.sin(phi)*vel[:,1] \
         +  np.cos(theta)*vel[:,2]

    v_theta = np.cos(theta)*np.cos(phi)*vel[:,0]\
              + np.cos(theta)*np.sin(phi)*vel[:,1]\
              - np.sin(theta)*vel[:,2]

    v_phi = -np.sin(phi)*vel[:,0] + np.cos(phi)*vel[:,1]

    return vr, v_theta, v_phi

def velocity_dispersion(pos, vel):
    """
    Computes the velocity dispersions in spherical coordinates.


    Parameters:
    ----------
    pos : 3d numpy array
        Array with the cartesian coordinates of the particles
    vel : 3d numpy array
        Array with the cartesian velocities of the particles

    Returns:
    --------
    sigma_r : float
        The value of sigma_r.
    sigma_theta : float
        The value of sigma_theta
    sigma_phi : float
        The value of sigma_phi
    n_part : float
        Number of particles used to compute the velocity dispersion.

    """

    vr, v_theta, v_phi = vel_cartesian_to_spherical(pos, vel)
    sigma_r = np.std(vr)
    sigma_theta = np.std(v_theta)
    sigma_phi = np.std(v_phi)

    n_part = len(vr)

    return sigma_r, sigma_theta, sigma_phi, n_part



def density_profile_octants(pos, vel, nbins, rmax, weights, weighted):

    """
    Computes the density profile in eight octants in the sky,
    defined in galactic coordinates as follows:

    octant 1 :
    octant 2 :
    octant 3 :
    octant 4 :
    octant 5 :
    octant 6 :
    octant 7 :
    octant 8 :

    Parameters:
    -----------

    Output:
    -------



    """

    ## Making the octants cuts:

    d_b_rads = np.linspace(-np.pi/2., np.pi/2., 5)
    d_l_rads = np.linspace(-np.pi, np.pi, 3)
    r_bins = np.linspace(0, 300, 31)

    ## Arrays to store the velocity dispersion profiles
    rho_octants = np.zeros((nbins-1, 8))
    r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
    ## Octants counter, k=0 is for the radial bins!
    k = 0

    l, b = pos_cartesian_to_galactic(pos, vel)

    for i in range(len(d_l_rads)-1):
        for j in range(len(d_b_rads)-1):
            index = np.where((l<d_l_rads[i+1]) & (l>d_l_rads[i]) &\
                             (b>d_b_rads[j]) & (b<d_b_rads[j+1]))

            if weighted==0:
                dr, rho_octants[:,k] \
                = den_profile(pos[index], vel[index], nbins, \
                                         rmax, weights)
            elif weighted==1:
                dr, rho_octants[:,k] \
                = den_profile(pos[index], vel[index], nbins, rmax,\
                                         weights, weighted=1)

            k+=1

    return dr, vr_octants, v_theta_octants, v_phi_octants



def mean_velocities(pos, vel, err_r=0, err_t=0, err_p=0):
    """
    Computes the velocity dispersions for stellar particles using
    the weights fro DM particles.

    Uses Eq: (3) in Laporte 13a to compute the velocity dispersion this is:

    \sigma_* =  \dfrac{\sum_{i}^N \omega_i (v_i - \bar{v_i})}{\sum_i^N \omega_i}

    N = number of particles, \omega_i the weights.

    It also assign errors to the velocities and distances by assigning Gaussian
    errors in both distances and velocities.

    Parameters:
    -----------

    pos : numpy.ndarray 
        Array with the Cartesian coordinates of the particles.

    vel : numpy.ndarray
        Array with the Cartesian velocities of the particles.


    err_r : float  
        Errors in v_r measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_r
        value. By default err_r = 0 means no error.

    err_theta : float
        Errors in v_theta measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_theta
        value. By default err_theta = 0 means no error.

   err_phi : float
        Errors in v_phi measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_phi
        value. By default err_phi = 0 means no error.


    err_p : float
        Errors in v_r measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_r
        value. By default err_r = 0 means no error.

    Returns:
    --------

    mean_r : float
        The value of sigma_r.
    mean_theta : float
        The value of sigma_theta
    mean_phi : float
        The value of sigma_phi

    """

    # error in position. 

    # Compute the velocity in spherical coordinates.
    vr, v_theta, v_phi = vel_cartesian_to_spherical(pos, vel)

    # Arrays for stellar velocity profiles. 
    vr_stellar = np.zeros(len(vr))
    vtheta_stellar = np.zeros(len(vr))
    vphi_stellar = np.zeros(len(vr))


    vr = np.random.normal(vr, err_r, size=len(vr))
    v_theta = np.random.normal(v_theta, err_t/np.sqrt(2), size=len(v_theta))
    v_phi = np.random.normal(v_phi, err_t/np.sqrt(2), size=len(v_phi))

    #vr = vr + np.random.randint(0, 1, size=len(vr))
    #v_theta = v_theta + np.random.randint(-err_t/np.sqrt(2), err_t/np.sqrt(2), size=len(v_theta))
    #v_phi = v_phi + np.random.randint(-err_t/np.sqrt(2), err_t/np.sqrt(2), size=len(v_phi))
    # mean values of the velocities.
    vr_mean = np.mean(vr)
    vtheta_mean = np.mean(v_theta)
    vphi_mean = np.mean(v_phi)

    # Number of particles : 
    n_part = len(vr)

    return vr_mean, vtheta_mean, vphi_mean, n_part



def mean_velocities_weights(pos, vel, weights, err_r=0, err_t=0, err_p=0):
    """
    Computes the mean velocities for stellar particles using
    the weights of the DM particles.

    Uses Eq: (3) in Laporte 13a to compute the velocity dispersion this is:

    \v_j* =  \dfrac{\sum_{i}^N \omega_i (v_{i,j})}{\sum_i^N \omega_i}

    N = number of particles, \omega_i the weights. 
    i = index of particles
    j = velocity component. 

    It also assign errors to the velocities and distances assigning Gaussian
    errors in both distances and velocities.

    Parameters:
    -----------

    pos : numpy.ndarray 
        Array with the Cartesian coordinates of the particles.

    vel : numpy.ndarray
        Array with the Cartesian velocities of the particles.


    err_r : float  
        Errors in v_r measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_r
        value. By default err_r = 0 means no error.

    err_theta : float
        Errors in v_theta measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_theta
        value. By default err_theta = 0 means no error.

    err_phi : float
        Errors in v_phi measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_phi
        value. By default err_phi = 0 means no error.


    err_p : float
        Errors in v_r measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_r
        value. By default err_r = 0 means no error.

    Returns:
    --------

    mean_r : float
        The value of sigma_r.
    mean_theta : float
        The value of sigma_theta
    mean_phi : float
        The value of sigma_phi

    """

    # error in position. 

    # Compute the velocity in spherical coordinates.
    vr, v_theta, v_phi = vel_cartesian_to_galactic(pos, vel)

    # Arrays for stellar velocity profiles. 
    vr_stellar = np.zeros(len(vr))
    vtheta_stellar = np.zeros(len(vr))
    vphi_stellar = np.zeros(len(vr))


    #r = np.random.normal(vr, err_r, size=len(vr))
    #_theta = np.random.normal(v_theta, err_t/np.sqrt(2), size=len(v_theta))
    #_phi = np.random.normal(v_phi, err_t/np.sqrt(2), size=len(v_phi))

    #vr = vr + np.random.randint(0, 1, size=len(vr))
    #v_theta = v_theta + np.random.randint(-err_t/np.sqrt(2), err_t/np.sqrt(2), size=len(v_theta))
    #v_phi = v_phi + np.random.randint(-err_t/np.sqrt(2), err_t/np.sqrt(2), size=len(v_phi))
    # mean values of the velocities.

    vr_mean = np.mean(vr*weights)/np.sum(weights)
    vtheta_mean = np.mean(v_theta*weights)/np.sum(weights)
    vphi_mean = np.mean(v_phi*weights)/np.sum(weights)


    return vr_mean, vtheta_mean, vphi_mean


def velocity_dispersion_weights(pos, vel, weights, err_r=0, err_t=0, err_p=0):
    """
    Computes the velocity dispersions for stellar particles using
    the weights fro DM particles.

    Uses Eq: (3) in Laporte 13a to compute the velocity dispersion this is:

    \sigma_* =  \dfrac{\sum_{i}^N \omega_i (v_i - \bar{v_i})}{\sum_i^N \omega_i}

    N = number of particles, \omega_i the weights.

    It also assign errors to the velocities and distances by assigning Gaussian
    errors in both distances and velocities.

    Parameters:
    -----------

    pos : numpy.ndarray 
        Array with the Cartesian coordinates of the particles.

    vel : numpy.ndarray
        Array with the Cartesian velocities of the particles.

    weights : numpy.ndarray
        Array with the weights of the DM particles.

    err_r : float  
        Errors in v_r measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_r
        value. By default err_r = 0 means no error.

    err_theta : float
        Errors in v_theta measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_theta
        value. By default err_theta = 0 means no error.

   err_phi : float
        Errors in v_phi measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_phi
        value. By default err_phi = 0 means no error.


    err_p : float
        Errors in v_r measurements. This value is assumed as the 
        standard deviation of a Gaussian distribution with mean the true v_r
        value. By default err_r = 0 means no error.

    Returns:
    --------

    sigma_r : float
        The value of sigma_r.
    sigma_theta : float
        The value of sigma_theta
    sigma_phi : float
        The value of sigma_phi

    """

    # error in position. 

    # Compute the velocity in spherical coordinates.
    #vr, v_theta, v_phi = vel_cartesian_to_spherical(pos, vel)
    vr, v_theta, v_phi = vel_cartesian_to_galactic(pos, vel)

    # Arrays for stellar velocity profiles. 
    vr_stellar = np.zeros(len(vr))
    vtheta_stellar = np.zeros(len(vr))
    vphi_stellar = np.zeros(len(vr))


    #vr = np.random.normal(vr, err_r, size=len(vr))
    #v_theta = np.random.normal(v_theta, err_t/np.sqrt(2), size=len(v_theta))
    #v_phi = np.random.normal(v_phi, err_t/np.sqrt(2), size=len(v_phi))

    vr = vr #+ np.random.randint(0, 1, size=len(vr))
    v_theta = v_theta #+ np.random.randint(-err_t/np.sqrt(2), err_t/np.sqrt(2), size=len(v_theta))
    v_phi = v_phi #+ np.random.randint(-err_t/np.sqrt(2), err_t/np.sqrt(2), size=len(v_phi))
    # mean values of the velocities.
    vr_mean = np.mean(vr)
    vtheta_mean = np.mean(v_theta)
    vphi_mean = np.mean(v_phi)

    # Number of particles : 
    n_part = len(vr)

    # stellar velocity dispersions 
    vr_stellar = weights*(vr-vr_mean)**2
    vtheta_stellar = weights*(v_theta-vtheta_mean)**2
    vphi_stellar = weights*(v_phi-vphi_mean)**2

    W_total = np.sum(weights)

    sigma_r = np.sqrt(np.abs(np.sum(vr_stellar))/W_total)
    sigma_theta = np.sqrt(np.sum(vtheta_stellar)/W_total)
    sigma_phi = np.sqrt(np.sum(vphi_stellar)/W_total)
   
    return sigma_r, sigma_theta, sigma_phi, n_part



def velocity_dispersions_r(pos, vel, n_bins, rmin, rmax, weights, weighted, err_r=0, err_t=0, err_p=0):
    """
    Compute the velocity dispersion in radial bins. 

    Parameters:
    ----------
    pos : numpy.ndarray
        Array with the Cartesian coordinates of the particles.
    vel : numpy.ndarray
        Array with the Cartesian velocities of the particles.
    n_bins : int
        Number of radial bins to compute the velocity dispersion.
    rmin : int
        Minimum radius to compute the velocity dispersion.
    rmax : int
        Maximum radius to compute the velocity dispersion.

    Returns:
    --------
    sigma_r : numpy array
    sigma_theta : numpy array
    sigma_phi : numpy array

    to-do:
    ------
    1. Put assert statements for the length of the weights and the positions and
       velocities.
    

    """
    print('rmin', rmin)
    print('rmax', rmax)

    dr = np.linspace(rmin, rmax, n_bins)
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    #r = pos
    vr_disp_r = np.zeros(n_bins-1)
    vtheta_disp_r = np.zeros(n_bins-1)
    vphi_disp_r = np.zeros(n_bins-1)
    n_part = np.zeros(n_bins-1)


    if weighted==1:
        print('Computing the velocity dispersion profile for the stellar halo!')
        for i in range(len(dr)-1):
            index = np.where((r<dr[i+1]) & (r>dr[i]))[0]
            vr_disp_r[i], vtheta_disp_r[i], vphi_disp_r[i], n_part[i]\
             = velocity_dispersion_weights(pos[index], vel[index]\
                                           ,weights[index], err_r, err_t,
                                           err_p)
    elif weighted==0:
        for i in range(len(dr)-1):
            index = np.where((r<dr[i+1]) & (r>dr[i]))[0]
            vr_disp_r[i], vtheta_disp_r[i], vphi_disp_r[i], n_part[i] = velocity_dispersion(pos[index], vel[index])

    return dr, vr_disp_r, vtheta_disp_r, vphi_disp_r, n_part


def anisotropy_parameter(sigma_theta, sigma_phi, sigma_r):
    sigma_t = np.sqrt(sigma_theta**2 + sigma_phi**2)
    beta = 1 - (sigma_t**2 / (2*sigma_r**2))
    return beta



def anisotropy_parameter_weights(pos, vel, n_bins, rmin, rmax, weights, weighted, err_r=0, err_t=0, err_p=0):
    dr, sigma_r, sigma_theta, sigma_phi, n_part = velocity_dispersions_r(pos, vel, n_bins, rmin, rmax, weights, weighted)
    sigma_t = np.sqrt(sigma_theta**2 + sigma_phi**2)
    beta = 1 - (sigma_t**2 / (2*sigma_r**2))
    return dr, beta


def velocity_dispersions_octants(pos, vel, nbins, rmin, rmax, weights, weighted,
                                 err_r=0, err_t=0, err_p=0, **kwargs):

    """
    Computes the velocity dispersion in eight Octants in the sky,
    defined in galactic coordinates as follows:

    octant 1 :
    octant 2 :
    octant 3 :
    octant 4 :
    octant 5 :
    octant 6 :
    octant 7 :
    octant 8 :

    Parameters:
    -----------

    Output:
    -------



    """

    ## Making the octants cuts:

    d_b_rads = np.linspace(-np.pi/2., np.pi/2., 3)
    d_l_rads = np.linspace(-np.pi, np.pi, 5)
    #r_bins = np.linspace(0, 300, 31)

    ## Arrays to store the velocity dispersion profiles
    vr_octants = np.zeros((nbins-1, 8))
    v_theta_octants = np.zeros((nbins-1, 8))
    v_phi_octants = np.zeros((nbins-1, 8))
    n_part_octants = np.zeros((nbins-1, 8))

    if 'beta' in kwargs:
        beta_octants = np.zeros((nbins-1, 8))

    ## Octants counter, k=0 is for the radial bins!
    k = 0

    l, b = pos_cartesian_to_galactic(pos, vel)



    for i in range(len(d_b_rads)-1):
        for j in range(len(d_l_rads)-1):
            index = np.where((l<d_l_rads[j+1]) & (l>d_l_rads[j]) &\
                             (b>d_b_rads[i]) & (b<d_b_rads[i+1]))


            if weighted==0:
                dr, vr_octants[:,k], v_theta_octants[:,k], v_phi_octants[:,k], \
                n_part_octants[:,k] = velocity_dispersions_r(pos[index], vel[index], nbins, \
                                         rmin, rmax, weights[index], 0, err_r, err_t)
            elif weighted==1:
                dr, vr_octants[:,k], v_theta_octants[:,k], v_phi_octants[:,k], \
                n_part_octants[:,k] = velocity_dispersions_r(pos[index], vel[index], nbins, rmin, rmax,\
                                                              weights[index], 1, err_r, err_t)
            if 'beta' in kwargs:
                beta_octants[:,k] = anisotropy_parameter(v_theta_octants[:,k], v_phi_octants[:,k], vr_octants[:,k])
            k+=1


    if 'beta' in kwargs:
        return dr, vr_octants, v_theta_octants, v_phi_octants, n_part_octants, beta_octants
    else:
        return dr, vr_octants, v_theta_octants, v_phi_octants, n_part_octants




def sigma2d_NN(pos, vel, lbins, bbins, n_n, d_slice, weights, err_r=0, err_t=0, err_p=0, relative=False, shell_width=5):
    """
    Returns a 2d histogram of the anisotropy parameter in galactic coordinates.

    Parameters:
    ----------
    pos : numpy ndarray
        3d array with the cartesian positions of the particles.
    vel : numpy.ndarray
        3d array with the cartesian velocoties of the particles.
    lbins : int
        Numer of bins to do the grid in latitude.
    bbins : int
        Number of bins to do the grid in logitude.
    n_n : int
        Number of neighbors.
    d_slice : float
        galactocentric distance to make the slice cut.
    weights : numpy.array
        Array with the weights for the particles
    relative :  If True, the velocity dispersion is computed relative to the
                mean. (default = False)

    Returns:
    --------

    sigma_r_grid : numpy ndarray
        2d array with the radial velocity dispersion.
    sigma_t_grid : numoy ndarray
        2d array with the tangential velocity dispersion.
    """

    ## Defining the grid in galactic coordinates.

    d_b_rads = np.linspace(-np.pi/2., np.pi/2., bbins)
    d_l_rads = np.linspace(-np.pi, np.pi, lbins)
    
    ## Defining the 2d arrays for the velocity dispersion.
    sigma_r_grid = np.zeros((lbins-1, bbins-1))
    sigma_t_grid = np.zeros((lbins-1, bbins-1))
    sigma_theta_grid = np.zeros((lbins-1, bbins-1))
    sigma_phi_grid = np.zeros((lbins-1, bbins-1))


    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5

    # Finding the NN nearest neighbors.
    k = 0
    neigh = NearestNeighbors(n_neighbors=n_n, radius=1, algorithm='ball_tree')
    ngbrs = neigh.fit(pos)

    # Computing mean velocity dispersion.
    if relative == True :
        print('Computing relative with respect to the mean changes in the velocity dispersion')
        index_cut = np.where((r<(d_slice+shell_width/2.)) & (r>(d_slice-shell_width/2.)))
        sigma_r_mean, sigma_theta_mean, sigma_phi_mean, n_p =  velocity_dispersion_weights(pos[index_cut], vel[index_cut], weights[index_cut], err_r, err_t, err_p)



    for i in range(len(d_l_rads)-1):
        for j in range(len(d_b_rads)-1):
            #print(i, j)
            gc = SkyCoord(l=d_l_rads[i]*u.radian, b=d_b_rads[j]*u.radian, frame='galactic', distance=d_slice*u.kpc)
            pos_grid = gc.cartesian.xyz.value

            # Finding the nearest neighbors.
            distances, indices = neigh.kneighbors([pos_grid])
            sigma_r, sigma_theta, sigma_phi, n_part = velocity_dispersion_weights(pos[indices[0,:]], vel[indices[0,:]],  weights[indices[0,:]], err_r, err_t, err_p)

            if relative==True :
                sigma_t_grid[i][j] = ((sigma_theta**2 + sigma_phi**2))**0.5 - (sigma_phi_mean**2 + sigma_theta_mean**2)**0.5
                sigma_r_grid[i][j] = sigma_r - sigma_r_mean
                sigma_theta_grid[i][j] = sigma_theta - sigma_theta_mean
                sigma_phi_grid[i][j] = sigma_phi - sigma_phi_mean

            else :
                sigma_t_grid[i][j] = ((sigma_theta**2 + sigma_phi**2))**0.5
                sigma_r_grid[i][j] = sigma_r
                sigma_theta_grid[i][j] = sigma_theta
                sigma_phi_grid[i][j] = sigma_phi

            k+=1

    return sigma_r_grid, sigma_t_grid, sigma_theta_grid, sigma_phi_grid, n_part


def velocities2d_NN(pos, vel, weights, lbins, bbins, n_n, d_slice,  err_r=0, err_t=0, err_p=0, relative=False, shell_width=5):
    """
    Returns a 2d histogram of the mean velocities in galactic coordinates.

    Parameters:
    ----------
    pos : numpy ndarray
        3d array with the cartesian positions of the particles.
    vel : numpy.ndarray
        3d array with the cartesian velocoties of the particles.
    weights : numpy.array
        1d array with the weights of the stellar particles.
    lbins : int
        Numer of bins to do the grid in latitude.
    bbins : int
        Number of bins to do the grid in logitude.
    n_n : int
        Number of neighbors.
    d_slice : float
        galactocentric distance to make the slice cut.
    weights : numpy.array
        Array with the weights for the particles
    relative :  If True, the velocity dispersion is computed relative to the
                mean. (default = False)

    Returns:
    --------

    vmean_r_grid : numpy ndarray
        2d array with the radial velocity dispersion.
    vmean_t_grid : numpy ndarray
        2d array with the tangential velocity dispersion.
    vmean_theta_grid : numpy ndarray
        2d array with the theta velocity dispersion.
    vmean_phi_grid : numpy ndarray
        2d array with the phi velocity dispersion.

    """

    ## Defining the grid in galactic coordinates.

    d_b_rads = np.linspace(-np.pi/2., np.pi/2., bbins)
    d_l_rads = np.linspace(-np.pi, np.pi, lbins)
    
    ## Defining the 2d arrays for the velocity dispersion.
    meanv_r_grid = np.zeros((lbins-1, bbins-1))
    meanv_t_grid = np.zeros((lbins-1, bbins-1))
    meanv_theta_grid = np.zeros((lbins-1, bbins-1))
    meanv_phi_grid = np.zeros((lbins-1, bbins-1))



    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5

    # Finding the NN nearest neighbors.
    neigh = NearestNeighbors(n_neighbors=n_n, radius=1, algorithm='ball_tree')
    ngbrs = neigh.fit(pos)

    # Computing mean velocity dispersion.
    for i in range(len(d_l_rads)-1):
        for j in range(len(d_b_rads)-1):
            #print(i, j)
            gc = SkyCoord(l=d_l_rads[i]*u.radian, b=d_b_rads[j]*u.radian, frame='galactic', distance=d_slice*u.kpc)
            pos_grid = gc.cartesian.xyz.value

            # Finding the nearest neighbors.
            distances, indices = neigh.kneighbors([pos_grid])
            mean_vr, mean_vtheta, mean_vphi = mean_velocities_weights(pos[indices[0,:]], vel[indices[0,:]], weights[indices], err_r, err_t, err_p)

            meanv_t_grid[i][j] = ((mean_vtheta**2 + mean_vphi**2))**0.5
            meanv_r_grid[i][j] = mean_vr
            meanv_theta_grid[i][j] = mean_vtheta
            meanv_phi_grid[i][j] = mean_vphi
      
    return meanv_r_grid, meanv_t_grid, meanv_theta_grid, meanv_phi_grid



