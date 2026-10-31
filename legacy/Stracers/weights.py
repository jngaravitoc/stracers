import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
from pygadgetreader import readsnap
from .tracers_dens import *
from .smooth import savitzky_golay


"""
To-do:

"""

def rho_tracers(r, M, profile, profile_params):
    """
    Density profiles for the

    Parameters:
    ----------
    r  : numpy.array
        distances to the DM particles.
    M : numpy.array
        masses to the particles.
    profile : str
        density profile of the stellar halo.
    profile_params : list
        parameters of the stellar halo density profile.
        for Hernquist : r_a (scale length)
            Plummer : r_a (scale length)
            NFW : c, rvir (concentration, virial radius)
            Einasto : n, r_eff 

    Returns:
    -------
    rho : numpy.array
        density profile.

    to-do:
    ------
    1. pass as an argument the density function instead of the if statements.
    
    """
    if profile == 'Plummer':
        rho = dens_plummer(r, M, profile_params[0])
    elif profile == 'Hernquist':
        rho = dens_hernquist(r, M, profile_params[0])
    elif profile == 'NFW':
        rho = dens_NFWnRvir(r, M, profile_params[0], profile_params[1])
    elif profile == 'Einasto':
        rho = dens_Einasto(r, M, profile_params[0], profile_params[1])

    return rho



#def interpolate_energies()

def energies(r_part, Ep, v, n_rbins, n_ebins, nbins_coarsed=100):
    r"""

    Parameters:
    -----------

    r_part: 'numpy.array`
        Galactocentric distances of the DM particles in [$kpc$].
    pmass: float
        Mass of a DM particle in units of solar masses.
    Ep : 'numpy.array`
        potential energy in units of
    v : 'numpy.array`
        velocities in km/s
    n_rbins : int
        number of bins, this will be used to bin the potential energy.

    n_ebins : int
        number of bins, this will be used to bin the energy.
    nbins_coarsed : int
        number of bins to bin the potential energy. before interpolation.

    Returns:
    --------

    $\Psi$ : $\Psi = -\Phi + \Phi_0$ Relative potential.
    $\epsilon$ : $\epsilon = -E + \phi_0$ relative total energy.


    """

    G = 4.30071e-6 # Gadget units
    #nbins_coarsed = 20
    E_k = 0.5*v**2 # Kinetic energy / mass.
    E = E_k + Ep
    epsilon = (-1.0)*E

    # Binning the data in logarithmic bins in radius
    # Should I need to worried about centering the bin?
    # yes if you are going to interpolate.

    rbins = np.logspace(np.min(np.log10(r_part)), np.max(np.log10(r_part)), nbins_coarsed)
    # spaced between r_bins
    dr = np.zeros(len(rbins)-1)
    for i in range(len(rbins)-1):
        dr[i] = rbins[i+1] - rbins[i]

    pot = np.zeros(nbins_coarsed-1)

    for i in range(len(rbins)-1):
        index_bins = np.where((r_part<rbins[i+1])
                              & (r_part>=rbins[i]))[0]

        if len(index_bins) == 0:
            pot[i] = 0
            print('Warning : No particles found at r={:0>2f} kpc'.format(rbins[i]))
        else:
            pot[i] = np.mean(Ep[index_bins])


    f_interp_pot = interp1d(rbins[:-1]+dr, pot, kind='cubic')
    r_interp = np.linspace(rbins[0]+dr[0], rbins[-2]+dr[-2], n_rbins)
    pot_interp = f_interp_pot(r_interp)
    psi = (-1.0)*pot_interp



    #Binning Energy for g(E) and f(E) (f(epsilon)) calculations
    Histo_E, Edges = np.histogram(E, bins=n_ebins)
    Histo_epsilon, epsedges = np.histogram(epsilon, bins=n_ebins)

    # Are these always positive?
    dE = Edges[1]-Edges[0]
    depsilon = epsedges[1]-epsedges[0]

    Edges = Edges + dE/2.
    epsedges = epsedges + depsilon/2.


    return r_interp, pot_interp, E, psi, Histo_E, Edges, Histo_epsilon, epsedges

def densities_derivatives(rbins, psi_bins, m_halo, interp_bins, profile, profile_params):
    """
    Computes the derivatives of

    rbins : number of radial bins.

    psi_bins : psi binned.

    inter_bins : int
        Values of the interpolated


    """
    assert len(rbins)==len(psi_bins), "Length of r and psi are not equal and I can't interpolate "

    print('Computing derivatives for the {} profile with parameters {}'.format(profile, profile_params))

    spl1 = InterpolatedUnivariateSpline(rbins, psi_bins)

    # interpolating in the radial bins.
    rbins_hr = np.linspace(min(rbins), max(rbins), interp_bins)
    #nu_tracer_hr = spl1(rbins_hr)
    nu_tracer=rho_tracers(rbins_hr, m_halo, profile, profile_params)/m_halo
    psi_hr = spl1(rbins_hr)

    # First derivative.
    dnu_dpsi = np.gradient(nu_tracer, psi_hr)
    #spl3 = interp1d(rbins, dnu_dpsi, kind='cubic')
    #dnu_dpsi_hr = spl3(rbins_hr)

    # second derivative
    #dnu2_dpsi2 = np.gradient(dnu_dpsi, psi2)

    # smoothing first derivative
    #dnu_dpsi_smooth = savitzky_golay(dnu_dpsi, 5, 3)
    dnu2_dpsi2 = np.gradient(dnu_dpsi, psi_hr)
    # smoothing second derivative
    #dnu2_dpsi2_smooth = savitzky_golay(dnu2_dpsi2, 5, 3)


    index_zeros = np.where(dnu2_dpsi2<0)[0]
    dnu2_dpsi2[index_zeros] = 0 

    plt.figure(figsize=(12, 6))
    plt.subplot(2, 1, 1)
    plt.loglog(rbins_hr, psi_hr, c='k')
    plt.xlabel(r'$r[kpc]$')
    plt.ylabel(r'$\Psi (r)$')
    plt.subplot(2, 1, 2)
    plt.plot(psi_hr, dnu2_dpsi2, c='k')
    #plt.plot(psi_hr[:-20], dnu2_dpsi2_smooth[:-20], c='C9')
    plt.xlabel(r'$\Psi$')
    plt.ylabel(r'$d^2\nu/d \Psi^2$')
    plt.savefig('density_profile_{}_derivaitves.png'.format(profile), bbox_inches='tight')

    # this is to prevent the second derivative to be negative. and avoid making the df negative.

    return rbins_hr, nu_tracer, psi_hr, dnu_dpsi, dnu2_dpsi2


def distribution_function(psi, dnu2_dpsi2, epsilon):
    """
    psi : relative potential

    dnu2_dpsi2 : second derivative of the tracers density with respect to
                 the relative potential.
    epsilon : Energy of the particles.


    return:
    -------

    df : numpy.array
        Distribution function.

    To-do:
    ------

    smooth df?

    """

    assert len(epsilon)<len(psi), 'The length of the energy bins has to be larger than those of the potential'

    factor = 1/(np.sqrt(8)*np.pi**2)
    dpsi = np.zeros(len(psi))
    for i in range(1,len(dpsi)):
        dpsi[i] = np.abs(psi[i]-psi[i-1])
    df = np.zeros(len(epsilon))

    for i in range(len(epsilon)):
        index = np.where(psi<epsilon[i])[0]
        #print(len(index))
        #print(len(index))
        if len(index)==0:
            df[i]=0
        else:
            df[i] = np.sum(dpsi[index]/(np.sqrt(epsilon[i] - psi[index])) * dnu2_dpsi2[index])
            if df[i]<0:
                # Add some check method!
                if np.sum(dpsi[index])<0:
                    print('Its \sum dpsi' )
                if np.sum(dnu2_dpsi2[index])<0:
                    print('Its \sum dnu2_dpsi2')
                if np.min(np.sqrt(epsilon[i] - psi[index]))<0:
                    print('Psi > epsilon')

                assert df[i]>=0, 'df with negative values, something is wrong.'

    df_smooth = savitzky_golay(df, 51, 3)
    return factor*df, factor*df_smooth

def density_of_states(rbins, E, pot):
    """
    Compute the density of states.

    g(E) = (4\pi)^2 \int_0^{r_E} r^2 \sqrt{2(E-\Phi(r))} dr

    Parameters:
    -----------
    rbins : numpy.array
        Array with
    E : numpy.array
        Total Energy.
    pot: numpy.array
        Potential Energy.

    Returns:
    --------

    g_E : numpy.array
        Density of states.

    """

    factor = (4*np.pi)**2
    g_E = np.zeros(len(E))

    dr = np.zeros(len(rbins))
    for i in range(1,len(dr)):
        dr[i] = rbins[i]-rbins[i-1]


    for i in range(len(E)):
        index = np.where(pot<=E[i])[0]
        if len(index)==0:
            g_E[i] = 0
            print('g_w==0 at E={:.2f}'.format(E[i]))
        else:
            r = rbins[index]
            g_E[i] = factor*np.sum(r**2 * np.sqrt(2*dr[index]*(E[i]-pot[index])))
    
            
    return g_E, savitzky_golay(g_E, 9, 3) # smoothing the curve

def differential_energy_distribution(hist_E, E_bins, m_part):
    """
    Differential Energy distribution.

    N(E) = n / dE

    n : number of particles with energy [E, E+dE].
    dE : Energy interval.

    Parameters:
    -----------
    hist_E :

    """
    N_E = np.zeros(len(hist_E))
    for i in range(len(hist_E)):
        dE = np.abs(E_bins[i+1]-E_bins[i])
        N_E[i] = hist_E[i]*m_part / dE

    return N_E, savitzky_golay(N_E, 7, 3) # smoothing the curve

def cast_weights(w, E_part, E_bins):
    """
    Assigns weights to each DM particle.
    For each energy bin it finds all the particles that
    have that energy and give the weight corresponding to that
    energy bin.


    """
    print('Casting particle weights')

    part_weights = np.zeros(len(E_part))
    print('array_weights')
    for i in range(1,len(E_bins)):
        index_part_E = np.where((E_part<E_bins[i]) & (E_part>=E_bins[i-1]))
        part_weights[index_part_E] = w[i]
        print('Done with bin {:1d}'.format(i))

    print('Done casting weights')

    return part_weights

def test_plots(E, NE, gE, df):
    print('Plotting test plot')
    print(min(E), max(E))
    print(min(NE), max(NE))
    print(min(gE), max(gE))
    print(min(df), max(df))

    plt.figure(figsize=(7, 6))
    plt.subplot(3, 1, 1)
    plt.loglog(E, NE, c='k')
    plt.ylabel(r'$N(E)$')
    plt.subplot(3, 1, 2)
    plt.loglog(E, gE, c='k')
    plt.ylabel(r'$g(E)$')
    plt.subplot(3, 1, 3)
    plt.loglog(E, df, c='k')
    plt.ylabel(r'$f(E)$')
    plt.xlabel(r'$-E$')
    plt.savefig('test_figure.png', bbox_inches='tight')
    plt.close()

    #return 0

def weights(r, Epp, v, mp, m_shalo, profiles, profile_params, interp_bins,\
            nr_bins, ne_bins, gsmooth=1, fsmooth=1, nsmooth=1):
    """
    Computes weights:

    r : numpy.array
        Positions 
    Epp: potential energy 
    
    v : 

    mp : float
        Mass of the particles in the simulations. It assumes that all the
        particles have the same mass.

    profiles : str
        Stellar halo density profile. (Hernquist, Plummer, Einasto, NFW)
    interp_bins : int 
        number of bins to do the interpolation for the derivatives of
        the stellar halo density profile.
    nr_bins : int
        number of radial bins.
    ne_bins : int
        number of energy bins.



    """
    n_interp = 10000 
    
    # used to interpolate the final results of g_E, N_E, f and
    #to cast the weights.

    print('Number of particles : ', len(r))
    # Computes energies!
    rbins, pot, E, psi, Histo_E, Edges, Histo_epsilon, eps_edges = energies(r, Epp, v, nr_bins, ne_bins)

    # Computes N_E
    N_E, N_E_smooth = differential_energy_distribution(Histo_E, Edges, mp)


    # Density of states. size = len(Edges)
    g_E, g_E_s = density_of_states(rbins, Edges, pot)

    #  Tracers densities derivatives.
    r_hr,  nu_tracer, psi_hr, dnu_dpsi_smooth, dnu2_dpsi2_smooth = densities_derivatives(rbins,
                                                                                         psi,
                                                                                         m_shalo,
                                                                                         interp_bins=interp_bins,
                                                                                         profile=profiles,
                                                                                         profile_params=profile_params)

    # Distribution function (f size = interp_bins)
    print(len(eps_edges), len(psi_hr))

    f, f_smooth = distribution_function(psi_hr, dnu2_dpsi2_smooth, eps_edges)

    # Interpolating g(E) and N(E)
    
    E_edges_inter = np.linspace(min(Edges), max(Edges[:-1]), n_interp)

    if gsmooth==1:
        print('Smoothing g(E)')
        g_E_interp = interp1d(Edges, g_E_s)
    elif gsmooth==0:
        g_E_interp = interp1d(Edges, g_E)

    g_E_I = g_E_interp(E_edges_inter)

    if nsmooth==1:
        print('Smoothing N(E)')
        N_E_interp = interp1d(Edges[:-1], N_E_smooth)
    elif nsmooth==0:
        N_E_interp = interp1d(Edges[:-1], N_E)

    N_E_I = N_E_interp(E_edges_inter)

    if fsmooth==1:
        print('Smoothing f(E)')
        f_E_interp = interp1d(-Edges, f_smooth)
    elif fsmooth==0:
        f_E_interp = interp1d(-Edges, f)

    f_E_I = f_E_interp(-E_edges_inter)


    test_plots(-E_edges_inter, N_E_I, g_E_I, f_E_I)
    
    # Weights
    w = f_E_I[::-1] * g_E_I / N_E_I

    print('Done computing weights')
    w_p = cast_weights(w, E, E_edges_inter)
    factor_w = np.sum(mp*len(E))
    alpha = factor_w/np.sum(w_p)

    print('Done')
    return w_p*alpha


def weights_snapshot(weights_snap1, ids_snap1, ids_snap2, pos_snap2, vel_snap2\
        ,mass_snap2, **kwargs):

    """
    Re-arrange the weights for a new snapshot (snap2) different from the one used to compute
    the weights (snap1). This is done because for computing the weights you have to
    truncate the halo. The particles (ids, positions, velocities, mass etc..)
    of the snap2 that are also in snap1 are
    selected. 

    Parameters:
    ----------

    weights_snap1 : numpy.array
        Weights computed in snap1.
    weights_ids : numpy.array.
        Ids corresponding to snap1 truncated.
    pos : numpy.array
        Positions of the particles in snap2.
    vel : numpy.array
        Velocities of the particles in snap2.
    mass : numpy.array
        Masses of the particles in snap2.
    kwargs:
        pot
        

    Returns:
    -------
    """

    assert len(ids_snap1) <= len(ids_snap2), 'Error: Length of weights ids larger than length of ids!'
    assert len(weights_snap1) == len(ids_snap1), 'Error: Length of ids different from length of weights.'    
    # Making copies of arrays.
    ## Ids from snap1 that are in snap1
    sorted_ids = np.argsort(ids_snap1)
    ids_sorted_snap1 = ids_snap1[sorted_ids]
    w_sorted = weights_snap1[sorted_ids]


    sorted_ids2 = np.argsort(ids_snap2)
    ids_sorted_snap2 = ids_snap2[sorted_ids2]
    pos_snap2_s = pos_snap2[sorted_ids2]
    vel_snap2_s = vel_snap2[sorted_ids2]
    mass_snap2_s = mass_snap2[sorted_ids2]

    common_ids = np.isin(ids_sorted_snap2, ids_sorted_snap1)


    return w_sorted, pos_snap2_s[common_ids], vel_snap2_s[common_ids], mass_snap2_s[common_ids]






if __name__ == "__main__":

    snapshot = sys.argv[1]


    pp = readsnap(snapshot, 'pos', 'dm')
    vv = readsnap(snapshot, 'vel', 'dm')
    mass = readsnap(snapshot, 'mass', 'dm')
    Epp = readsnap(snapshot, 'pot', 'dm')
    ids = readsnap(snapshot, 'pid', 'dm')


    rr = np.sqrt(pp[:,0]**2+pp[:,1]**2+pp[:,2]**2)

    # truncating the halo

    r_cut = index = np.where(rr<100)[0]

    pp = pp[r_cut]
    rr = rr[r_cut]
    vv = vv[r_cut]
    mass = mass[r_cut]
    Epp = Epp[r_cut]
    ids = ids[r_cut]

    # Energies

    partmass=mass[3]*1e10 #generated the halo particles as "bulge"-type in Gadget file
    v2=vv[:,0]**2+vv[:,1]**2+vv[:,2]**2
    Ekk=0.5*v2
    #weight, p_ids = weight_triaxial(rr,Ekk,Epp,ids,partmass,0.01,100,1, profile)

    rbins, nu_tracer, psi2, dnu_dpsi, dnu2_dpsi2, rbins_hr, nu_tracer_hr, psi2_hr, dnu_dpsi_hr, dnu2_dpsi2_hr = weight_triaxial(rr, Ekk, Epp, ids, partmass, 0.01, 100, 1, 'Hernquist', [40.82])
    print(len(rbins))


    plt.figure(figsize=(6, 15))
    plt.subplot(3,1,1)
    plt.plot(np.log(rbins_hr), np.log(nu_tracer_hr), lw=1, c='C9')
    plt.scatter(np.log(rbins), np.log(nu_tracer), s=1, c='k')
    plt.subplot(3,1,2)
    plt.plot(np.log(rbins_hr), np.log(dnu_dpsi_hr), lw=1, c='C9')
    plt.scatter(np.log(rbins), np.log(dnu_dpsi), s=1, c='k')
    plt.subplot(3,1,3)
    plt.plot(np.log(rbins_hr), np.log(dnu2_dpsi2_hr), lw=1, c='C9')
    plt.scatter(np.log(rbins), np.log(dnu2_dpsi2), s=1, c='k')
    plt.savefig('nu_tracers.png', bbox_inches='tight', dpi=150)
