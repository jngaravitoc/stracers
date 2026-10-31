import numpy as np
import matplotlib.pyplot as plt

"""
Units are in Msun and Kpc

"""

def dens_plummer(r, M, a):
    rho = 3*M / (4 *np.pi * a**3) * (1 + r**2/a**2)**(-5/2.)
    return rho


def dens_hernquist(r, M, a):
    rho = M / (2 * np.pi) * a / (r*(r+a)**3)
    print('Computing Hernquist density profile for the tracers with parameters:')
    print('a = {:.3f}'.format(a))
    return rho

def dens_NFWnRvir(r, M, c, Rvir):
    a = Rvir / c
    M = M
    f = np.log(1.0 + Rvir/a) - (Rvir/a / (1.0 + Rvir/a))
    rho = M / ((4.0 * np.pi * a**3.0 * f) * (r / a) * (1.0 + (r/a))**2.0)
    return rho

def dens_Einasto(r, M, n, r_eff):
    """
    Input:
    ------
    r :
    M :
    n :
    r_eff

    Output:
    -------
    rho :

    See Merrit 06: https://arxiv.org/abs/astro-ph/0509417

    """
    assert n>0.5, 'Einasto profile no implemented for n<0.5, see for more details https://arxiv.org/abs/astro-ph/0509417'

    print('Computing Einasto density profile for the tracers with parameters:')
    print('n = {:.3f}, r_eff = {:.3f}'.format(n, r_eff))

    d_n = 3*n - 1/3. + 0.0079/n

    rho = np.exp(-d_n*((r/r_eff)**(1/n)-1))
    return rho
