import numpy as np
import matplotlib.pyplot as plt
import sys
from pygadgetreader import readsnap
from weights2 import weight_triaxial

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

    rbins, nu_tracer, psi2_hr, dnu_dpsi_smooth, dnu2_dpsi2_smooth = weight_triaxial(rr, Ekk, Epp, ids, partmass, 0.01, 100, 1, 'Hernquist', [40.82])
    print(len(rbins))
    
    import matplotlib.pyplot as plt
    
    
    plt.figure(figsize=(6, 15))
    plt.subplot(3,1,1)
    plt.plot(np.log(psi2_hr), np.log(nu_tracer), c='k')
    #plt.scatter(np.log(rbins), np.log(nu_tracer), s=1, c='k')
    plt.subplot(3,1,2)
    plt.plot(np.log(psi2_hr), np.log(dnu_dpsi_smooth), c='k')
    #plt.scatter(np.log(rbins), np.log(dnu_dpsi), s=1, c='k')
    plt.subplot(3,1,3)
    plt.plot(np.log(psi2_hr), np.log(dnu2_dpsi2_smooth), c='k')
    #plt.scatter(np.log(rbins), np.log(dnu2_dpsi2), s=1, c='k')
    plt.savefig('nu_tracers.png', bbox_inches='tight', dpi=150)
