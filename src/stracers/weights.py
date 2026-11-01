r"""
Particle-tagging weights for stellar tracers in (triaxial) N-body halos.

Implements the method of Laporte et al. (2013), following ``weight_triaxial``
in the reference notebook (chervs/tracers_triaxial). Each DM particle gets a
weight that depends only on its energy,

.. math::
    w(E) = \frac{N_*(E)}{N(E)} = \frac{f_*(E)\, g(E)}{N(E)},

where :math:`f_*` is the Eddington distribution function of a tracer with
density :math:`\nu(r)` embedded in the spherically averaged halo potential,
:math:`g(E)` is the density of states and :math:`N(E)` the differential energy
distribution of the DM.

Units: kpc, km/s and Msun; potentials in (km/s)^2 (Gadget convention).

The main entry point is :func:`compute_weights`; the other functions are the
individual steps of the algorithm and are documented for reference.
"""

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.special import betainc, betaln, hyp2f1

from .profiles import mass_fraction_within

__all__ = ["compute_weights", "WeightsResult"]

G_GADGET = 4.30071e-6  # kpc (km/s)^2 / Msun


@dataclass
class WeightsResult:
    """
    Output of :func:`compute_weights`.

    ``weights`` has the same length as the input arrays; particles that were
    not used (outside ``rmax`` or at the centre) have zero weight. The stellar
    mass carried by particle i is ``weights[i] * mass[i]``.

    Attributes
    ----------
    weights : ndarray, shape (N,)
        Per-particle weights, normalised so that ``sum(weights * mass)``
        equals ``stellar_mass``.
    mask : ndarray of bool, shape (N,)
        Particles used in the calculation (0 < r < rmax).
    r_bins : ndarray
        Radii [kpc] at which the spherically averaged potential is evaluated.
    psi : ndarray
        Relative potential psi = -<Phi>(r) at ``r_bins``.
    nu : ndarray
        Tracer density at ``r_bins``.
    dnu_dpsi : ndarray
        First derivative of nu with respect to psi.
    d2nu_dpsi2 : ndarray
        Second derivative of nu with respect to psi.
    E_edges : ndarray
        Edges of the energy bins.
    E_centers : ndarray
        Centres of the energy bins.
    N_E : ndarray
        Differential energy distribution of the DM [Msun / (km/s)^2].
    g_E : ndarray
        Density of states.
    f_E : ndarray
        Tracer distribution function at epsilon = -E_centers (negative values
        clipped to 0).
    w_bins : ndarray
        Weight of each energy bin, with the same normalisation as
        ``weights``.
    n_negative_f : int
        Number of energy bins where the Eddington integral was negative.
    stellar_mass : float
        Stellar mass carried by the particles inside ``rmax``,
        ``sum(weights * mass)`` [Msun].
    r_range : tuple of float
        Radial range [kpc] covered by the radial bins kept, over which the
        density of states is integrated (plus the sphere inside
        ``r_range[0]``). N(E) counts the particles inside ``r_range[1]``,
        except, with ``derivative="gradient"``, those in dropped bins.
    bin_edges : ndarray, shape (n_kept, 2)
        Inner and outer edges [kpc] of the radial bins kept.
    """

    weights: np.ndarray
    mask: np.ndarray
    # radial diagnostics
    r_bins: np.ndarray
    psi: np.ndarray
    nu: np.ndarray
    dnu_dpsi: np.ndarray
    d2nu_dpsi2: np.ndarray
    # energy diagnostics
    E_edges: np.ndarray
    E_centers: np.ndarray
    N_E: np.ndarray
    g_E: np.ndarray
    f_E: np.ndarray
    w_bins: np.ndarray
    n_negative_f: int
    stellar_mass: float
    r_range: tuple
    bin_edges: np.ndarray


def radial_potential(r, pot, rmin, rmax, dlogr, min_per_bin):
    """
    Spherically averaged potential <Phi>(r) in logarithmic radial bins.

    Parameters
    ----------
    r, pot : ndarray
        Radii and potentials of the particles.
    rmin, rmax : float
        Range of the radial bins.
    dlogr : float
        Bin width in log10(r).
    min_per_bin : int
        Bins with fewer particles are discarded.

    Returns
    -------
    r_bins : ndarray
        Bin centres (midpoints in log r).
    dr : ndarray
        Bin widths.
    phi : ndarray
        Mean potential in each bin.
    phi_err : ndarray
        Standard error of the mean potential.
    bin_edges : ndarray, shape (n_kept, 2)
        Inner and outer edges of each bin kept. The density of states must be
        integrated over the same volume as the particles counted in N(E): the
        sphere inside the first edge plus these shells.
    """
    log_min, log_max = np.log10(rmin), np.log10(rmax)
    n_bins = max(int(np.ceil((log_max - log_min) / dlogr)), 1)
    log_edges = np.linspace(log_min, log_max, n_bins + 1)

    idx = np.clip(np.digitize(np.log10(r), log_edges) - 1, 0, n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins)
    pot_sum = np.bincount(idx, weights=pot, minlength=n_bins)
    pot_sq = np.bincount(idx, weights=pot**2, minlength=n_bins)

    keep = counts >= min_per_bin
    if keep.sum() < 3:
        raise ValueError(
            "Fewer than 3 radial bins with >= {} particles; decrease dlogr or "
            "min_per_bin".format(min_per_bin)
        )
    r_bins = 10 ** (0.5 * (log_edges[:-1] + log_edges[1:]))
    dr = 10 ** log_edges[1:] - 10 ** log_edges[:-1]
    n = counts[keep]
    phi = pot_sum[keep] / n
    var = np.clip(pot_sq[keep] / n - phi**2, 0, None)
    phi_err = np.sqrt(var / n)
    bin_edges = np.column_stack([10 ** log_edges[:-1], 10 ** log_edges[1:]])[keep]
    return r_bins[keep], dr[keep], phi, phi_err, bin_edges


def smooth_potential(r_bins, phi, phi_err, n_fine, r_range=None):
    """
    Smooth <Phi>(r) with a spline and resample it on a fine radial grid.

    A quartic smoothing spline is fitted to Phi(ln r) with weights
    1 / ``phi_err`` and smoothing factor equal to the number of bins (i.e. a
    chi^2 of one per bin), so that it follows the potential within its
    statistical errors.

    The fine grid is made of ``n_fine`` logarithmically spaced cells that
    cover ``r_range`` exactly, with the radii at the cell centres (in log r)
    and ``dr`` the cell widths. Passing the range of the ``bin_edges`` of
    :func:`radial_potential`, ``(bin_edges[0, 0], bin_edges[-1, 1])``, keeps
    the radial range of the coarse bins, half a bin beyond the first and last
    bin centres, so that the density of states covers the same volume as the
    particles in N(E). The spline is extrapolated over that half bin and
    interpolates across any bin dropped in between.

    Parameters
    ----------
    r_bins, phi, phi_err : ndarray
        Output of :func:`radial_potential`.
    n_fine : int
        Number of cells in the output grid.
    r_range : tuple of float, optional
        Radial range of the output grid, normally that of the ``bin_edges``
        returned by :func:`radial_potential`. Defaults to the first and last
        bin centres.

    Returns
    -------
    r, dr, phi : ndarray
        Radii, cell widths and smoothed potential on the fine grid.
    """
    lnr = np.log(r_bins)
    # floor the errors so that bins with ~exact potentials do not force interpolation
    err = np.maximum(phi_err, 1e-6 * np.abs(phi).max())
    spl = UnivariateSpline(lnr, phi, w=1.0 / err, k=4, s=len(lnr))
    if r_range is None:
        r_range = (r_bins[0], r_bins[-1])
    ln_edges = np.linspace(np.log(r_range[0]), np.log(r_range[1]), n_fine + 1)
    lnr_fine = 0.5 * (ln_edges[1:] + ln_edges[:-1])
    return np.exp(lnr_fine), np.diff(np.exp(ln_edges)), spl(lnr_fine)


def tracer_derivatives(nu, psi):
    """
    First and second derivatives of nu with respect to psi.

    Computed with second-order finite differences (:func:`numpy.gradient`)
    on the, possibly non-uniform, psi grid.
    """
    dnu_dpsi = np.gradient(nu, psi)
    d2nu_dpsi2 = np.gradient(dnu_dpsi, psi)
    return dnu_dpsi, d2nu_dpsi2


def outer_eddington(epsilon, psi_min, nu_min, dnu_dpsi_min):
    r"""
    Contribution to the Eddington integral from psi < psi_min (r > rmax).

    Computes

    .. math::
        \int_0^{\min(\mathcal{E}, \Psi_{min})} \frac{d^2\nu}{d\Psi^2}
        \frac{d\Psi}{\sqrt{\mathcal{E} - \Psi}}
        + \frac{1}{\sqrt{\mathcal{E}}}
        \left.\frac{d\nu}{d\Psi}\right|_{\Psi=0}

    assuming that beyond the data :math:`\nu` continues as a power law
    :math:`\nu = \nu_{min} (\Psi/\Psi_{min})^k`, with the logarithmic slope
    :math:`k` measured at :math:`\Psi_{min}`. For :math:`k > 1` the boundary
    term vanishes and the integral has the closed form

    .. math::
        \nu_{min} \frac{k (k-1)}{\Psi_{min}^k}\, \mathcal{E}^{k - 3/2}\,
        B\left(x; k-1, \tfrac{1}{2}\right), \qquad
        x = \min(1, \Psi_{min}/\mathcal{E}),

    with :math:`B(x; a, b)` the incomplete beta function. For :math:`k \le 1`
    :math:`\nu` is extrapolated linearly instead, which gives
    :math:`(d\nu/d\Psi)_{min} / \sqrt{\mathcal{E}}`.

    Parameters
    ----------
    epsilon : ndarray
        Relative energies; the result is 0 for epsilon <= 0.
    psi_min, nu_min, dnu_dpsi_min : float
        psi, nu and dnu/dpsi at the outermost radius.
    """
    eps = np.asarray(epsilon, dtype=float)
    out = np.zeros_like(eps)
    pos = eps > 0
    if nu_min <= 0:
        # the profile vanishes at psi_min (e.g. sharply tapered): nothing beyond
        return out
    k = dnu_dpsi_min * psi_min / nu_min
    if k > 1:
        e = eps[pos]
        x = np.minimum(psi_min / e, 1.0)
        ib = betainc(k - 1.0, 0.5, x)
        # evaluated in log space: eps^(k - 3/2) / psi_min^k overflows for steep profiles
        with np.errstate(divide="ignore"):
            log_val = (
                np.log(nu_min * k * (k - 1.0))
                + (k - 1.5) * np.log(e)
                - k * np.log(psi_min)
                + np.log(ib)
                + betaln(k - 1.0, 0.5)
            )
        val = np.exp(log_val)
        # where the incomplete beta underflows (small x) use the equivalent
        # form (nu_min k / psi_min) eps^(-1/2) 2F1(1/2, k-1; k; x), stable there
        small = ib == 0
        val[small] = nu_min * k / psi_min / np.sqrt(e[small]) * hyp2f1(0.5, k - 1.0, k, x[small])
        out[pos] = val
    else:
        out[pos] = dnu_dpsi_min / np.sqrt(eps[pos])
    return out


def eddington(epsilon, psi, d2nu_dpsi2, nu=None, dnu_dpsi=None, outer_slope=None):
    r"""
    Eddington inversion evaluated at the relative energies ``epsilon``.

    .. math::
        f(\mathcal{E}) = \frac{1}{\sqrt{8}\pi^2}
        \int_{\Psi_{min}}^{\mathcal{E}} \frac{d^2\nu}{d\Psi^2}
        \frac{d\Psi}{\sqrt{\mathcal{E} - \Psi}}

    :math:`d^2\nu/d\Psi^2` is taken as piecewise constant (the mean of the two
    nodes) on each segment :math:`[\Psi_j, \Psi_{j+1}]`, and the kernel is
    integrated exactly,

    .. math::
        \int_{\Psi_j}^{\min(\Psi_{j+1}, \mathcal{E})}
        \frac{d\Psi}{\sqrt{\mathcal{E} - \Psi}}
        = 2\left(\sqrt{\mathcal{E} - \Psi_j}
        - \sqrt{\mathcal{E} - \min(\Psi_{j+1}, \mathcal{E})}\right),

    which handles the integrable singularity at :math:`\Psi = \mathcal{E}`
    without noise.

    Parameters
    ----------
    epsilon : ndarray
        Relative energies at which f is evaluated.
    psi, d2nu_dpsi2 : ndarray
        Relative potential and second derivative of nu on the radial grid.
    nu, dnu_dpsi : ndarray, optional
        If both are given, the part of the integral at psi < psi_min is added
        with :func:`outer_eddington`.
    outer_slope : float, optional
        Logarithmic slope dln(nu)/dln(psi) used for that extension. Defaults
        to the local derivative at psi_min, which is the least accurate value
        on the grid; :func:`compute_weights` passes a secant slope over the
        outermost radial bin instead. nu is then continuous at psi_min but
        dnu/dpsi jumps by D = dnu/dpsi(psi_min) - outer_slope nu_min / psi_min,
        and the term D / sqrt(epsilon - psi_min) that this jump contributes to
        the Eddington integral, for epsilon > psi_min, is added.

    Returns
    -------
    f : ndarray
        Distribution function at ``epsilon`` (may be negative).
    """
    order = np.argsort(psi)
    psi = psi[order]
    d2 = d2nu_dpsi2[order]
    d2_seg = 0.5 * (d2[1:] + d2[:-1])
    lo, hi = psi[:-1], psi[1:]

    eps = np.asarray(epsilon, dtype=float)[:, None]
    hi_c = np.minimum(hi[None, :], eps)
    active = lo[None, :] < eps
    kernel = 2.0 * (
        np.sqrt(np.clip(eps - lo[None, :], 0, None))
        - np.sqrt(np.clip(eps - hi_c, 0, None))
    )
    f = np.sum(np.where(active, kernel * d2_seg[None, :], 0.0), axis=1)

    if nu is not None and dnu_dpsi is not None:
        nu_min, dnu_min = nu[order][0], dnu_dpsi[order][0]
        if outer_slope is not None:
            dnu_out = outer_slope * nu_min / psi[0]
            above = eps[:, 0] > psi[0]
            f[above] += (dnu_min - dnu_out) / np.sqrt(eps[above, 0] - psi[0])
            dnu_min = dnu_out
        f += outer_eddington(eps[:, 0], psi[0], nu_min, dnu_min)

    return f / (np.sqrt(8.0) * np.pi**2)


def density_of_states(E, r, dr, phi, r_in=0.0):
    r"""
    Density of states,

    .. math::
        g(E) = (4\pi)^2 \int_0^{r_E} r^2 \sqrt{2 (E - \Phi(r))}\, dr,

    evaluated as a sum over the radial bins with :math:`\Phi(r_k) < E`. The
    sum stops at the outermost bin, consistent with N(E) being measured only
    for particles inside rmax. The sphere :math:`r < r_{in}` inside the first
    bin is added assuming :math:`\Phi(r) = \Phi(r_1)` there, which contributes
    :math:`(4\pi)^2 \frac{r_{in}^3}{3} \sqrt{2 (E - \Phi(r_1))}`.

    Parameters
    ----------
    E : ndarray
        Energies at which g is evaluated.
    r, dr, phi : ndarray
        Radial grid, bin widths and potential.
    r_in : float
        Inner edge of the first radial bin.
    """
    E = np.asarray(E, dtype=float)
    diff = E[:, None] - phi[None, :]
    integrand = r**2 * dr * np.sqrt(2.0 * np.clip(diff, 0, None))
    inner = r_in**3 / 3.0 * np.sqrt(2.0 * np.clip(E - phi[0], 0, None))
    return (4.0 * np.pi) ** 2 * (np.sum(integrand, axis=1) + inner)


def _stellar_mass_inside(
    nu, rmax, stellar_mass, stellar_fraction, within, dm_inside, dm_total, truncated=False
):
    """Stellar mass to be carried by the particles inside rmax."""
    if within not in ("rmax", "total"):
        raise ValueError("within must be 'rmax' or 'total'")
    if stellar_mass is not None and stellar_fraction is not None:
        raise ValueError("Give either stellar_mass or stellar_fraction, not both")
    if stellar_fraction is not None:
        m_star = stellar_fraction * (dm_inside if within == "rmax" else dm_total)
    else:
        m_star = 1.0 if stellar_mass is None else stellar_mass
    if m_star <= 0:
        raise ValueError("The stellar mass must be positive")
    if within == "total" and not truncated:
        m_star *= mass_fraction_within(nu, rmax)
    return m_star


def compute_weights(
    pos,
    vel,
    pot,
    mass,
    tracer_profile,
    *,
    stellar_mass=None,
    stellar_fraction=None,
    within="rmax",
    r_trunc=None,
    rmax=300.0,
    dlogr=0.1,
    min_per_bin=20,
    n_ebins=100,
    center=None,
    derivative="spline",
    n_fine=1000,
    boundary_term=True,
    truncation_correction=False,
    beyond_rmax="zero",
    G=G_GADGET,
):
    """
    Compute stellar-tracer weights for the DM particles of a halo.

    Parameters
    ----------
    pos, vel : (N, 3) arrays
        Positions [kpc] and velocities [km/s] of the DM particles.
    pot : (N,) array
        Potential of each particle [(km/s)^2], e.g. from the snapshot.
    mass : float or (N,) array
        Particle mass(es) [Msun]. The method assumes equal-mass particles.
    tracer_profile : callable
        Tracer density nu(r); its normalisation is irrelevant. See
        :mod:`stracers.profiles`.
    stellar_mass : float, optional
        Stellar mass of the tracer population [Msun]. Defaults to 1 if
        neither ``stellar_mass`` nor ``stellar_fraction`` is given.
    stellar_fraction : float, optional
        Stellar mass as a fraction of the DM mass, M_* = f M_DM. Cannot be
        combined with ``stellar_mass``.
    within : {"rmax", "total"}
        What the stellar mass (or fraction) refers to.

        - "rmax" (default): the stellar mass inside ``rmax``, i.e.
          ``sum(weights * mass) == stellar_mass``. For ``stellar_fraction``
          the DM mass is that of the particles inside ``rmax``.
        - "total": the total stellar mass of the profile, integrated to
          infinity. Only the part inside ``rmax``,
          M_* x :func:`~stracers.profiles.mass_fraction_within`, is carried
          by the particles. For ``stellar_fraction`` the DM mass is that of
          all the particles passed in. Fails for profiles with divergent
          mass (e.g. NFW). With ``r_trunc`` all stars are inside ``rmax``,
          so the total and the mass inside ``rmax`` coincide.
    r_trunc : float, optional
        Outer radius of the stellar halo [kpc], applied as an energy
        truncation: stars with E > <Phi>(r_trunc) get zero weight, so in the
        spherical approximation no star can reach beyond ``r_trunc`` and the
        truncated halo stays in equilibrium. The density then falls to zero
        at ``r_trunc`` and lies below nu(r) close to it. Must be within the
        radial range of the selected particles. For a profile that follows
        nu(r) with a smooth edge instead, taper the profile with
        :func:`~stracers.profiles.truncate`.
    rmax : float
        Only particles with r < rmax are used [kpc].
    dlogr : float
        Width of the logarithmic radial bins used to average the potential.
    min_per_bin : int
        Radial bins with fewer particles are discarded.
    n_ebins : int
        Number of energy bins.
    center : (3,) array, optional
        Halo centre; defaults to the origin.
    derivative : {"spline", "gradient"}
        "spline" (default) fits a smoothing spline to Phi(ln r) and resamples
        it on ``n_fine`` radii before differentiating nu with respect to psi.
        "gradient" differentiates directly on the binned potential, as in the
        reference notebook; it is noisier, especially for cored tracers.
    boundary_term : bool
        Add the part of the Eddington integral at psi < psi(rmax), assuming
        nu(psi) continues as a power law beyond rmax. It is needed to recover
        the tracer profile near rmax and gives non-zero f(E) to particles
        with -E < psi(rmax). The reference notebook omits it.
    truncation_correction : bool
        Subtract G M_tot / max(r) from the potential, as in the reference
        notebook. This is a constant shift of Phi and E; it only changes the
        weights through ``boundary_term``, which assumes Phi(infinity) = 0.
    beyond_rmax : {"zero", "energy"}
        Weights of the particles beyond ``rmax``. "zero" (default): they get
        no weight. "energy": they get the weight of their energy bin, as the
        particles inside ``rmax`` with the same energy (zero if their energy
        is outside the energy grid). The weights, computed from the particles
        inside ``rmax``, are then a function of energy only for every
        particle, so the stellar halo is stationary in a static potential:
        use it to follow the stars in time, when particles beyond ``rmax``
        can orbit back inside. The normalisation (``stellar_mass``) still
        refers to the particles inside ``rmax``.
    G : float
        Gravitational constant in the units of the inputs.

    Returns
    -------
    WeightsResult
    """
    pos = np.asarray(pos, dtype=float)
    vel = np.asarray(vel, dtype=float)
    pot = np.asarray(pot, dtype=float)
    n_part = len(pos)
    mass = np.broadcast_to(np.asarray(mass, dtype=float), (n_part,))

    if center is not None:
        pos = pos - np.asarray(center, dtype=float)
    r_all = np.sqrt(np.sum(pos**2, axis=1))

    if truncation_correction:
        pot = pot - G * np.sum(mass) / np.max(r_all)

    mask = (r_all > 0) & (r_all < rmax)
    if not np.any(mask):
        raise ValueError("No particles with 0 < r < rmax")
    if beyond_rmax not in ("zero", "energy"):
        raise ValueError("beyond_rmax must be 'zero' or 'energy'")
    if r_trunc is not None and r_trunc > rmax:
        raise ValueError("r_trunc must not be larger than rmax")
    m_star = _stellar_mass_inside(
        tracer_profile, rmax, stellar_mass, stellar_fraction, within,
        dm_inside=np.sum(mass[mask]), dm_total=np.sum(mass),
        truncated=r_trunc is not None,
    )
    r = r_all[mask]
    m = mass[mask]
    E = 0.5 * np.sum(vel[mask] ** 2, axis=1) + pot[mask]

    # Spherically averaged potential
    r_bins, dr, phi, phi_err, bin_edges = radial_potential(
        r, pot[mask], r.min(), r.max(), dlogr, min_per_bin
    )
    r_lo, r_hi = bin_edges[0, 0], bin_edges[-1, 1]
    if derivative == "spline":
        r_bins, dr, phi = smooth_potential(r_bins, phi, phi_err, n_fine, r_range=(r_lo, r_hi))
    elif derivative != "gradient":
        raise ValueError("derivative must be 'gradient' or 'spline'")
    psi = -phi
    if np.any(np.diff(psi) >= 0):
        warnings.warn(
            "psi(r) is not strictly decreasing; derivatives of nu with respect "
            "to psi will be noisy. Try a larger dlogr or derivative='spline'.",
            RuntimeWarning,
        )

    # Tracer density and its derivatives
    nu = np.asarray(tracer_profile(r_bins), dtype=float)
    dnu_dpsi, d2nu_dpsi2 = tracer_derivatives(nu, psi)

    # Energy grid shared by N(E), g(E) and f(E)
    E_edges = np.linspace(E.min(), E.max(), n_ebins + 1)
    E_centers = 0.5 * (E_edges[1:] + E_edges[:-1])
    dE = np.diff(E_edges)
    e_idx = np.clip(np.digitize(E, E_edges) - 1, 0, n_ebins - 1)

    # N(E) and g(E) over the same volume: g(E) covers the sphere inside r_lo
    # plus the radial grid. The spline grid is continuous up to r_hi; the
    # gradient grid only has the bins kept, so particles in a dropped
    # intermediate shell are left out too. Particles left out of N(E) still
    # get weights below.
    if derivative == "spline":
        in_volume = r <= r_hi
    else:
        j = np.searchsorted(bin_edges[:, 0], r, side="right") - 1
        in_volume = (r < r_lo) | ((j >= 0) & (r <= bin_edges[np.maximum(j, 0), 1]))
    N_E = np.bincount(e_idx[in_volume], weights=m[in_volume], minlength=n_ebins) / dE
    g_E = density_of_states(E_centers, r_bins, dr, phi, r_in=r_lo)
    if boundary_term:
        # slope of nu(psi) for the extension beyond the grid: a secant over the
        # outermost radial bin (dlogr), robust to the noise of a derivative at
        # the end of the grid
        i_in = min(np.searchsorted(r_bins, r_bins[-1] / 10**dlogr), len(r_bins) - 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            k_out = np.log(nu[-1] / nu[i_in]) / np.log(psi[-1] / psi[i_in])
        f_E = eddington(-E_centers, psi, d2nu_dpsi2, nu, dnu_dpsi,
                        outer_slope=k_out if np.isfinite(k_out) else None)
    else:
        f_E = eddington(-E_centers, psi, d2nu_dpsi2)

    n_negative = int(np.sum(f_E < 0))
    if n_negative:
        warnings.warn(
            "{} energy bins have f(E) < 0; their weights are set to 0. This "
            "is expected in a few poorly sampled, most bound bins; many such "
            "bins mean the tracer profile is inconsistent with an isotropic "
            "DF in this potential.".format(n_negative),
            RuntimeWarning,
        )
    f_E = np.clip(f_E, 0, None)

    unbound_t = np.zeros(len(E), dtype=bool)
    if r_trunc is not None:
        if not r_bins[0] < r_trunc <= r_bins[-1]:
            raise ValueError(
                "r_trunc must lie within the radial range of the particles "
                "({:.3g} to {:.3g} kpc)".format(r_bins[0], r_bins[-1])
            )
        phi_t = -np.interp(np.log(r_trunc), np.log(r_bins), psi)
        f_E[E_centers > phi_t] = 0.0
        # per particle as well, so that the bin containing phi_t is cut exactly
        unbound_t = E > phi_t

    w_bins = np.zeros(n_ebins)
    good = N_E > 0
    w_bins[good] = f_E[good] * g_E[good] / N_E[good]
    if not np.all(np.isfinite(w_bins)):
        raise FloatingPointError(
            "Non-finite weights; the tracer profile may be too steep or noisy "
            "for the radial grid. Try derivative='spline' or a larger r_core."
        )

    w_part = w_bins[e_idx]
    w_part[unbound_t] = 0.0
    total = np.sum(w_part * m)
    if total <= 0:
        raise ValueError("All weights are zero; check the inputs and tracer profile")
    w_part *= m_star / total

    weights = np.zeros(n_part)
    weights[mask] = w_part
    w_bins *= m_star / total
    if beyond_rmax == "energy":
        out = ~mask & (r_all > 0)
        E_out = 0.5 * np.sum(vel[out] ** 2, axis=1) + pot[out]
        j = np.digitize(E_out, E_edges) - 1
        j[E_out == E_edges[-1]] = n_ebins - 1
        valid = (j >= 0) & (j < n_ebins)
        if r_trunc is not None:
            valid &= E_out <= phi_t
        w_out = np.zeros(len(E_out))
        w_out[valid] = w_bins[j[valid]]
        weights[out] = w_out

    return WeightsResult(
        weights=weights,
        mask=mask,
        r_bins=r_bins,
        psi=psi,
        nu=nu,
        dnu_dpsi=dnu_dpsi,
        d2nu_dpsi2=d2nu_dpsi2,
        E_edges=E_edges,
        E_centers=E_centers,
        N_E=N_E,
        g_E=g_E,
        f_E=f_E,
        w_bins=w_bins,
        n_negative_f=n_negative,
        stellar_mass=m_star,
        r_range=(r_lo, r_hi),
        bin_edges=bin_edges,
    )
