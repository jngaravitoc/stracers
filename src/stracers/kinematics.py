"""
Kinematics of the stellar tracers.

The stars painted onto dark matter particles move with them, so the
kinematics of the stellar halo are weighted averages over the particles, with
the weights of :func:`stracers.compute_weights` (or the stellar masses
``weights * mass``; the normalisation does not matter):

    mean:        v_bar = sum_i w_i v_i / sum_i w_i
    dispersion:  sigma^2 = sum_i w_i (v_i - v_bar)^2 / sum_i w_i

Spherical components follow the physics (ISO) convention: theta is the polar
angle from the z axis and phi the azimuth, so v_theta is positive towards
increasing theta (away from +z). Units are those of the inputs (kpc, km/s).
"""

from dataclasses import dataclass

import numpy as np

__all__ = [
    "spherical_velocities",
    "weighted_mean",
    "weighted_dispersion",
    "bulk_velocity",
    "anisotropy",
    "velocity_profile",
    "VelocityProfile",
]


def spherical_velocities(pos, vel):
    """
    Radial, polar and azimuthal velocity components.

    Parameters
    ----------
    pos, vel : (N, 3) arrays
        Positions and velocities relative to the centre of the halo.

    Returns
    -------
    v_r, v_theta, v_phi : ndarray
        Velocity components along r, theta (polar angle) and phi (azimuth).
        Particles exactly at the centre get v_r = v_theta = v_phi = 0.
    """
    pos = np.asarray(pos, dtype=float)
    vel = np.asarray(vel, dtype=float)
    x, y, z = pos.T
    vx, vy, vz = vel.T
    R = np.hypot(x, y)
    r = np.hypot(R, z)
    with np.errstate(invalid="ignore", divide="ignore"):
        v_r = np.where(r > 0, (x * vx + y * vy + z * vz) / r, 0.0)
        v_theta = np.where((r > 0) & (R > 0), (z * (x * vx + y * vy) / R - R * vz) / r, 0.0)
        v_phi = np.where(R > 0, (x * vy - y * vx) / R, 0.0)
    return v_r, v_theta, v_phi


def weighted_mean(values, weights=None, axis=0):
    """
    Weighted mean, ``sum(w * x) / sum(w)`` along ``axis``.

    With ``weights=None`` all the values have the same weight. Returns `nan`
    if the weights sum to zero.
    """
    values = np.asarray(values, dtype=float)
    if weights is None:
        return np.mean(values, axis=axis)
    w = np.asarray(weights, dtype=float)
    if values.ndim > w.ndim:  # e.g. (N, 3) velocities with (N,) weights
        w = np.expand_dims(w, tuple(range(w.ndim, values.ndim)))
    total = np.sum(w * np.ones_like(values), axis=axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sum(w * values, axis=axis) / total


def weighted_dispersion(values, weights=None, axis=0):
    """
    Weighted dispersion, ``sqrt(sum(w * (x - x_bar)^2) / sum(w))``, with
    ``x_bar`` the weighted mean along ``axis``.
    """
    values = np.asarray(values, dtype=float)
    mean = weighted_mean(values, weights, axis=axis)
    return np.sqrt(weighted_mean((values - np.expand_dims(mean, axis)) ** 2, weights, axis=axis))


def bulk_velocity(vel, weights=None):
    """Weighted mean velocity vector of the tracers, shape (3,)."""
    return weighted_mean(np.asarray(vel, dtype=float), weights, axis=0)


def anisotropy(sigma_r, sigma_theta, sigma_phi):
    r"""
    Velocity anisotropy,
    :math:`\beta = 1 - (\sigma_\theta^2 + \sigma_\phi^2) / (2 \sigma_r^2)`.

    :math:`\beta = 0` for an isotropic system, 1 for purely radial and
    :math:`-\infty` for purely tangential orbits.
    """
    sigma_r = np.asarray(sigma_r, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - (np.asarray(sigma_theta) ** 2 + np.asarray(sigma_phi) ** 2) / (2.0 * sigma_r**2)


@dataclass
class VelocityProfile:
    """
    Output of :func:`velocity_profile`; all arrays have one value per radial
    bin, and bins without tracers have `nan`.

    Attributes
    ----------
    r : ndarray
        Bin centres (geometric mean of the edges) [kpc].
    edges : ndarray
        Bin edges [kpc].
    n : ndarray of int
        Number of particles with non-zero weight in each bin.
    n_eff : ndarray
        Effective number of particles, ``sum(w)**2 / sum(w**2)``: the
        statistical uncertainty of a dispersion is about
        ``sigma / sqrt(2 n_eff)``.
    mass : ndarray
        Sum of the weights in each bin (the stellar mass if the weights are
        stellar masses).
    mean_vr, mean_vtheta, mean_vphi : ndarray
        Weighted mean velocities [km/s]; ``mean_vphi`` measures rotation.
    sigma_r, sigma_theta, sigma_phi : ndarray
        Weighted velocity dispersions [km/s].
    sigma_3d : ndarray
        Total dispersion, ``sqrt(sigma_r**2 + sigma_theta**2 + sigma_phi**2)``.
    beta : ndarray
        Velocity anisotropy (:func:`anisotropy`).
    """

    r: np.ndarray
    edges: np.ndarray
    n: np.ndarray
    n_eff: np.ndarray
    mass: np.ndarray
    mean_vr: np.ndarray
    mean_vtheta: np.ndarray
    mean_vphi: np.ndarray
    sigma_r: np.ndarray
    sigma_theta: np.ndarray
    sigma_phi: np.ndarray
    sigma_3d: np.ndarray
    beta: np.ndarray


def velocity_profile(pos, vel, weights=None, bins=20, rmin=None, rmax=None,
                     center=None, vcenter=None):
    """
    Radial profiles of the mean velocities and velocity dispersions of the
    tracers.

    Parameters
    ----------
    pos, vel : (N, 3) arrays
        Positions [kpc] and velocities [km/s] of the dark matter particles.
    weights : (N,) array, optional
        Weights of the particles, e.g. ``res.weights`` from
        :func:`stracers.compute_weights` or the stellar masses
        ``res.weights * mass``. Without weights, the profiles are those of the
        dark matter.
    bins : int or array_like
        Radial bin edges [kpc], or a number of logarithmic bins between
        ``rmin`` and ``rmax``.
    rmin, rmax : float, optional
        Radial range for an integer ``bins``; by default that of the
        particles with non-zero weight.
    center : (3,) array, optional
        Position of the halo centre; defaults to the origin.
    vcenter : (3,) array, optional
        Velocity of the halo centre, subtracted from the velocities; defaults
        to zero. Use e.g. the centre-of-mass velocity of the halo (not the
        bulk velocity of the tracers, which would remove any net motion of the
        stellar halo relative to it).

    Returns
    -------
    VelocityProfile
    """
    pos = np.asarray(pos, dtype=float)
    vel = np.asarray(vel, dtype=float)
    if center is not None:
        pos = pos - np.asarray(center, dtype=float)
    if vcenter is not None:
        vel = vel - np.asarray(vcenter, dtype=float)
    w = np.ones(len(pos)) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != (len(pos),):
        raise ValueError("weights must have one value per particle")
    if np.any(w < 0):
        raise ValueError("weights must be non-negative")

    r = np.linalg.norm(pos, axis=1)
    use = (w > 0) & (r > 0)
    if np.ndim(bins) == 0:
        lo = r[use].min() if rmin is None else rmin
        hi = r[use].max() if rmax is None else rmax
        edges = np.geomspace(lo, hi, int(bins) + 1)
    else:
        edges = np.asarray(bins, dtype=float)
    n_bins = len(edges) - 1

    idx = np.digitize(r, edges) - 1
    idx[r == edges[-1]] = n_bins - 1  # include the outer edge
    inside = use & (idx >= 0) & (idx < n_bins)
    idx, w = idx[inside], w[inside]
    v_r, v_theta, v_phi = spherical_velocities(pos[inside], vel[inside])

    def binned(x):
        return np.bincount(idx, weights=x, minlength=n_bins)

    mass = binned(w)
    with np.errstate(invalid="ignore", divide="ignore"):
        n_eff = mass**2 / binned(w**2)
        means, sigmas = [], []
        for v in (v_r, v_theta, v_phi):
            mean = binned(w * v) / mass
            # about the mean of each bin (two passes, for numerical accuracy)
            var = binned(w * (v - mean[idx]) ** 2) / mass
            means.append(mean)
            sigmas.append(np.sqrt(var))

    sigma_r, sigma_theta, sigma_phi = sigmas
    return VelocityProfile(
        r=np.sqrt(edges[1:] * edges[:-1]),
        edges=edges,
        n=np.bincount(idx, minlength=n_bins),
        n_eff=n_eff,
        mass=mass,
        mean_vr=means[0],
        mean_vtheta=means[1],
        mean_vphi=means[2],
        sigma_r=sigma_r,
        sigma_theta=sigma_theta,
        sigma_phi=sigma_phi,
        sigma_3d=np.sqrt(sigma_r**2 + sigma_theta**2 + sigma_phi**2),
        beta=anisotropy(sigma_r, sigma_theta, sigma_phi),
    )
