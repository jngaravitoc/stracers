"""Diagnostic plots for :class:`stracers.WeightsResult`."""

import numpy as np

__all__ = ["plot_diagnostics"]


def plot_diagnostics(result, fig=None):
    """
    Plot psi(r), d2nu/dpsi2, the energy distributions and the binned weights.

    Returns the matplotlib figure.
    """
    import matplotlib.pyplot as plt

    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    axes = fig.subplots(2, 2)

    ax = axes[0, 0]
    ax.loglog(result.r_bins, result.psi, c="k")
    ax.set_xlabel(r"$r\,[\mathrm{kpc}]$")
    ax.set_ylabel(r"$\Psi(r)$")

    ax = axes[0, 1]
    ax.plot(result.psi, result.d2nu_dpsi2, c="k")
    ax.set_xlabel(r"$\Psi$")
    ax.set_ylabel(r"$d^2\nu/d\Psi^2$")

    eps = -result.E_centers
    ax = axes[1, 0]
    for y, label in ((result.N_E, "N(E)"), (result.g_E, "g(E)"), (result.f_E, "f(E)")):
        good = y > 0
        if np.any(good):
            ax.semilogy(eps[good], y[good] / y[good].max(), label=label)
    ax.set_xlabel(r"$\epsilon = -E$")
    ax.set_ylabel("normalised")
    ax.legend()

    ax = axes[1, 1]
    good = result.w_bins > 0
    ax.semilogy(eps[good], result.w_bins[good], c="k")
    ax.set_xlabel(r"$\epsilon = -E$")
    ax.set_ylabel(r"$w(E)$")

    fig.tight_layout()
    return fig
