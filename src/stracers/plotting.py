"""
Plots: diagnostics of :class:`stracers.WeightsResult` and the Milky Way
stellar halo profiles of :mod:`stracers.mw`.
"""

import numpy as np

__all__ = ["plot_diagnostics", "plot_mw_profiles", "LATEX_STYLE"]

#: Matplotlib settings for LaTeX-like figures (Computer Modern, all text in
#: math mode) without a LaTeX installation. Use them with
#: ``matplotlib.pyplot.rc_context(LATEX_STYLE)`` around the plotting *and*
#: saving of a figure. Where LaTeX is installed, ``text.usetex = True`` with
#: ``text.latex.preamble = r"\usepackage{amsmath}"`` renders the labels of
#: :func:`plot_mw_profiles` with LaTeX itself.
LATEX_STYLE = {
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["cmr10", "DejaVu Serif"],
    "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False,
}


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


def _math_text(text):
    """
    Text as upright LaTeX math, e.g. 'Bell et al. 2008' ->
    $\\mathrm{Bell\\ et\\ al.\\ 2008}$. "&" becomes "and" (matplotlib's
    mathtext has no \\& and a bare & is not valid in LaTeX math), and hyphens
    are set as text so that they do not become minus signs.
    """
    text = text.replace(", &", " and").replace(" & ", " and ")
    text = text.replace(" ", r"\ ").replace("-", r"\text{-}")
    return r"$\mathrm{" + text + "}$"


def plot_mw_profiles(keys=None, tracers=None, ax=None, r_norm=30.0, r_extrap=200.0,
                     r_plot_min=1.0, legend=True):
    """
    Plot Milky Way stellar halo density profiles from :mod:`stracers.mw`.

    Each profile is drawn in its published form (no core, sharp breaks),
    normalised to 1 at ``r_norm``: a solid line over the radial range of the
    fit, and a dashed line extrapolating it from the end of that range out to
    ``r_extrap``. Fits quoted only with an outer radius ("<40 kpc") are drawn
    solid from ``r_plot_min``. All the text is in LaTeX math mode; use
    :data:`LATEX_STYLE` for Computer Modern fonts.

    Parameters
    ----------
    keys : list of str, optional
        Profiles to plot; defaults to all of :func:`stracers.mw.available`.
    tracers : list of str, optional
        Only plot the profiles traced by these populations, as in the
        ``tracer`` column of the table (e.g. ``["RRLs", "RRab"]``).
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is made if not given.
    r_norm : float
        Radius [kpc] where every profile is normalised to 1.
    r_extrap : float
        Radius [kpc] up to which the profiles are extrapolated.
    r_plot_min : float
        Inner radius [kpc] of the plot.
    legend : bool
        Add a legend with the references.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    from . import mw

    if keys is None:
        keys = mw.available()
    entries = [mw.get_entry(k) for k in keys]
    if tracers is not None:
        entries = [e for e in entries if e["tracer"] in tracers]
    if ax is None:
        _, ax = plt.subplots(figsize=(8.0, 5.0))

    refs = [(e["reference"], e["model"]) for e in entries]
    colors = plt.get_cmap("tab20" if len(entries) > 10 else "tab10")
    r_right = r_extrap
    for i, e in enumerate(entries):
        nu = mw.mw_profile(e["key"], r_core=0.0, smoothness=0.0)
        norm = float(nu(r_norm))
        lo = e["r_min"] if e["r_min"] is not None else r_plot_min
        hi = e["r_max"]
        r_right = max(r_right, hi)
        label = "{} ({})".format(e["reference"], e["model"].upper())
        if refs.count((e["reference"], e["model"])) > 1:
            label = "{} ({}, {})".format(e["reference"], e["model"].upper(), e["tracer"])
        color = colors(i % colors.N)
        r = np.geomspace(lo, hi, 200)
        ax.plot(r, nu(r) / norm, color=color, lw=1.6, label=_math_text(label))
        if hi < r_extrap:
            r = np.geomspace(hi, r_extrap, 100)
            ax.plot(r, nu(r) / norm, color=color, lw=1.2, ls="--")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(r_plot_min, 1.1 * r_right)
    ax.set_xlabel(r"$r\ [\mathrm{kpc}]$")
    ax.set_ylabel(r"$\nu(r)\, /\, \nu(%g\ \mathrm{kpc})$" % r_norm)
    ax.axvline(r_norm, color="0.6", lw=0.8, ls=":")
    if legend:
        ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    return ax
