"""
Saving, loading and carrying weights across snapshots.

Weights are computed on one snapshot (usually the initial conditions) and
then followed to later snapshots through the particle IDs.
"""

import numpy as np

__all__ = ["save_weights", "load_weights", "match_weights"]


def save_weights(file_name, ids, weights, **attrs):
    """
    Save particle IDs and weights to an HDF5 file.

    Extra keyword arguments are stored as file attributes (e.g. the tracer
    profile parameters or the stellar mass).
    """
    import h5py

    ids = np.asarray(ids)
    weights = np.asarray(weights)
    if ids.shape != weights.shape:
        raise ValueError("ids and weights must have the same shape")
    with h5py.File(file_name, "w") as f:
        f.create_dataset("ids", data=ids)
        f.create_dataset("weights", data=weights)
        for key, value in attrs.items():
            f.attrs[key] = value


def load_weights(file_name):
    """Load ``(ids, weights, attrs)`` from a file written by :func:`save_weights`."""
    import h5py

    with h5py.File(file_name, "r") as f:
        return f["ids"][()], f["weights"][()], dict(f.attrs)


def match_weights(ids_ref, weights_ref, ids_new, fill_value=0.0):
    """
    Map weights computed for particles ``ids_ref`` onto the particles
    ``ids_new`` of another snapshot.

    Returns an array aligned with ``ids_new``; particles that are not in
    ``ids_ref`` get ``fill_value``.
    """
    ids_ref = np.asarray(ids_ref)
    weights_ref = np.asarray(weights_ref)
    ids_new = np.asarray(ids_new)
    if len(ids_ref) != len(weights_ref):
        raise ValueError("ids_ref and weights_ref must have the same length")
    if len(np.unique(ids_ref)) != len(ids_ref):
        raise ValueError("ids_ref contains duplicated IDs")

    out = np.full(len(ids_new), fill_value, dtype=float)
    if len(ids_ref) == 0:
        return out
    sorter = np.argsort(ids_ref)
    pos = np.clip(np.searchsorted(ids_ref, ids_new, sorter=sorter), 0, len(ids_ref) - 1)
    found = ids_ref[sorter[pos]] == ids_new
    out[found] = weights_ref[sorter[pos[found]]]
    return out
