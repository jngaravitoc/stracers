"""
Mock catalogues of tracer populations.

:func:`write_catalogue` writes the stars of a :class:`stracers.TracerPopulation`
to an HDF5 or CSV file, in chunks, so that catalogues of ~10^8 stars can be
written with the memory of a single chunk. :func:`read_catalogue` reads them
back. Two kinds of catalogue:

- ``kind="stars"``: one row per mock star, drawn as a Poisson number of stars
  on each particle (the stars of a particle share its position and velocity);
- ``kind="particles"``: one row per particle in the radial range, with its
  expected number of stars ``n_stars``.

Positions are relative to the halo centre of the population, and velocities
relative to ``vcenter``.
"""

import datetime
import json
import pathlib
import warnings

import numpy as np

from .kinematics import spherical_velocities

__all__ = ["write_catalogue", "read_catalogue", "catalogue_columns", "COLUMNS"]

#: Columns of a catalogue, in order. ``id`` needs the particle IDs of the
#: population; ``n_stars`` is only in ``kind="particles"`` catalogues.
COLUMNS = ("id", "x", "y", "z", "vx", "vy", "vz", "r", "v_r", "v_theta", "v_phi", "n_stars")

UNITS = {"x": "kpc", "y": "kpc", "z": "kpc", "r": "kpc", "vx": "km/s", "vy": "km/s", "vz": "km/s",
         "v_r": "km/s", "v_theta": "km/s", "v_phi": "km/s", "n_stars": "", "id": ""}

_CSV_WARN_ROWS = 10_000_000


def _resolve_columns(population, kind, columns):
    if kind not in ("stars", "particles"):
        raise ValueError("kind must be 'stars' or 'particles'")
    available = [c for c in COLUMNS
                 if not (c == "id" and population.ids is None)
                 and not (c == "n_stars" and kind == "stars")]
    if columns is None:
        return available
    unknown = [c for c in columns if c not in available]
    if unknown:
        raise ValueError("Unknown or unavailable columns {}; available: {}".format(unknown, available))
    return list(columns)


def _rows(population, idx, columns, vcenter, dtype):
    """The catalogue columns for the particles (or stars) ``idx``."""
    pos = population.pos[idx]
    vel = population.vel[idx]
    if vcenter is not None:
        vel = vel - vcenter
    out = {}
    if "id" in columns:
        out["id"] = population.ids[idx]
    for i, c in enumerate("xyz"):
        if c in columns:
            out[c] = pos[:, i].astype(dtype)
        if "v" + c in columns:
            out["v" + c] = vel[:, i].astype(dtype)
    if "r" in columns:
        out["r"] = np.linalg.norm(pos, axis=1).astype(dtype)
    if {"v_r", "v_theta", "v_phi"} & set(columns):
        for name, v in zip(("v_r", "v_theta", "v_phi"), spherical_velocities(pos, vel)):
            if name in columns:
                out[name] = v.astype(dtype)
    if "n_stars" in columns:
        out["n_stars"] = population.n_stars[idx].astype(np.float64)
    return {c: out[c] for c in columns}


def _chunks(population, kind, rng, chunk_size):
    """Indices (particles, or particles repeated once per star) in chunks of ~chunk_size rows."""
    if kind == "particles":
        sel = np.flatnonzero(population.in_range)
        for start in range(0, len(sel), chunk_size):
            yield sel[start:start + chunk_size]
        return
    sel = np.flatnonzero(population.in_range & (population.n_stars > 0))
    if len(sel) == 0:
        return
    # chunk boundaries from the expected number of stars, so that a chunk has
    # ~chunk_size rows even when particles carry many stars each
    expected = np.cumsum(population.n_stars[sel])
    bounds = np.searchsorted(expected, np.arange(chunk_size, expected[-1], chunk_size), side="right")
    edges = np.unique(np.concatenate([[0], bounds, [len(sel)]]))
    for a, b in zip(edges[:-1], edges[1:]):
        for start in range(a, b, chunk_size):  # also at most chunk_size particles
            part = sel[start:min(b, start + chunk_size)]
            counts = rng.poisson(population.n_stars[part])
            yield np.repeat(part, counts)


def catalogue_columns(population, kind="stars", rng=None, columns=None, vcenter=None, dtype="float32"):
    """
    A mock catalogue in memory, as a dict of column arrays (see
    :func:`write_catalogue` for the arguments). For large catalogues use
    :func:`write_catalogue`, which writes them in chunks.
    """
    population._require_particles()
    columns = _resolve_columns(population, kind, columns)
    rng = np.random.default_rng(rng)
    vcenter = None if vcenter is None else np.asarray(vcenter, dtype=float)
    parts = [_rows(population, idx, columns, vcenter, dtype)
             for idx in _chunks(population, kind, rng, np.iinfo(np.int64).max)]
    if not parts:
        return {c: _rows(population, np.zeros(0, dtype=int), columns, vcenter, dtype)[c] for c in columns}
    return {c: np.concatenate([p[c] for p in parts]) for c in columns}


def _metadata(population, kind, columns, seed, vcenter, n_rows):
    def tolist(x):
        return None if x is None else np.asarray(x, dtype=float).tolist()

    return {
        "tracer": population.tracer, "key": population.key, "reference": population.reference,
        "url": population.url, "r_range": list(population.r_range),
        "fit_range": list(population.fit_range), "kind": kind,
        "expected_stars": population.total, "n_rows": int(n_rows), "seed": seed,
        "snapshot": population.snapshot, "center": tolist(population.center),
        "vcenter": tolist(vcenter), "columns": columns,
        "units": {c: UNITS[c] for c in columns},
        "frame": "galactocentric, relative to the halo centre",
        "stracers_version": _version(),
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }


def _version():
    from . import __version__

    return __version__


def write_catalogue(population, path, kind="stars", rng=None, columns=None, vcenter=None,
                    chunk_size=10_000_000, dtype="float32", compression="gzip", overwrite=False):
    """
    Write a mock catalogue of a tracer population to an HDF5 or CSV file.

    The catalogue is written in chunks of about ``chunk_size`` rows, so its
    size is not limited by memory.

    Parameters
    ----------
    population : TracerPopulation
        The population (:func:`stracers.paint_tracer`, possibly placed on
        another snapshot with :meth:`~stracers.TracerPopulation.at_snapshot`).
    path : str or path
        Output file: ``.h5``/``.hdf5`` (needs h5py) or ``.csv``.
    kind : {"stars", "particles"}
        One row per mock star (Poisson draw), or one row per particle in the
        radial range with its expected ``n_stars``.
    rng : int or numpy.random.Generator, optional
        Random seed or generator for the Poisson draw. With the same integer
        seed the catalogue does not depend on ``chunk_size``.
    columns : list of str, optional
        Columns to write (default: all of :data:`COLUMNS` that apply). With
        all the columns a star takes ~45 bytes in float32.
    vcenter : (3,) array, optional
        Velocity of the halo centre, subtracted from the velocities.
    chunk_size : int
        Rows (and at most particles) per chunk.
    dtype : str
        Float type of the columns (``id`` keeps the type of the IDs, and
        ``n_stars`` is float64).
    compression : str or None
        HDF5 compression filter, e.g. ``"gzip"``, ``"lzf"`` or `None`.
    overwrite : bool
        Replace an existing file.

    Returns
    -------
    int
        Number of rows written.
    """
    population._require_particles()
    path = pathlib.Path(path)
    fmt = {".h5": "hdf5", ".hdf5": "hdf5", ".csv": "csv"}.get(path.suffix.lower())
    if fmt is None:
        raise ValueError("Unknown catalogue format {!r}: use .h5, .hdf5 or .csv".format(path.suffix))
    if path.exists() and not overwrite:
        raise FileExistsError("{} exists; pass overwrite=True to replace it".format(path))
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    columns = _resolve_columns(population, kind, columns)
    seed = int(rng) if isinstance(rng, (int, np.integer)) else None
    rng = np.random.default_rng(rng)
    vcenter = None if vcenter is None else np.asarray(vcenter, dtype=float)
    chunks = (_rows(population, idx, columns, vcenter, dtype)
              for idx in _chunks(population, kind, rng, int(chunk_size)))
    writer = _write_hdf5 if fmt == "hdf5" else _write_csv
    n_rows = writer(path, chunks, population, columns, dtype, compression, chunk_size)
    meta = _metadata(population, kind, columns, seed, vcenter, n_rows)
    _store_metadata(path, fmt, meta)
    return n_rows


def _write_hdf5(path, chunks, population, columns, dtype, compression, chunk_size):
    import h5py

    with h5py.File(path, "w") as f:
        g = f.create_group("catalogue")
        empty = _rows(population, np.zeros(0, dtype=int), columns, None, dtype)
        for c in columns:
            g.create_dataset(c, shape=(0,), maxshape=(None,), dtype=empty[c].dtype,
                             chunks=(int(min(chunk_size, 1_000_000)),), compression=compression)
        n = 0
        for rows in chunks:
            m = len(rows[columns[0]])
            for c in columns:
                g[c].resize((n + m,))
                g[c][n:n + m] = rows[c]
            n += m
    return n


def _write_csv(path, chunks, population, columns, dtype, compression, chunk_size):
    fmts = ["%d" if c == "id" else "%.10g" if c == "n_stars" else "%.7g" for c in columns]
    n = 0
    warned = False
    with open(path, "w") as f:
        f.write("# metadata: {}\n")  # placeholder, rewritten by _store_metadata
        f.write(",".join(columns) + "\n")
        for rows in chunks:
            m = len(rows[columns[0]])
            if m:
                # float64 holds the IDs (< 2^53) exactly, and %d formats them
                np.savetxt(f, np.column_stack([rows[c].astype(np.float64) for c in columns]),
                           fmt=fmts, delimiter=",")
            n += m
            if n > _CSV_WARN_ROWS and not warned:
                warnings.warn("CSV catalogues of more than {:.0e} rows are large and slow; "
                              "HDF5 is recommended".format(_CSV_WARN_ROWS), RuntimeWarning)
                warned = True
    return n


def _store_metadata(path, fmt, meta):
    text = json.dumps(meta)
    if fmt == "hdf5":
        import h5py

        with h5py.File(path, "a") as f:
            f.attrs["metadata"] = text
            f.attrs["format"] = "stracers.catalogue"
        return
    # CSV: replace the first line (the placeholder) with the metadata, without
    # reading the whole file into memory
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(path) as src, open(tmp, "w") as dst:
        src.readline()
        dst.write("# metadata: " + text + "\n")
        for line in src:
            dst.write(line)
    tmp.replace(path)


def read_catalogue(path):
    """
    Read a catalogue written by :func:`write_catalogue`.

    Returns
    -------
    columns : dict of ndarray
        The catalogue columns.
    metadata : dict
        Tracer, fit, ranges, kind, number of rows, seed, snapshot, centre,
        units, version and creation date.
    """
    path = pathlib.Path(path)
    if path.suffix.lower() in (".h5", ".hdf5"):
        import h5py

        with h5py.File(path, "r") as f:
            meta = json.loads(f.attrs["metadata"])
            return {c: f["catalogue"][c][()] for c in meta["columns"]}, meta
    with open(path) as f:
        first = f.readline()
        if not first.startswith("# metadata: "):
            raise ValueError("{} is not a stracers catalogue".format(path))
        meta = json.loads(first[len("# metadata: "):])
        f.readline()  # column names
        data = np.loadtxt(f, delimiter=",", ndmin=2)
    if data.size == 0:
        data = np.zeros((0, len(meta["columns"])))
    cols = {}
    for i, c in enumerate(meta["columns"]):
        cols[c] = data[:, i].astype(np.int64) if c == "id" else data[:, i]
    return cols, meta
