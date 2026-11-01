"""
Write a 1M-particle random sample of the dark matter halo of an isolated
Milky Way simulation, used by tests/test_simulation.py and the simulation
tutorial.

The snapshots (Gadget-4 HDF5, 10^8 dark matter particles each) are read in
chunks, so the peak memory stays near 1 GB. The same particles, chosen by ID
from the first snapshot, are kept in every snapshot and stored sorted by ID,
with the halo centre found with nba's shrinking sphere.

Usage::

    python tests/make_sim_sample.py [snapshot_dir] [output_file]

The output defaults to $STRACERS_SIM_SAMPLE, or
~/.cache/stracers/mw_iso_b1_1M.hdf5.
"""

import os
import pathlib
import sys
import time

import h5py
import numpy as np

SIM_DIR = "/home/ngc/Work/research/mwlmc_sims"
SNAP = "MW_100M_beta1_vir_OM3_G4_{}.hdf5"
SNAPSHOTS = ["000", "001", "010", "044", "110", "250"]
N_SAMPLE = 1_000_000
CHUNK = 10_000_000
SEED = 2024
MASS_UNIT = 1e10  # Msun: the snapshots are in kpc, km/s and 1e10 Msun


def default_output():
    return pathlib.Path(os.environ.get(
        "STRACERS_SIM_SAMPLE", pathlib.Path.home() / ".cache" / "stracers" / "mw_iso_b1_1M.hdf5"))


def read_sample(path, ids_sample):
    """Positions, velocities, potential and IDs of the sampled particles, sorted by ID."""
    out = {k: [] for k in ("pos", "vel", "pot", "ids")}
    with h5py.File(path, "r") as f:
        dm = f["PartType1"]
        n = dm["ParticleIDs"].shape[0]
        for i in range(0, n, CHUNK):
            ids = dm["ParticleIDs"][i:i + CHUNK]
            sel = np.isin(ids, ids_sample, assume_unique=True)
            out["ids"].append(ids[sel])
            out["pos"].append(dm["Coordinates"][i:i + CHUNK][sel])
            out["vel"].append(dm["Velocities"][i:i + CHUNK][sel])
            out["pot"].append(dm["Potential"][i:i + CHUNK][sel])
        header = dict(f["Header"].attrs)
        m_particle = float(dm["Masses"][0]) * MASS_UNIT
    out = {k: np.concatenate(v) for k, v in out.items()}
    order = np.argsort(out["ids"])
    out = {k: v[order] for k, v in out.items()}
    return out, float(header["Time"]), m_particle, n


def centre(pos, vel, mass):
    from nba.com import CenterHalo

    halo = {"pos": pos.astype(np.float64), "vel": vel.astype(np.float64),
            "mass": np.full(len(pos), mass)}
    com, vcom = CenterHalo(halo).shrinking_sphere_numba()
    return np.asarray(com, dtype=float), np.asarray(vcom, dtype=float)


def main(sim_dir=SIM_DIR, output=None):
    output = pathlib.Path(output) if output else default_output()
    output.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    with h5py.File(os.path.join(sim_dir, SNAP.format(SNAPSHOTS[0])), "r") as f:
        all_ids = f["PartType1/ParticleIDs"][()]
    ids_sample = np.sort(rng.choice(all_ids, N_SAMPLE, replace=False))
    del all_ids

    with h5py.File(output, "w") as out:
        out.attrs["simulation"] = SNAP.format("XXX")
        out.attrs["n_sample"] = N_SAMPLE
        out.attrs["seed"] = SEED
        out.attrs["units"] = "kpc, km/s, Msun, (km/s)^2"
        for snap in SNAPSHOTS:
            t0 = time.time()
            data, t, m_particle, n_total = read_sample(os.path.join(sim_dir, SNAP.format(snap)), ids_sample)
            assert np.array_equal(data["ids"], ids_sample), "sampled IDs missing in snapshot " + snap
            mass = m_particle * n_total / N_SAMPLE  # the sample keeps the mass of the halo
            com, vcom = centre(data["pos"], data["vel"], mass)
            g = out.create_group(snap)
            for k, v in data.items():
                g.create_dataset(k, data=v)
            g.attrs.update(time=t, mass=mass, center=com, vcenter=vcom, n_total=n_total)
            print(f"{snap}: t = {t:.2f}, centre = {np.round(com, 2)} kpc, "
                  f"vcentre = {np.round(vcom, 2)} km/s ({time.time() - t0:.0f} s)", flush=True)
    print("wrote", output)


if __name__ == "__main__":
    main(*sys.argv[1:])
