import numpy as np
import pytest

from stracers import load_weights, match_weights, save_weights


def test_save_load_roundtrip(tmp_path):
    ids = np.arange(10, 20)
    w = np.linspace(0, 1, 10)
    fname = tmp_path / "w.h5"
    save_weights(fname, ids, w, stellar_mass=1e8, profile="hernquist")
    ids2, w2, attrs = load_weights(fname)
    np.testing.assert_array_equal(ids2, ids)
    np.testing.assert_array_equal(w2, w)
    assert attrs["stellar_mass"] == 1e8
    assert attrs["profile"] == "hernquist"


def test_match_weights_shuffled_and_subset():
    rng = np.random.default_rng(0)
    ids_ref = rng.permutation(100) + 1
    w_ref = ids_ref * 0.5
    # later snapshot: shuffled, some particles missing, some new ones
    ids_new = np.concatenate([rng.permutation(ids_ref)[:80], [1000, 1001]])
    w_new = match_weights(ids_ref, w_ref, ids_new)
    np.testing.assert_allclose(w_new[:80], ids_new[:80] * 0.5)
    np.testing.assert_array_equal(w_new[80:], 0.0)


def test_match_weights_rejects_duplicates():
    with pytest.raises(ValueError):
        match_weights([1, 1, 2], [0.1, 0.2, 0.3], [1, 2])
