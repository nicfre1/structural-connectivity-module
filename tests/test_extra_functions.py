import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra-functions"))

import compute_covariance_gradients as ccg  # noqa: E402
import compute_covariance_gradients_vertex as ccgv  # noqa: E402
import drawem32_to_mcribs as d2m  # noqa: E402


def test_region_basename_strips_prefix_and_suffix():
    assert d2m.region_basename("lh_dwm_07_thickness") == "dwm_07"
    assert d2m.region_basename("rh_bankssts_thickness") == "bankssts"
    assert d2m.region_basename("lh_dwm_01_area") == "dwm_01"


def test_convert_hemisphere_mean_aggregation():
    df = pd.DataFrame({
        "sample": ["s1", "s2"],
        "lh_a_thickness": [1.0, 2.0],
        "lh_b_thickness": [3.0, 4.0],
    })
    mapping = {"regions": {"target": ["a", "b"]}}
    out, report = d2m.convert_hemisphere(
        df, ["lh_a_thickness", "lh_b_thickness"], "sample", "lh", mapping, "mean"
    )
    assert list(out["lh_target_thickness"]) == [2.0, 3.0]
    assert (report["status"] == "OK").all()


def test_convert_hemisphere_flags_missing_sources():
    df = pd.DataFrame({"sample": ["s1"], "lh_a_thickness": [1.0]})
    mapping = {"regions": {"target": ["a"], "orphan": ["missing"]}}
    out, report = d2m.convert_hemisphere(
        df, ["lh_a_thickness"], "sample", "lh", mapping, "mean"
    )
    assert "lh_orphan_thickness" not in out.columns
    orphan_row = report.loc[report["mcribs_region"] == "orphan"].iloc[0]
    assert orphan_row["status"] == "NO_SOURCE"


def test_individualized_covariance_from_z_is_symmetric_with_unit_diagonal():
    z = np.array([0.5, -1.2, 0.0, 2.1])
    matrix = ccg.individualized_covariance_from_z(z)
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 1.0)
    assert (matrix > 0).all() and (matrix <= 1).all()


def test_rescale_to_minus1_1_bounds():
    x = np.array([1.0, 5.0, 3.0, -2.0])
    rescaled = ccg.rescale_to_minus1_1(x)
    assert np.isclose(rescaled.min(), -1.0)
    assert np.isclose(rescaled.max(), 1.0)


def test_scatter_to_full_places_values_at_mask_and_nan_elsewhere():
    values = np.array([[1.0, 2.0], [3.0, 4.0]])
    mask = np.array([True, False, True, False])
    full = ccgv.scatter_to_full(values, mask)
    assert full.shape == (4, 2)
    assert np.array_equal(full[mask], values)
    assert np.isnan(full[~mask]).all()


def test_load_cohort_roundtrips_generated_example(tmp_path):
    import subprocess
    import sys

    script_dir = Path(__file__).resolve().parents[1] / "extra-functions"
    subprocess.run(
        [sys.executable, str(script_dir / "generate_example_thickness_vertex.py"),
         "--out_dir", str(tmp_path), "--n_subjects", "6", "--n_per_hemi", "20"],
        check=True, capture_output=True,
    )
    npz_path = tmp_path / "cohort_thickness_vertex_dens-demo.npz"
    X, mask, n_per_hemi, subjects = ccgv.load_cohort(npz_path)
    assert X.shape == (6, int(mask.sum()))
    assert n_per_hemi == 20
    assert len(subjects) == 6
