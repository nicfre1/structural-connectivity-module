import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra-functions"))

import compute_covariance_gradients as ccg  # noqa: E402
import compute_covariance_gradients_vertex as ccgv  # noqa: E402
import drawem32_to_mcribs as d2m  # noqa: E402
import plot_group_gradients_brain as atlas  # noqa: E402


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
        [sys.executable, str(script_dir / "generate_example_data.py"), "vertex",
         "--out_dir", str(tmp_path), "--n_subjects", "6", "--n_per_hemi", "20"],
        check=True, capture_output=True,
    )
    npz_path = tmp_path / "cohort_thickness_vertex_dens-demo.npz"
    X, mask, n_per_hemi, subjects = ccgv.load_cohort(npz_path)
    assert X.shape == (6, int(mask.sum()))
    assert n_per_hemi == 20
    assert len(subjects) == 6


def test_normalize_name_handles_hemisphere_and_ctx_tokens():
    assert atlas.normalize_name("ctx-lh-bankssts") == ("lh", "bankssts")
    assert atlas.normalize_name("rh_Bankssts") == ("rh", "bankssts")
    assert atlas.normalize_name("bankssts") == (None, "bankssts")
    assert atlas.normalize_name("SFG_L_6_1") == ("lh", "sfg_6_1")
    assert atlas.normalize_name("Sch200_RH_007") == ("rh", "sch200_007")


def test_paint_hemisphere_uses_hemisphere_specific_rows_and_ignores_unknown():
    annot_labels = np.array([0, 0, 1, 2, -1])
    annot_names = ["bankssts", "insula", "unknown"]
    values = {"ctx-lh-bankssts": 1.0, "ctx-rh-bankssts": 9.0, "ctx-lh-insula": 2.0}

    lh, report = atlas.paint_hemisphere(values, annot_labels, annot_names, "lh", {})
    assert np.allclose(lh[:3], [1.0, 1.0, 2.0])
    assert np.isnan(lh[3]) and np.isnan(lh[4])
    assert report["matched"] == ["bankssts", "insula"]
    assert report["unmatched_labels"] == []

    rh, report = atlas.paint_hemisphere(values, annot_labels, annot_names, "rh", {})
    assert np.allclose(rh[:2], [9.0, 9.0])
    assert np.isnan(rh[2])
    assert report["unmatched_labels"] == ["insula"]


def test_paint_hemisphere_paints_hemisphere_free_rows_on_both_sides():
    annot_labels = np.array([0, 1])
    annot_names = ["SFG_L_6_1", "SFG_L_6_2"]
    values = {"SFG_6_1": 0.5, "SFG_6_2": -0.5}
    lh, _ = atlas.paint_hemisphere(values, annot_labels, annot_names, "lh", {})
    rh, _ = atlas.paint_hemisphere(values, annot_labels, annot_names, "rh", {})
    assert np.allclose(lh, [0.5, -0.5]) and np.allclose(rh, [0.5, -0.5])


def test_paint_hemisphere_name_map_reconciles_names_and_reports_unused():
    annot_labels = np.array([0, 1])
    annot_names = ["precentral", "insula"]
    values = {"foo": 0.3, "never_used": 0.9}
    painted, report = atlas.paint_hemisphere(
        values, annot_labels, annot_names, "lh", {"foo": "precentral"}
    )
    assert painted[0] == 0.3 and np.isnan(painted[1])
    assert report["unused_regions"] == ["never_used"]
