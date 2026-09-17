"""
Native-surface resampling chain (vertex_resample_utils.py and the scripts
built on it: build_vertex_grid_assets.py, resample_subject_thickness_vertex.py,
build_cohort_thickness_vertex.py, generate_example_native_surfaces.py).

Kept in its own module (not test_extra_functions.py) because it needs
nibabel -- not in ../requirements.txt, only in
extra-functions/requirements-brain.txt -- and the whole-chain test also
needs Connectome Workbench (`wb_command`, an external tool, not
pip-installable). `pytest.importorskip` at module level skips this entire
file when nibabel is missing (e.g. the standard CI runner); a module-level
skip must not take the unrelated DK-68/vertex-gradient tests in
test_extra_functions.py down with it, hence the split.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra-functions"))

pytest.importorskip("nibabel")

import vertex_resample_utils as vru  # noqa: E402

requires_workbench = pytest.mark.skipif(
    shutil.which("wb_command") is None,
    reason="Connectome Workbench (wb_command) not installed",
)


def test_session_dirs_returns_anat_and_xfm():
    session = Path("/data/sub-001/ses-1")
    anat_dir, xfm_dir = vru.session_dirs(session)
    assert anat_dir == session / "anat"
    assert xfm_dir == session / "xfm"


def test_find_one_raises_on_no_match_and_on_multiple_matches(tmp_path):
    with pytest.raises(FileNotFoundError):
        vru.find_one(tmp_path, "*_thickness.shape.gii")

    (tmp_path / "a_thickness.shape.gii").touch()
    (tmp_path / "b_thickness.shape.gii").touch()
    with pytest.raises(RuntimeError):
        vru.find_one(tmp_path, "*_thickness.shape.gii")


@requires_workbench
def test_native_surface_resampling_chain_runs_end_to_end(tmp_path):
    """Full smoke test with real wb_command and fully synthetic data:
    generate_example_native_surfaces -> build_vertex_grid_assets ->
    build_cohort_thickness_vertex -> compute_covariance_gradients_vertex.
    """
    script_dir = Path(__file__).resolve().parents[1] / "extra-functions"
    wb_command = shutil.which("wb_command")

    def run_script(name, *args):
        subprocess.run(
            [sys.executable, str(script_dir / name), *map(str, args)],
            check=True, capture_output=True,
        )

    native_dir = tmp_path / "native"
    assets_dir = tmp_path / "assets"
    cohort_dir = tmp_path / "cohort"
    results_dir = tmp_path / "results"

    run_script(
        "generate_example_native_surfaces.py",
        "--out_dir", native_dir, "--wb_command", wb_command,
        "--n_subjects", 2, "--n_native", 42, "--seed", 0,
    )
    run_script(
        "build_vertex_grid_assets.py",
        "--template_dir", native_dir / "template", "--wb_command", wb_command,
        "--output_dir", assets_dir, "--n_vertices", 642,
    )
    run_script(
        "build_cohort_thickness_vertex.py",
        "--derivatives_root", native_dir / "derivatives", "--output_dir", cohort_dir,
        "--assets_dir", assets_dir, "--template_dir", native_dir / "template",
        "--wb_command", wb_command, "--lowres", "0.6k", "--fwhm", 5, "--n_subjects", 2,
    )
    cohort_npz = cohort_dir / "cohort_thickness_vertex_dens-0.6k.npz"
    assert cohort_npz.exists()

    run_script(
        "compute_covariance_gradients_vertex.py",
        cohort_npz, "--out_dir", results_dir, "--no_brain_plot",
    )
    assert (results_dir / "group_covariance" / "group_gradients.npz").exists()
