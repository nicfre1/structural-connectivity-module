#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Shared resampling helpers: native vertex-wise cortical thickness -> a common
low-resolution surface grid, via Connectome Workbench (`wb_command`).

Used by build_vertex_grid_assets.py, resample_subject_thickness_vertex.py and
build_cohort_thickness_vertex.py, so that a single subject and a whole cohort
go through EXACTLY the same resampling chain. See
others/vertex_resample_utils.md for details.

Resampling chain, per hemisphere:

    native thickness (subject's own native mesh, e.g. dHCP ~88k vertices)
      |   wb_command -metric-resample  (ADAP_BARY_AREA, area-corrected)
      |   sphere: the subject's own native-to-common-space registration
      v
    thickness on the common template mesh (inter-subject correspondence)
      |   wb_command -metric-resample  (BARYCENTRIC)
      |   sphere: common template mesh -> low-resolution icosphere (build_vertex_grid_assets.py)
      v
    thickness on the low-resolution grid  (what gets stacked for the gradients)

No study data, no gradient equations here -- geometry only. This is a
generic port of a resampling scheme validated on dHCP neonatal data; it
depends only on filenames/roles (native thickness, native medial-wall mask,
native midthickness, native-to-template registration sphere, template
sphere/area) that follow common structural-pipeline (BIDS-derivatives-style)
conventions, not on any particular study or cohort.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

HEMISPHERES = (("left", "L"), ("right", "R"))


# ---------------------------------------------------------------------------
# small file utilities
# ---------------------------------------------------------------------------
def find_one(directory: Path, pattern: str) -> Path:
    matches = sorted(
        p for p in directory.glob(pattern) if not p.name.startswith("._")
    )
    if not matches:
        raise FileNotFoundError(f"No {pattern!r} in {directory}")
    if len(matches) > 1:
        raise RuntimeError(f"More than one {pattern!r}: {matches}")
    return matches[0]


def first_array(path: Path) -> np.ndarray:
    image = nib.load(str(path))
    if not image.darrays:
        raise RuntimeError(f"No data array in {path}")
    return np.asarray(image.darrays[0].data).squeeze()


def save_shape(values: np.ndarray, path: Path, name: str = "data") -> Path:
    array = nib.gifti.GiftiDataArray(
        np.asarray(values, dtype=np.float32),
        intent="NIFTI_INTENT_SHAPE",
        meta={"Name": name},
    )
    nib.save(nib.gifti.GiftiImage(darrays=[array]), str(path))
    return path


def run(*command: object) -> None:
    printable = " ".join(str(item) for item in command)
    print(f"  + {printable}", flush=True)
    subprocess.run([str(item) for item in command], check=True)


# ---------------------------------------------------------------------------
# core: one subject, one hemisphere -> thickness on the low-resolution grid
# ---------------------------------------------------------------------------
def resample_hemisphere(
    *,
    hemisphere: str,
    anat_dir: Path,
    xfm_dir: Path,
    template_dir: Path,
    assets_dir: Path,
    work_dir: Path,
    wb_command: str,
    lowres_label: str,
    thickness_desc: str,
) -> tuple[Path, Path]:
    """Returns (thickness .shape.gii, roi .shape.gii) on the low-res grid."""
    short = dict(HEMISPHERES)[hemisphere]
    work_dir.mkdir(parents=True, exist_ok=True)

    # --- native inputs --------------------------------------------------
    if thickness_desc in ("", "thickness"):
        thickness_path = find_one(anat_dir, f"*_hemi-{hemisphere}_thickness.shape.gii")
    else:
        thickness_path = find_one(
            anat_dir, f"*_hemi-{hemisphere}_{thickness_desc}.shape.gii"
        )
    medialwall_path = find_one(
        anat_dir, f"*_hemi-{hemisphere}_desc-medialwall_mask.shape.gii"
    )
    native_midthickness = find_one(
        anat_dir, f"*_hemi-{hemisphere}_midthickness.surf.gii"
    )
    native_registered_sphere = find_one(
        xfm_dir,
        f"*_hemi-{hemisphere}_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii",
    )

    thickness = first_array(thickness_path).astype(np.float32)
    medialwall = first_array(medialwall_path).astype(np.float32)  # 1 = cortex
    if thickness.shape != medialwall.shape:
        raise RuntimeError(
            f"hemi-{hemisphere}: thickness {thickness.shape} vs mask {medialwall.shape}"
        )

    native_thick = save_shape(thickness, work_dir / f"native_thick_{short}.shape.gii", "thick")
    native_roi = save_shape(medialwall, work_dir / f"native_roi_{short}.shape.gii", "roi")

    native_area = work_dir / f"native_va_{short}.shape.gii"
    run(wb_command, "-surface-vertex-areas", native_midthickness, native_area)

    # --- template references ---------------------------------------------
    dhcp32k_sphere = (
        template_dir
        / f"week-40_hemi-{hemisphere}_space-dhcpSym_dens-32k_sphere.surf.gii"
    )
    dhcp32k_area = assets_dir / f"area_hemi-{short}_dhcpSym_dens-32k.shape.gii"
    lowres_sphere = assets_dir / f"sphere_hemi-{short}_dhcpSym_dens-{lowres_label}.surf.gii"
    for required in (dhcp32k_sphere, dhcp32k_area, lowres_sphere):
        if not required.exists():
            raise FileNotFoundError(
                f"{required} missing -- run build_vertex_grid_assets.py first"
            )

    # --- step A: native -> common template mesh (area-corrected) ---------
    thick_32k = work_dir / f"thick_dhcpSym32k_{short}.shape.gii"
    roi_32k = work_dir / f"roi_dhcpSym32k_{short}.shape.gii"
    run(
        wb_command, "-metric-resample",
        native_thick, native_registered_sphere, dhcp32k_sphere, "ADAP_BARY_AREA",
        thick_32k,
        "-area-metrics", native_area, dhcp32k_area,
        "-current-roi", native_roi,
        "-valid-roi-out", roi_32k,
    )

    # --- step B: common template mesh -> low-resolution grid -------------
    thick_low = work_dir / f"thick_dhcpSym{lowres_label}_{short}.shape.gii"
    roi_low = work_dir / f"roi_dhcpSym{lowres_label}_{short}.shape.gii"
    run(
        wb_command, "-metric-resample",
        thick_32k, dhcp32k_sphere, lowres_sphere, "BARYCENTRIC",
        thick_low,
        "-current-roi", roi_32k,
        "-valid-roi-out", roi_low,
    )
    return thick_low, roi_low


def resample_subject(
    *,
    anat_dir: Path,
    xfm_dir: Path,
    template_dir: Path,
    assets_dir: Path,
    work_dir: Path,
    wb_command: str,
    lowres_label: str,
    thickness_desc: str = "thickness",
) -> dict[str, tuple[Path, Path]]:
    """{'L': (thick, roi), 'R': (thick, roi)} on the low-resolution grid."""
    out: dict[str, tuple[Path, Path]] = {}
    for hemisphere, short in HEMISPHERES:
        out[short] = resample_hemisphere(
            hemisphere=hemisphere,
            anat_dir=anat_dir,
            xfm_dir=xfm_dir,
            template_dir=template_dir,
            assets_dir=assets_dir,
            work_dir=work_dir,
            wb_command=wb_command,
            lowres_label=lowres_label,
            thickness_desc=thickness_desc,
        )
    return out


def stack_lowres(thick_roi: dict[str, tuple[Path, Path]]) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate L then R -> (thickness (2n,), roi (2n,) bool)."""
    thick_l = first_array(thick_roi["L"][0]).astype(np.float32)
    thick_r = first_array(thick_roi["R"][0]).astype(np.float32)
    roi_l = first_array(thick_roi["L"][1]) > 0.5
    roi_r = first_array(thick_roi["R"][1]) > 0.5
    return (
        np.concatenate([thick_l, thick_r]),
        np.concatenate([roi_l, roi_r]),
    )


def session_dirs(subject_session_dir: Path) -> tuple[Path, Path]:
    """(anat_dir, xfm_dir) for a .../sub-XXX/ses-YYY directory."""
    return subject_session_dir / "anat", subject_session_dir / "xfm"


def validate_sessions(
    sessions: list[Path], thickness_desc: str = "thickness"
) -> tuple[list[Path], list[tuple[str, str]]]:
    """Filter out sessions missing a required file or with an unreadable one
    (a truncated GIFTI = a partial download). Returns (good, [(id, reason)])."""
    thick_pat = (
        "*_hemi-{h}_thickness.shape.gii"
        if thickness_desc in ("", "thickness")
        else f"*_hemi-{{h}}_{thickness_desc}.shape.gii"
    )
    required = [
        ("anat", thick_pat),
        ("anat", "*_hemi-{h}_desc-medialwall_mask.shape.gii"),
        ("anat", "*_hemi-{h}_midthickness.surf.gii"),
        ("xfm", "*_hemi-{h}_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii"),
    ]
    good: list[Path] = []
    bad: list[tuple[str, str]] = []
    for session_dir in sessions:
        anat_dir, xfm_dir = session_dirs(session_dir)
        subject = f"{session_dir.parent.name}_{session_dir.name}"
        problem: str | None = None
        for sub_dir_name, pattern in required:
            directory = anat_dir if sub_dir_name == "anat" else xfm_dir
            for hemisphere in ("left", "right"):
                try:
                    path = find_one(directory, pattern.format(h=hemisphere))
                    _ = nib.load(str(path)).darrays[0].data.shape
                except Exception as error:  # noqa: BLE001
                    problem = f"{pattern.format(h=hemisphere)} -- {type(error).__name__}: {str(error)[:90]}"
                    break
            if problem:
                break
        if problem:
            bad.append((subject, problem))
        else:
            good.append(session_dir)
    return good, bad
