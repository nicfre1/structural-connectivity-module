#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Vertex-wise cortical thickness for ONE subject, native surface -> the common
low-resolution grid (see build_vertex_grid_assets.py).

See others/resample_subject_thickness_vertex.md for details.

Input : a subject's session directory, containing anat/ and xfm/
        (standard structural-pipeline / BIDS-derivatives layout).
Output: an .npz with
          thickness : (2N,)  vertex thickness, hemi L then R, low-res grid
          roi       : (2N,)  bool, True = valid cortical vertex
          n_per_hemi: int
          subject   : str
          lowres    : str (e.g. "10k")
        + the intermediate L/R .shape.gii files, kept in --work_dir
          (build_cohort_thickness_vertex.py reuses them for smoothing on
          the group surface).

No smoothing here -- surface smoothing happens once, in
build_cohort_thickness_vertex.py, on the group-average midthickness.

Example:
resample_subject_thickness_vertex.py \\
    /path/to/derivatives/sub-XXX/ses-YYY \\
    --assets_dir ./assets --template_dir /path/to/dhcp-symmetric-template \\
    --wb_command /Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command \\
    --lowres 10k --output vertex_thickness/sub-XXX_ses-YYY_thickness_vertex.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from vertex_resample_utils import resample_subject, session_dirs, stack_lowres


def build_subject_vertex_npz(
    *,
    session_dir: Path,
    assets_dir: Path,
    template_dir: Path,
    wb_command: str,
    lowres: str,
    output_npz: Path,
    work_dir: Path | None = None,
    thickness_desc: str = "thickness",
) -> dict:
    anat_dir, xfm_dir = session_dirs(session_dir)
    subject = f"{session_dir.parent.name}_{session_dir.name}"
    work_dir = work_dir or (output_npz.parent / "_work" / subject)

    thick_roi = resample_subject(
        anat_dir=anat_dir,
        xfm_dir=xfm_dir,
        template_dir=template_dir,
        assets_dir=assets_dir,
        work_dir=work_dir,
        wb_command=wb_command,
        lowres_label=lowres,
        thickness_desc=thickness_desc,
    )
    thickness, roi = stack_lowres(thick_roi)
    n_per_hemi = thickness.size // 2

    output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_npz,
        thickness=thickness.astype(np.float32),
        roi=roi,
        n_per_hemi=np.int64(n_per_hemi),
        subject=subject,
        lowres=lowres,
        thick_shape_L=str(thick_roi["L"][0]),
        thick_shape_R=str(thick_roi["R"][0]),
        roi_shape_L=str(thick_roi["L"][1]),
        roi_shape_R=str(thick_roi["R"][1]),
    )
    valid = int(roi.sum())
    print(
        f"{subject}: {thickness.size} vertices ({n_per_hemi}/hemi), "
        f"{valid} valid cortical, "
        f"thickness {np.nanmin(thickness[roi]):.2f}-{np.nanmax(thickness[roi]):.2f} mm"
    )
    return {"subject": subject, "thickness": thickness, "roi": roi, "work_dir": work_dir}


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("session_dir", type=Path,
                    help=".../sub-XXX/ses-YYY (must contain anat/ and xfm/)")
    p.add_argument("--assets_dir", type=Path, default=Path("assets"))
    p.add_argument("--template_dir", required=True, type=Path)
    p.add_argument("--wb_command", required=True, type=str)
    p.add_argument("--lowres", default="10k", help="grid label (see build_vertex_grid_assets.py)")
    p.add_argument("--thickness_desc", default="thickness",
                    help="'thickness' (default) or e.g. 'desc-corr_thickness'")
    p.add_argument("--output", required=True, type=Path)
    return p


def main():
    args = _build_arg_parser().parse_args()

    if not args.session_dir.exists():
        raise FileNotFoundError(args.session_dir)

    build_subject_vertex_npz(
        session_dir=args.session_dir,
        assets_dir=args.assets_dir,
        template_dir=args.template_dir,
        wb_command=args.wb_command,
        lowres=args.lowres,
        output_npz=args.output,
        thickness_desc=args.thickness_desc,
    )


if __name__ == "__main__":
    main()
