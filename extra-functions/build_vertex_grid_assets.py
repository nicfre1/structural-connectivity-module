#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build the low-resolution icosphere grid assets used by the vertex-wise
pipeline, from a dHCP symmetric surface template.

See others/build_vertex_grid_assets.md for details, and
generate_example_native_surfaces.py to exercise this without a real
template (a tiny synthetic one, for a smoke test of the whole chain).

Needs Connectome Workbench (`wb_command`) and a local copy of the dHCP
symmetric surface template (week-40, publicly released by the developing
Human Connectome Project) -- neither is shipped here; --template-dir must
point at your own copy.

Writes, in --output_dir (default ./assets):

  sphere_hemi-{L,R}_dhcpSym_dens-{N}k.surf.gii
      icosphere at ~N thousand vertices (wb_command -surface-create-sphere
      then -surface-flip-lr, so L and R are in vertex-to-vertex
      correspondence). The resampling target of the whole chain.
  area_hemi-{L,R}_dhcpSym_dens-32k.shape.gii
      per-vertex surface area of the template's 32k mesh, one copy per
      hemisphere (fixed GIFTI structure) -- the <new-area> of
      resample_subject_thickness_vertex.py's ADAP_BARY_AREA correction.
  sulc_hemi-{L,R}_dhcpSym_dens-{N}k.shape.gii
      template sulcal depth resampled onto the ~N k grid, for shading a
      cortical-surface render.

The template's "dens-32k" mesh IS the fs_LR 32k mesh (same faces, spheres
within ~1 mm): a generic icosphere is therefore already in the right frame,
exactly like HCP's standard downsampled spheres.

Example:
build_vertex_grid_assets.py --template-dir /path/to/dhcp-symmetric-template \\
    --wb-command /Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command \\
    --output_dir assets --n_vertices 10242
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from vertex_resample_utils import HEMISPHERES, run

VALID_ICO = {642: "0.6k", 2562: "2.5k", 10242: "10k", 40962: "40k"}


def lowres_label(n_vertices: int) -> str:
    if n_vertices not in VALID_ICO:
        raise SystemExit(
            "--n_vertices must be a valid icosahedral tessellation "
            f"(4^k*10+2): {sorted(VALID_ICO)}. Got: {n_vertices}"
        )
    return VALID_ICO[n_vertices]


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--template_dir", required=True, type=Path,
                    help="Directory with the dHCP symmetric template files "
                         "(week-40_*_sphere/area/sulc...).")
    p.add_argument("--wb_command", required=True, type=str,
                    help="Path to (or name of) the wb_command executable.")
    p.add_argument("--output_dir", type=Path, default=Path("assets"))
    p.add_argument("--n_vertices", type=int, default=10242,
                    help="Vertices per hemisphere (642 / 2562 / 10242 / 40962).")
    return p


def main():
    args = _build_arg_parser().parse_args()

    label = lowres_label(args.n_vertices)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. ~N k icosphere, L and R vertex-matched --------------------------
    tmp_right = args.output_dir / f"_tmp_sphere_R_{label}.surf.gii"
    sphere_r = args.output_dir / f"sphere_hemi-R_dhcpSym_dens-{label}.surf.gii"
    sphere_l = args.output_dir / f"sphere_hemi-L_dhcpSym_dens-{label}.surf.gii"

    run(args.wb_command, "-surface-create-sphere", args.n_vertices, tmp_right)
    run(args.wb_command, "-surface-flip-lr", tmp_right, sphere_l)
    shutil.move(str(tmp_right), str(sphere_r))
    run(args.wb_command, "-set-structure", sphere_r, "CORTEX_RIGHT")
    run(args.wb_command, "-set-structure", sphere_l, "CORTEX_LEFT")

    # --- 2. week-40 32k area, one copy per hemisphere -----------------------
    # the symmetric template ships one shared "hemi-LR" file.
    area_lr = next(
        p for p in args.template_dir.glob(
            "week-40_hemi-LR_space-dhcpSym_dens-32k_surface_area.shape.gii"
        )
        if not p.name.startswith("._")
    )
    for hemisphere, short in HEMISPHERES:
        structure = "CORTEX_LEFT" if hemisphere == "left" else "CORTEX_RIGHT"
        area_hemi = args.output_dir / f"area_hemi-{short}_dhcpSym_dens-32k.shape.gii"
        shutil.copyfile(area_lr, area_hemi)
        run(args.wb_command, "-set-structure", area_hemi, structure)

    # --- 3. week-40 sulc -> ~N k grid (render shading) -----------------------
    sulc_lr = next(
        p for p in args.template_dir.glob(
            "week-40_hemi-LR_space-dhcpSym_dens-32k_sulc.shape.gii"
        )
        if not p.name.startswith("._")
    )
    for hemisphere, short in HEMISPHERES:
        structure = "CORTEX_LEFT" if hemisphere == "left" else "CORTEX_RIGHT"
        dhcp32k_sphere = (
            args.template_dir
            / f"week-40_hemi-{hemisphere}_space-dhcpSym_dens-32k_sphere.surf.gii"
        )
        sulc_32k = args.output_dir / f"_tmp_sulc_{short}_32k.shape.gii"
        shutil.copyfile(sulc_lr, sulc_32k)
        run(args.wb_command, "-set-structure", sulc_32k, structure)
        sulc_low = args.output_dir / f"sulc_hemi-{short}_dhcpSym_dens-{label}.shape.gii"
        target_sphere = args.output_dir / f"sphere_hemi-{short}_dhcpSym_dens-{label}.surf.gii"
        run(
            args.wb_command, "-metric-resample",
            sulc_32k, dhcp32k_sphere, target_sphere, "BARYCENTRIC", sulc_low,
        )
        sulc_32k.unlink(missing_ok=True)

    print(f"\ndHCPsym dens-{label} grid ready in: {args.output_dir}")
    for path in sorted(args.output_dir.glob(f"*dens-{label}*")):
        print(f"  {path.name}")
    for path in sorted(args.output_dir.glob("area_hemi-*_dens-32k*")):
        print(f"  {path.name}")


if __name__ == "__main__":
    main()
