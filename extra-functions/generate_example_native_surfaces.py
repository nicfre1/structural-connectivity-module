#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a fully synthetic "template" + cohort "derivatives" tree, in the
exact file layout build_vertex_grid_assets.py / resample_subject_thickness_
vertex.py / build_cohort_thickness_vertex.py expect, so the whole native
surface -> low-resolution grid resampling chain can be smoke-tested with
Connectome Workbench and no real anatomy at all.

The surfaces are plain icospheres (wb_command -surface-create-sphere) and
the metrics are random values with the right shape -- geometrically valid
inputs for wb_command, but not anatomically meaningful. This only proves the
resampling *plumbing* is correct; run the same three scripts against a real
dHCP symmetric template and real subject derivatives to get a real result.

Needs Connectome Workbench (`wb_command`).

See others/generate_example_native_surfaces.md for details.

Example:
generate_example_native_surfaces.py --out_dir example_native --n_subjects 3 \\
    --wb_command /Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command
"""

import argparse
from pathlib import Path

import numpy as np

from vertex_resample_utils import HEMISPHERES, run, save_shape

TEMPLATE_32K = 162  # stand-in for the real dHCP template's 32k mesh


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--out_dir", default="example_native",
                    help="Directory where template/ and derivatives/ are written.")
    p.add_argument("--wb_command", required=True, type=str)
    p.add_argument("--n_subjects", type=int, default=3)
    p.add_argument("--n_native", type=int, default=42,
                    help="Vertices on the synthetic native mesh per hemisphere.")
    p.add_argument("--seed", type=int, default=0)
    return p


def _make_hemisphere_spheres(wb_command, n_vertices, right_path, left_path):
    run(wb_command, "-surface-create-sphere", n_vertices, right_path)
    run(wb_command, "-surface-flip-lr", right_path, left_path)
    run(wb_command, "-set-structure", right_path, "CORTEX_RIGHT")
    run(wb_command, "-set-structure", left_path, "CORTEX_LEFT")


def main():
    args = _build_arg_parser().parse_args()
    rng = np.random.default_rng(args.seed)

    out_dir = Path(args.out_dir)
    template_dir = out_dir / "template"
    derivatives_root = out_dir / "derivatives"
    template_dir.mkdir(parents=True, exist_ok=True)
    derivatives_root.mkdir(parents=True, exist_ok=True)

    # --- synthetic "week-40 dHCPsym 32k" template ---------------------------
    print("=== synthetic template ===")
    sphere_r32k = template_dir / "week-40_hemi-right_space-dhcpSym_dens-32k_sphere.surf.gii"
    sphere_l32k = template_dir / "week-40_hemi-left_space-dhcpSym_dens-32k_sphere.surf.gii"
    _make_hemisphere_spheres(args.wb_command, TEMPLATE_32K, sphere_r32k, sphere_l32k)

    area_lr = rng.uniform(0.5, 1.5, size=TEMPLATE_32K)
    save_shape(area_lr, template_dir / "week-40_hemi-LR_space-dhcpSym_dens-32k_surface_area.shape.gii", "area")
    sulc_lr = rng.normal(0.0, 1.0, size=TEMPLATE_32K)
    save_shape(sulc_lr, template_dir / "week-40_hemi-LR_space-dhcpSym_dens-32k_sulc.shape.gii", "sulc")

    # --- synthetic subject derivatives --------------------------------------
    print("\n=== synthetic subjects ===")
    for i in range(args.n_subjects):
        subject = f"sub-EX{i + 1:03d}"
        session = "ses-1"
        anat_dir = derivatives_root / subject / session / "anat"
        xfm_dir = derivatives_root / subject / session / "xfm"
        anat_dir.mkdir(parents=True, exist_ok=True)
        xfm_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"{subject}_{session}"

        for hemisphere, short in HEMISPHERES:
            native_sphere = anat_dir / f"{prefix}_hemi-{hemisphere}_midthickness.surf.gii"
            if hemisphere == "right":
                run(args.wb_command, "-surface-create-sphere", args.n_native, native_sphere)
                run(args.wb_command, "-set-structure", native_sphere, "CORTEX_RIGHT")
            else:
                tmp_right = anat_dir / f"_tmp_native_R_{prefix}.surf.gii"
                run(args.wb_command, "-surface-create-sphere", args.n_native, tmp_right)
                run(args.wb_command, "-surface-flip-lr", tmp_right, native_sphere)
                run(args.wb_command, "-set-structure", native_sphere, "CORTEX_LEFT")
                tmp_right.unlink(missing_ok=True)

            # a synthetic "already on the target sphere" registration: the
            # native mesh itself is a valid sphere, so reusing it satisfies
            # the vertex-count / sphere-validity requirements of the
            # resampling step without claiming any real registration.
            registered_sphere = xfm_dir / (
                f"{prefix}_hemi-{hemisphere}_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii"
            )
            registered_sphere.write_bytes(native_sphere.read_bytes())

            thickness = rng.uniform(1.5, 4.0, size=args.n_native)
            save_shape(thickness, anat_dir / f"{prefix}_hemi-{hemisphere}_thickness.shape.gii", "thickness")
            medialwall = np.ones(args.n_native)
            save_shape(medialwall, anat_dir / f"{prefix}_hemi-{hemisphere}_desc-medialwall_mask.shape.gii", "roi")

        print(f"  {subject}/{session}")

    print(f"\nWrote synthetic template -> {template_dir}")
    print(f"Wrote {args.n_subjects} synthetic subject(s) -> {derivatives_root}")


if __name__ == "__main__":
    main()
