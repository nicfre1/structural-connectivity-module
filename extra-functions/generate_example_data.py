#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate synthetic example data to try the extra-functions pipelines without
any real data. Three kinds:

  drawem32  synthetic DrawEM32 thickness TSVs (DrawEM32 -> MCRIBS -> gradients)
  vertex    synthetic vertex-wise cohort .npz (compute_covariance_gradients_vertex.py)
  native    synthetic template + subject surfaces for the native-surface
            resampling chain (needs Connectome Workbench)

See others/generate_example_data.md for details.

Examples:
generate_example_data.py drawem32 --out_dir example_data --n_subjects 25
generate_example_data.py vertex --out_dir example_data_vertex --n_subjects 25
generate_example_data.py native --out_dir example_native --wb_command /path/to/wb_command
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_PARCELS = 40  # drawem32: parcels per hemisphere
N_PER_HEMI = 162  # vertex: vertices per hemisphere
TEMPLATE_32K = 162  # native: stand-in for the real dHCP template's 32k mesh


def _latent_axes(n):
    position = np.linspace(0.0, 1.0, n)
    return np.cos(np.pi * position), np.sin(2.0 * np.pi * position)


# ---------------------------------------------------------------------------
# drawem32
# ---------------------------------------------------------------------------
def _hemisphere_frame(rng, subjects, id_column, hemi, suffix,
                      baseline, scale_1, scale_2, noise_sd,
                      score_1, score_2, subject_offset):
    axis_1, axis_2 = _latent_axes(N_PARCELS)
    n_subjects = len(subjects)

    data = {id_column: subjects}
    for idx in range(N_PARCELS):
        label = f"dwm_{idx + 1:02d}"
        values = (
            baseline
            + subject_offset
            + scale_1 * axis_1[idx] * score_1
            + scale_2 * axis_2[idx] * score_2
            + rng.normal(0.0, noise_sd, size=n_subjects)
        )
        data[f"{hemi}_{label}{suffix}"] = np.round(values, 4)
    return pd.DataFrame(data)


def generate_drawem32(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    subjects = [f"sub-{i + 1:03d}" for i in range(args.n_subjects)]
    n = args.n_subjects

    score_1 = rng.normal(1.0, 0.25, size=n)
    score_2 = rng.normal(1.0, 0.35, size=n)
    subject_offset = rng.normal(0.0, 0.12, size=n)

    for hemi in ("lh", "rh"):
        thickness = _hemisphere_frame(
            rng, subjects, args.id_column, hemi, "_thickness",
            baseline=2.5, scale_1=0.55, scale_2=0.22, noise_sd=0.08,
            score_1=score_1, score_2=score_2, subject_offset=subject_offset,
        )
        thickness_path = out_dir / f"drawem32_thickness_{hemi}_stats.tsv"
        thickness.to_csv(thickness_path, sep="\t", index=False)
        print(f"Wrote {thickness_path}  ({thickness.shape[0]} subjects, "
              f"{thickness.shape[1] - 1} parcels)")

        if args.with_area:
            area = _hemisphere_frame(
                rng, subjects, args.id_column, hemi, "_area",
                baseline=900.0, scale_1=260.0, scale_2=110.0, noise_sd=40.0,
                score_1=score_1, score_2=score_2,
                subject_offset=subject_offset * 200.0,
            )
            area_path = out_dir / f"drawem32_area_{hemi}_stats.tsv"
            area.to_csv(area_path, sep="\t", index=False)
            print(f"Wrote {area_path}")


# ---------------------------------------------------------------------------
# vertex
# ---------------------------------------------------------------------------
def generate_vertex(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    n = args.n_subjects
    n_per_hemi = args.n_per_hemi
    subject_ids = np.array([f"sub-{i + 1:03d}" for i in range(n)])

    axis_1, axis_2 = _latent_axes(n_per_hemi)
    score_1 = rng.normal(1.0, 0.25, size=n)
    score_2 = rng.normal(1.0, 0.35, size=n)
    subject_offset = rng.normal(0.0, 0.12, size=n)

    hemi = np.array(["L"] * n_per_hemi + ["R"] * n_per_hemi)
    thickness_full = np.empty((n, 2 * n_per_hemi), dtype=np.float32)
    for side_index in (0, 1):
        sl = slice(side_index * n_per_hemi, (side_index + 1) * n_per_hemi)
        thickness_full[:, sl] = (
            2.5
            + subject_offset[:, None]
            + 0.55 * axis_1[None, :] * score_1[:, None]
            + 0.22 * axis_2[None, :] * score_2[:, None]
            + rng.normal(0.0, 0.08, size=(n, n_per_hemi))
        )

    # medial-wall-like mask: a contiguous low-index block per hemisphere is
    # excluded, same on every subject (as the real group cortical mask is).
    n_excluded = int(round(args.medial_wall_fraction * n_per_hemi))
    mask = np.ones(2 * n_per_hemi, dtype=bool)
    if n_excluded > 0:
        mask[:n_excluded] = False
        mask[n_per_hemi:n_per_hemi + n_excluded] = False

    out_npz = out_dir / "cohort_thickness_vertex_dens-demo.npz"
    np.savez_compressed(
        out_npz,
        thickness=thickness_full[:, mask],
        mask=mask,
        hemi=hemi,
        subject_ids=subject_ids,
        fwhm=np.float64(0.0),
        lowres="demo",
        min_coverage=np.float64(1.0),
        n_per_hemi=np.int64(n_per_hemi),
    )
    print(f"Wrote {out_npz}  ({n} subjects, {int(mask.sum())} valid vertices "
          f"of {mask.size}, {n_per_hemi} per hemisphere)")


# ---------------------------------------------------------------------------
# native
# ---------------------------------------------------------------------------
def generate_native(args):
    # nibabel/Workbench helpers are only needed for this kind.
    from vertex_resample_utils import HEMISPHERES, run, save_shape

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
    run(args.wb_command, "-surface-create-sphere", TEMPLATE_32K, sphere_r32k)
    run(args.wb_command, "-surface-flip-lr", sphere_r32k, sphere_l32k)
    run(args.wb_command, "-set-structure", sphere_r32k, "CORTEX_RIGHT")
    run(args.wb_command, "-set-structure", sphere_l32k, "CORTEX_LEFT")

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

        for hemisphere, _short in HEMISPHERES:
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


# ---------------------------------------------------------------------------
def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    kinds = p.add_subparsers(dest="kind", required=True, metavar="{drawem32,vertex,native}")

    d = kinds.add_parser("drawem32", help="synthetic DrawEM32 thickness TSVs")
    d.add_argument("--out_dir", default="example_data",
                   help="Directory where the synthetic TSV files are written.")
    d.add_argument("--n_subjects", type=int, default=25,
                   help="Number of synthetic subjects.")
    d.add_argument("--id_column", default="sample",
                   help="Subject identifier column.")
    d.add_argument("--seed", type=int, default=0,
                   help="Random seed for reproducibility.")
    d.add_argument("--with_area", action="store_true",
                   help="Also write matching surface-area files.")
    d.set_defaults(func=generate_drawem32)

    v = kinds.add_parser("vertex", help="synthetic vertex-wise cohort .npz")
    v.add_argument("--out_dir", default="example_data_vertex",
                   help="Directory where the synthetic cohort .npz is written.")
    v.add_argument("--n_subjects", type=int, default=25,
                   help="Number of synthetic subjects.")
    v.add_argument("--n_per_hemi", type=int, default=N_PER_HEMI,
                   help="Vertices per hemisphere on the synthetic grid.")
    v.add_argument("--medial_wall_fraction", type=float, default=0.1,
                   help="Fraction of vertices per hemisphere excluded as "
                        "non-cortical (mask=False), same role as the medial "
                        "wall in the real pipeline.")
    v.add_argument("--seed", type=int, default=0,
                   help="Random seed for reproducibility.")
    v.set_defaults(func=generate_vertex)

    n = kinds.add_parser("native", help="synthetic template + subject surfaces (needs Workbench)")
    n.add_argument("--out_dir", default="example_native",
                   help="Directory where template/ and derivatives/ are written.")
    n.add_argument("--wb_command", required=True, type=str)
    n.add_argument("--n_subjects", type=int, default=3)
    n.add_argument("--n_native", type=int, default=42,
                   help="Vertices on the synthetic native mesh per hemisphere.")
    n.add_argument("--seed", type=int, default=0)
    n.set_defaults(func=generate_native)
    return p


def main():
    args = _build_arg_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
