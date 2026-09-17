#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a synthetic cohort vertex-thickness .npz to test the vertex-wise
gradient pipeline without real surface data.

Same schema as the real vertex pipeline's cohort file (see
compute_covariance_gradients_vertex.py): thickness (n_subjects, n_valid),
mask (2*n_per_hemi,) bool, hemi (2*n_per_hemi,) 'L'/'R', subject_ids,
n_per_hemi.

See others/generate_example_thickness_vertex.md for details.

Example:
generate_example_thickness_vertex.py --out_dir example_data_vertex --n_subjects 25
"""

import argparse
from pathlib import Path

import numpy as np

N_PER_HEMI = 162


def _latent_axes(n_per_hemi):
    position = np.linspace(0.0, 1.0, n_per_hemi)
    axis_1 = np.cos(np.pi * position)
    axis_2 = np.sin(2.0 * np.pi * position)
    return axis_1, axis_2


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--out_dir", default="example_data_vertex",
                   help="Directory where the synthetic cohort .npz is written.")
    p.add_argument("--n_subjects", type=int, default=25,
                   help="Number of synthetic subjects.")
    p.add_argument("--n_per_hemi", type=int, default=N_PER_HEMI,
                   help="Vertices per hemisphere on the synthetic grid.")
    p.add_argument("--medial_wall_fraction", type=float, default=0.1,
                   help="Fraction of vertices per hemisphere excluded as "
                        "non-cortical (mask=False), same role as the medial "
                        "wall in the real pipeline.")
    p.add_argument("--seed", type=int, default=0,
                   help="Random seed for reproducibility.")
    return p


def main():
    args = _build_arg_parser().parse_args()

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


if __name__ == "__main__":
    main()
