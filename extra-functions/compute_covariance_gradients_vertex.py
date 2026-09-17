#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute individual and group covariance matrices, group and individual
gradients, and lambdas from a cohort vertex-thickness file (vertex-wise,
instead of the 68 Desikan-Killiany regions of compute_covariance_gradients.py).

Same equations as compute_covariance_gradients.py, one row per cortical
vertex of a common surface grid instead of one row per region. See
others/compute_covariance_gradients_vertex.md for details, and
generate_example_thickness_vertex.py to build a synthetic cohort .npz.

The cohort .npz must contain: thickness (n_subjects, n_valid_vertices),
mask (2*n_per_hemi,) bool, hemi (2*n_per_hemi,) 'L'/'R', subject_ids,
n_per_hemi. This is the schema written by the real dHCP vertex pipeline's
cohort-assembly step (native surface -> template grid, resampled with
Connectome Workbench) -- that resampling itself is out of scope here, since
it needs real per-subject surface data and is not data-free.

Example:
compute_covariance_gradients_vertex.py cohort_thickness_vertex_dens-10k.npz
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from brainspace.gradient import GradientMaps

from compute_covariance_gradients import (
    individualized_covariance_from_z,
    rescale_to_minus1_1,
    save_scree_figure,
)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        "cohort_npz",
        type=Path,
        help="Cohort vertex-thickness .npz (see generate_example_thickness_vertex.py).",
    )
    p.add_argument(
        "--out_dir",
        default="results_vertex",
        help="Directory where CSV/NPZ files and figures are written.",
    )
    p.add_argument(
        "--n_components",
        type=int,
        default=10,
        help="Number of gradient components.",
    )
    p.add_argument(
        "--sparsity",
        type=float,
        default=0.9,
        help="Sparsity passed to GradientMaps (BrainSpace default, same as "
             "the real vertex pipeline).",
    )
    p.add_argument(
        "--save_matrices",
        action="store_true",
        help="Also write each subject's V x V individualized covariance "
             "matrix as .npy float32 (large at real vertex counts).",
    )
    p.add_argument(
        "--figure_max_vertices",
        type=int,
        default=2000,
        help="Regular subsampling of the individualized-covariance heatmaps.",
    )
    p.add_argument(
        "--left_surface",
        type=Path,
        default=None,
        help="Group midthickness .surf.gii (left hemisphere), for the "
             "cortical-surface render. Omit to skip the render.",
    )
    p.add_argument(
        "--right_surface",
        type=Path,
        default=None,
        help="Group midthickness .surf.gii (right hemisphere).",
    )
    p.add_argument(
        "--no_brain_plot",
        action="store_true",
        help="Skip the cortical-surface render even if surfaces are given.",
    )
    p.add_argument(
        "-v",
        default="WARNING",
        const="INFO",
        nargs="?",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        dest="verbose",
        help="Verbosity level. Default: WARNING, or INFO when -v is given.",
    )
    return p


def load_cohort(cohort_npz):
    data = np.load(cohort_npz, allow_pickle=True)
    X = data["thickness"].astype(np.float64)
    mask = data["mask"].astype(bool)
    n_per_hemi = int(data["n_per_hemi"])
    subject_ids = list(map(str, data["subject_ids"]))
    return X, mask, n_per_hemi, subject_ids


def scatter_to_full(values, mask):
    """(n_valid, k) restricted to mask -> (mask.size, k) with NaN elsewhere."""
    full = np.full((mask.size, values.shape[1]), np.nan, dtype=np.float32)
    full[mask] = values
    return full


def save_individualized_covariance_heatmap(matrix, subject_id, n_subjects, n_vertices, step, path):
    small = matrix[::step, ::step]
    figure, axis = plt.subplots(figsize=(9, 9))
    image = axis.imshow(
        small, cmap="viridis", vmin=0, vmax=1,
        interpolation="nearest", aspect="equal",
    )
    figure.colorbar(image, ax=axis, label="Individualized covariance")
    axis.set_title(
        f"{subject_id}\nV={n_vertices} (shown 1/{step}), n_cohort={n_subjects}"
    )
    axis.set_xticks([])
    axis.set_yticks([])
    figure.tight_layout()
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def render_group_gradients_brain_vertex(gradients_npz, left_surface, right_surface, out_png):
    script = Path(__file__).with_name("plot_group_gradients_brain_vertex.py")
    try:
        completed = subprocess.run(
            [sys.executable, str(script),
             "--gradients-npz", str(gradients_npz),
             "--left-surface", str(left_surface),
             "--right-surface", str(right_surface),
             "--output", str(out_png)],
            capture_output=True, text=True, timeout=600,
        )
        if completed.returncode == 0 and out_png.exists():
            logging.warning("Saved cortical-surface render to %s", out_png)
            return
        logging.warning(
            "Cortical-surface render failed (exit %s). Last output: %s",
            completed.returncode,
            (completed.stderr or completed.stdout or "").strip().splitlines()[-1:]
            or "<none>",
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        logging.warning("Cortical-surface render could not run (%s).", exc)


def main():
    args = _build_arg_parser().parse_args()
    logging.getLogger().setLevel(logging.getLevelName(args.verbose))

    X, mask, n_per_hemi, subjects = load_cohort(args.cohort_npz)
    n_subjects, n_vertices = X.shape
    if n_subjects < 2:
        raise SystemExit(
            f"Cohort z-scoring needs at least 2 subjects (got {n_subjects})."
        )
    logging.warning(
        "Cohort: %d subjects, %d valid vertices (%d per hemisphere on the full grid)",
        n_subjects, n_vertices, n_per_hemi,
    )

    out_dir = Path(args.out_dir)
    indiv_cov_png_dir = out_dir / "individualized_covariance" / "png"
    group_dir = out_dir / "group_covariance"
    indiv_grad_dir = out_dir / "individual_gradients"
    for directory in (indiv_cov_png_dir, group_dir, indiv_grad_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if args.save_matrices:
        indiv_cov_npy_dir = out_dir / "individualized_covariance" / "npy"
        indiv_cov_npy_dir.mkdir(parents=True, exist_ok=True)

    X_z = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    logging.info(
        "z-score check: vertex means ~0 (%s), vertex stds ~1 (%s)",
        np.allclose(X_z.mean(axis=0), 0, atol=1e-10),
        np.allclose(X_z.std(axis=0), 1, atol=1e-10),
    )

    step = max(1, n_vertices // args.figure_max_vertices)
    individual_covariances = []
    for index, subject_id in enumerate(subjects):
        matrix = individualized_covariance_from_z(X_z[index])
        individual_covariances.append(matrix)
        if args.save_matrices:
            np.save(
                indiv_cov_npy_dir / f"{subject_id}_individualized_covariance.npy",
                matrix.astype(np.float32),
            )
        save_individualized_covariance_heatmap(
            matrix, subject_id, n_subjects, n_vertices, step,
            indiv_cov_png_dir / f"{subject_id}_individualized_covariance.png",
        )
    logging.warning(
        "Wrote %d individualized covariance heatmaps%s",
        n_subjects, " (+ .npy matrices)" if args.save_matrices else "",
    )

    covariance_group = np.corrcoef(X_z, rowvar=False)

    gm_group = GradientMaps(
        n_components=args.n_components, random_state=0, approach="dm", kernel=None,
    )
    gm_group.fit(covariance_group, sparsity=args.sparsity)
    group_grad = gm_group.gradients_

    full_group = scatter_to_full(group_grad, mask)
    group_gradients_npz = group_dir / "group_gradients.npz"
    np.savez_compressed(
        group_gradients_npz,
        gradients=full_group, gradients_rescaled=np.column_stack(
            [rescale_to_minus1_1(full_group[:, k]) for k in range(min(4, full_group.shape[1]))]
        ),
        lambdas=gm_group.lambdas_, mask=mask, n_per_hemi=np.int64(n_per_hemi),
    )
    pd.DataFrame(
        {f"Gradient_{i + 1}": group_grad[:, i] for i in range(group_grad.shape[1])}
    ).to_csv(group_dir / "group_gradients_valid_vertices.csv", index_label="valid_vertex")
    pd.DataFrame(
        gm_group.lambdas_,
        index=[f"Gradient_{i + 1}" for i in range(len(gm_group.lambdas_))],
        columns=["lambda"],
    ).to_csv(group_dir / "group_lambdas.csv")
    save_scree_figure(gm_group.lambdas_, group_dir / "group_lambdas.png")
    logging.warning("Wrote group gradients (npz + csv) and lambdas")

    failed = []
    for index, subject_id in enumerate(subjects):
        gm_indiv = GradientMaps(
            n_components=args.n_components, random_state=0,
            approach="dm", kernel=None, alignment="procrustes",
        )
        try:
            gm_indiv.fit(individual_covariances[index], reference=group_grad, sparsity=args.sparsity)
        except Exception as error:  # noqa: BLE001
            failed.append(subject_id)
            logging.warning("Skipping %s: %s", subject_id, error)
            continue
        aligned = np.asarray(gm_indiv.aligned_)
        if aligned.ndim == 3:
            aligned = aligned[0]

        full = scatter_to_full(aligned, mask)
        rescaled = np.column_stack(
            [rescale_to_minus1_1(full[:, k]) for k in range(min(4, full.shape[1]))]
        )
        np.savez_compressed(
            indiv_grad_dir / f"{subject_id}_gradients_aligned.npz",
            gradients=full, gradients_rescaled=rescaled, mask=mask,
            n_per_hemi=np.int64(n_per_hemi),
        )
        pd.DataFrame(
            {f"Gradient_{j + 1}": aligned[:, j] for j in range(aligned.shape[1])}
        ).to_csv(
            indiv_grad_dir / f"{subject_id}_gradients_aligned_valid_vertices.csv",
            index_label="valid_vertex",
        )
    if failed:
        logging.warning("%d subject(s) failed Procrustes alignment: %s", len(failed), failed)
    logging.warning("Wrote aligned individual gradients (npz + csv)")

    if not args.no_brain_plot and args.left_surface and args.right_surface:
        render_group_gradients_brain_vertex(
            group_gradients_npz, args.left_surface, args.right_surface,
            group_dir / "group_gradients_brain.png",
        )
    elif not args.no_brain_plot:
        logging.warning(
            "No --left_surface/--right_surface given: skipping the cortical-surface "
            "render (unlike the 68-region pipeline, there is no bundled generic "
            "vertex mesh to fall back on)."
        )

    logging.warning("Done. Outputs written to %s", out_dir)


if __name__ == "__main__":
    main()
