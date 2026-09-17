#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Render the first four vertex-wise gradients on a group cortical surface.

Vertex-wise counterpart of plot_group_gradients_brain.py: instead of
enigmatoolbox's bundled fsa5 (parcel-based) surface, this paints values
directly on a group midthickness mesh supplied by the caller (built by the
user's own vertex pipeline -- there is no bundled generic vertex mesh here).

Needs a VTK stack (see others/plot_group_gradients_brain_vertex.md and
requirements-brain.txt); run with a VTK-capable interpreter.

Example:
plot_group_gradients_brain_vertex.py \\
    --gradients-npz group_gradients.npz \\
    --left-surface groupmid_hemi-L.surf.gii --right-surface groupmid_hemi-R.surf.gii \\
    --output group_gradients_brain.png
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
from brainspace.mesh.mesh_creation import build_polydata
from brainspace.plotting import plot_hemispheres


def load_surface(path):
    image = nib.load(str(path))
    return (
        np.asarray(image.darrays[0].data, dtype=float),
        np.asarray(image.darrays[1].data, dtype=np.int32),
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--gradients-npz", required=True, type=Path)
    parser.add_argument("--left-surface", required=True, type=Path)
    parser.add_argument("--right-surface", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--n-gradients", type=int, default=4)
    parser.add_argument(
        "--label-text", nargs="+", default=["Grad1", "Grad2", "Grad3", "Grad4"]
    )
    parser.add_argument("--zoom", type=float, default=1.55)
    parser.add_argument("--size", type=int, nargs=2, default=(1200, 400))
    parser.add_argument("--cmap", default="viridis_r")
    args = parser.parse_args()

    for path in (args.gradients_npz, args.left_surface, args.right_surface):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    data = np.load(args.gradients_npz, allow_pickle=True)
    values = data["gradients_rescaled"] if "gradients_rescaled" in data else data["gradients"]
    n_per_hemi = int(data["n_per_hemi"])
    n_grad = min(args.n_gradients, values.shape[1])

    left_points, left_faces = load_surface(args.left_surface)
    right_points, right_faces = load_surface(args.right_surface)
    if left_points.shape[0] != n_per_hemi or right_points.shape[0] != n_per_hemi:
        raise RuntimeError(
            f"surface has {left_points.shape[0]}/{right_points.shape[0]} vertices, "
            f"expected {n_per_hemi} per hemisphere -- wrong grid?"
        )
    surf_lh = build_polydata(left_points, left_faces)
    surf_rh = build_polydata(right_points, right_faces)

    array_list = [
        np.concatenate([values[:n_per_hemi, k], values[n_per_hemi:, k]])
        for k in range(n_grad)
    ]

    plot_hemispheres(
        surf_lh, surf_rh,
        array_name=array_list,
        cmap=args.cmap,
        color_bar=True,
        label_text=args.label_text[:n_grad],
        zoom=args.zoom,
        size=tuple(args.size),
        interactive=False,
        offscreen=True,
        screenshot=True,
        filename=str(args.output),
        transparent_bg=False,
        background=(1, 1, 1),
    )
    print("Saved brain plot to:", args.output)


if __name__ == "__main__":
    main()
