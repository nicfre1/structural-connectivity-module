#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Render the first gradients on a cortical surface. Three modes:

  dk      Desikan-Killiany regions on fsaverage5 (ENIGMA's aparc_fsa5);
          nothing to provide but the gradient CSV. This is the default mode.
  atlas   region-wise gradients for ANY atlas: you provide the surface and the
          atlas annotation (FreeSurfer .annot or GIFTI labels) for it.
  vertex  vertex-wise gradients (.npz from compute_covariance_gradients_vertex.py)
          on a group midthickness surface you provide.

Needs a VTK stack (see others/plot_group_gradients_brain.md and
requirements-brain.txt); run with a VTK-capable interpreter.

Examples:
plot_group_gradients_brain.py --gradient-csv group_gradients_rescaled_minus1_1.csv --output brain.png
plot_group_gradients_brain.py atlas --gradient-csv g.csv --left-surface lh.surf.gii --right-surface rh.surf.gii \\
    --left-annot lh.myatlas.annot --right-annot rh.myatlas.annot --output brain.png
plot_group_gradients_brain.py vertex --gradients-npz group_gradients.npz \\
    --left-surface groupmid_L.surf.gii --right-surface groupmid_R.surf.gii --output brain.png
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODES = ("dk", "atlas", "vertex")
DEFAULT_MODE = "dk"

HEMI_TOKENS = {"lh": "lh", "l": "lh", "left": "lh", "rh": "rh", "r": "rh", "right": "rh"}
IGNORED_LABELS = {"unknown", "medialwall", "medial_wall", "corpuscallosum", "???", "background"}


# ---------------------------------------------------------------------------
# atlas matching (numpy/pandas only: importable and testable without VTK)
# ---------------------------------------------------------------------------
def normalize_name(name):
    """-> (hemisphere or None, normalised name without hemisphere/ctx tokens)."""
    tokens = [t for t in re.split(r"[_\-.\s]+", str(name).strip().lower()) if t]
    hemi, kept = None, []
    for token in tokens:
        if token in HEMI_TOKENS and hemi is None:
            hemi = HEMI_TOKENS[token]
        elif token != "ctx":
            kept.append(token)
    return hemi, "_".join(kept)


def load_surface(path):
    import nibabel as nib

    path = Path(path)
    if path.suffix == ".gii":
        image = nib.load(str(path))
        points = faces = None
        for darray in image.darrays:
            if darray.intent == nib.nifti1.intent_codes["NIFTI_INTENT_POINTSET"]:
                points = np.asarray(darray.data, dtype=float)
            elif darray.intent == nib.nifti1.intent_codes["NIFTI_INTENT_TRIANGLE"]:
                faces = np.asarray(darray.data, dtype=np.int32)
        if points is None or faces is None:
            raise ValueError(f"{path}: no POINTSET/TRIANGLE arrays found")
        return points, faces
    points, faces = nib.freesurfer.read_geometry(str(path))
    return np.asarray(points, dtype=float), np.asarray(faces, dtype=np.int32)


def load_annotation(path):
    """-> (per-vertex label index, label names); index -1 = no label."""
    import nibabel as nib

    path = Path(path)
    if path.suffix == ".annot":
        labels, _ctab, names = nib.freesurfer.read_annot(str(path))
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in names]
        return np.asarray(labels, dtype=int), names
    image = nib.load(str(path))
    table = image.labeltable.get_labels_as_dict()
    codes = np.asarray(image.darrays[0].data).ravel().astype(int)
    keys = sorted(table)
    index_of = {key: i for i, key in enumerate(keys)}
    labels = np.array([index_of.get(c, -1) for c in codes], dtype=int)
    return labels, [table[k] for k in keys]


def paint_hemisphere(values, annot_labels, annot_names, hemi, name_map):
    """Map region values onto vertices; returns (vertex values, report dict)."""
    lookup, region_hemis = {}, {}
    for region, value in values.items():
        region_hemi, key = normalize_name(name_map.get(region, region))
        lookup[(region_hemi, key)] = (region, value)
        region_hemis[region] = region_hemi

    label_values = np.full(len(annot_names), np.nan)
    matched, unmatched_labels, used = [], [], set()
    for i, name in enumerate(annot_names):
        _h, key = normalize_name(name)
        if key in IGNORED_LABELS or not key:
            continue
        hit = lookup.get((hemi, key)) or lookup.get((None, key))
        if hit is None:
            unmatched_labels.append(name)
            continue
        region, value = hit
        label_values[i] = value
        used.add(region)
        matched.append(name)

    vertex_values = np.full(annot_labels.shape, np.nan)
    valid = annot_labels >= 0
    vertex_values[valid] = label_values[annot_labels[valid]]
    unused_regions = [
        r for r in values
        if r not in used and region_hemis[r] in (hemi, None)
    ]
    return vertex_values, {
        "matched": matched,
        "unmatched_labels": unmatched_labels,
        "unused_regions": unused_regions,
    }


# ---------------------------------------------------------------------------
# rendering (VTK)
# ---------------------------------------------------------------------------
def render(surf_lh, surf_rh, array_list, labels, args, **extra):
    from brainspace.plotting import plot_hemispheres

    plot_hemispheres(
        surf_lh, surf_rh,
        array_name=array_list,
        cmap=args.cmap,
        color_bar=True,
        label_text=labels,
        zoom=args.zoom,
        size=tuple(args.size),
        interactive=False,
        offscreen=True,
        screenshot=True,
        filename=str(args.output),
        transparent_bg=False,
        background=(1, 1, 1),
        **extra,
    )
    print("Saved brain plot to:", args.output)


def run_dk(args):
    import enigmatoolbox
    import nibabel as nib
    from brainspace.mesh.mesh_creation import build_polydata
    from enigmatoolbox.utils.parcellation import parcel_to_surface

    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.gradient_csv, index_col=0)

    array_list = [
        parcel_to_surface(frame[column].to_numpy(dtype=float), "aparc_fsa5")
        for column in args.gradient_columns
    ]

    surf_dir = os.path.join(
        os.path.dirname(enigmatoolbox.__file__), "datasets", "surfaces"
    )
    lh = nib.load(os.path.join(surf_dir, "fsa5_lh.surf.gii"))
    rh = nib.load(os.path.join(surf_dir, "fsa5_rh.surf.gii"))
    surf_lh = build_polydata(
        np.asarray(lh.darrays[0].data, dtype=float),
        np.asarray(lh.darrays[1].data, dtype=np.int32),
    )
    surf_rh = build_polydata(
        np.asarray(rh.darrays[0].data, dtype=float),
        np.asarray(rh.darrays[1].data, dtype=np.int32),
    )
    labels = args.label_text or [f"Grad{i + 1}" for i in range(len(args.gradient_columns))]
    render(surf_lh, surf_rh, array_list, labels, args)


def run_vertex(args):
    from brainspace.mesh.mesh_creation import build_polydata

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
    labels = (args.label_text or [f"Grad{i + 1}" for i in range(n_grad)])[:n_grad]
    render(surf_lh, surf_rh, array_list, labels, args)


def run_atlas(args):
    from brainspace.mesh.mesh_creation import build_polydata

    args.output.parent.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.gradient_csv, index_col=0)
    name_map = json.loads(args.name_map.read_text()) if args.name_map else {}
    labels = args.label_text or [f"Grad{i + 1}" for i in range(len(args.gradient_columns))]

    hemispheres = {}
    for hemi, surf_path, annot_path in (
        ("lh", args.left_surface, args.left_annot),
        ("rh", args.right_surface, args.right_annot),
    ):
        points, faces = load_surface(surf_path)
        annot_labels, annot_names = load_annotation(annot_path)
        if annot_labels.size != points.shape[0]:
            raise SystemExit(
                f"{hemi}: annotation has {annot_labels.size} vertices but the "
                f"surface has {points.shape[0]} -- they must be on the same mesh."
            )
        hemispheres[hemi] = (points, faces, annot_labels, annot_names)

    arrays = {column: [] for column in args.gradient_columns}
    for hemi, (points, faces, annot_labels, annot_names) in hemispheres.items():
        report = None
        for column in args.gradient_columns:
            values = frame[column].to_dict()
            vertex_values, report = paint_hemisphere(
                values, annot_labels, annot_names, hemi, name_map
            )
            arrays[column].append(vertex_values)
        n_labels = len(report["matched"]) + len(report["unmatched_labels"])
        print(f"{hemi}: matched {len(report['matched'])}/{n_labels} annotation labels "
              f"to CSV regions")
        if not report["matched"]:
            raise SystemExit(
                f"{hemi}: no CSV region matched any annotation label. Example CSV "
                f"names: {list(frame.index[:3])}; example annotation names: "
                f"{annot_names[:3]}. Use --name-map."
            )
        if report["unmatched_labels"]:
            print(f"  annotation labels with no CSV region (left grey): "
                  f"{report['unmatched_labels'][:8]}"
                  f"{' ...' if len(report['unmatched_labels']) > 8 else ''}")
        if report["unused_regions"]:
            print(f"  CSV regions not found in the annotation (not drawn): "
                  f"{report['unused_regions'][:8]}"
                  f"{' ...' if len(report['unused_regions']) > 8 else ''}")

    surf_lh = build_polydata(hemispheres["lh"][0], hemispheres["lh"][1])
    surf_rh = build_polydata(hemispheres["rh"][0], hemispheres["rh"][1])
    array_list = [np.concatenate(arrays[column]) for column in args.gradient_columns]
    render(
        surf_lh, surf_rh, array_list, labels, args,
        nan_color=(*args.nan_color, 1.0),
        color_range=tuple(args.color_range) if args.color_range else None,
    )


# ---------------------------------------------------------------------------
def _add_common(p):
    p.add_argument("--output", required=True, type=Path, help="PNG to write.")
    p.add_argument("--label-text", nargs="+", default=None,
                   help="Row labels (default: Grad1 Grad2 ...).")
    p.add_argument("--cmap", default="viridis_r")
    p.add_argument("--zoom", type=float, default=1.55)
    p.add_argument("--size", type=int, nargs=2, default=(1200, 400))


def _add_gradient_columns(p):
    p.add_argument("--gradient-columns", nargs="+",
                   default=["Gradient_1", "Gradient_2", "Gradient_3", "Gradient_4"])


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    modes = p.add_subparsers(
        dest="mode", required=True, metavar="{dk,atlas,vertex}",
        help="Default mode when none is given: dk.",
    )

    d = modes.add_parser("dk", help="Desikan-Killiany on fsaverage5")
    d.add_argument("--gradient-csv", required=True, type=Path,
                   help="CSV, first column = the 68 DK regions (lh then rh), then Gradient_1 ...")
    _add_gradient_columns(d)
    _add_common(d)
    d.set_defaults(func=run_dk)

    a = modes.add_parser("atlas", help="any atlas: your surface + annotation")
    a.add_argument("--gradient-csv", required=True, type=Path,
                   help="CSV, first column = region names, then Gradient_1 ...")
    a.add_argument("--left-surface", required=True, type=Path)
    a.add_argument("--right-surface", required=True, type=Path)
    a.add_argument("--left-annot", required=True, type=Path)
    a.add_argument("--right-annot", required=True, type=Path)
    _add_gradient_columns(a)
    a.add_argument("--name-map", type=Path, default=None,
                   help="JSON {csv region name: annotation label name} for "
                        "names the automatic matching cannot reconcile.")
    a.add_argument("--nan-color", type=float, nargs=3, default=(0.85, 0.85, 0.85),
                   metavar=("R", "G", "B"),
                   help="RGB (0-1) for vertices with no value (unlabelled, "
                        "medial wall, unmatched regions). Default light grey.")
    a.add_argument("--color-range", type=float, nargs=2, default=None,
                   metavar=("MIN", "MAX"),
                   help="Fixed colour scale shared by all rows (default: each "
                        "row autoscaled to its data).")
    _add_common(a)
    a.set_defaults(func=run_atlas)

    v = modes.add_parser("vertex", help="vertex-wise gradients on your group surface")
    v.add_argument("--gradients-npz", required=True, type=Path)
    v.add_argument("--left-surface", required=True, type=Path)
    v.add_argument("--right-surface", required=True, type=Path)
    v.add_argument("--n-gradients", type=int, default=4)
    _add_common(v)
    v.set_defaults(func=run_vertex)
    return p


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] not in MODES and argv[0] not in ("-h", "--help"):
        argv = [DEFAULT_MODE, *argv]
    args = _build_arg_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
