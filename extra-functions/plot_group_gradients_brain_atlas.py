#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Render region-wise gradients on a cortical surface, for ANY atlas: the caller
supplies the surface and the atlas annotation for that surface.

Atlas-agnostic counterpart of plot_group_gradients_brain.py (which only knows
Desikan-Killiany on fsaverage5): the values come from a gradient CSV (one row
per region, columns Gradient_1 ...), the geometry from --left-surface /
--right-surface, and the vertex -> region assignment from --left-annot /
--right-annot (FreeSurfer .annot, or a GIFTI label file). Region names in the
CSV are matched to the annotation's label names after normalisation (case,
hemisphere tokens such as lh/rh/L/R, "ctx" prefix); use --name-map for
anything the normalisation cannot reconcile.

Needs a VTK stack (see others/plot_group_gradients_brain_atlas.md and
requirements-brain.txt); run with a VTK-capable interpreter.

Example:
plot_group_gradients_brain_atlas.py --gradient-csv group_gradients_rescaled_minus1_1.csv \\
    --left-surface lh.surf.gii --right-surface rh.surf.gii \\
    --left-annot lh.myatlas.annot --right-annot rh.myatlas.annot \\
    --output group_gradients_brain.png
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

HEMI_TOKENS = {"lh": "lh", "l": "lh", "left": "lh", "rh": "rh", "r": "rh", "right": "rh"}
IGNORED_LABELS = {"unknown", "medialwall", "medial_wall", "corpuscallosum", "???", "background"}


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


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    p.add_argument("--gradient-csv", required=True, type=Path,
                   help="CSV, first column = region names, then Gradient_1 ...")
    p.add_argument("--left-surface", required=True, type=Path)
    p.add_argument("--right-surface", required=True, type=Path)
    p.add_argument("--left-annot", required=True, type=Path)
    p.add_argument("--right-annot", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--gradient-columns", nargs="+",
                   default=["Gradient_1", "Gradient_2", "Gradient_3", "Gradient_4"])
    p.add_argument("--label-text", nargs="+", default=None,
                   help="Row labels (default: Grad1 Grad2 ...).")
    p.add_argument("--name-map", type=Path, default=None,
                   help="JSON {csv region name: annotation label name} for "
                        "names the automatic matching cannot reconcile.")
    p.add_argument("--cmap", default="viridis_r")
    p.add_argument("--nan-color", type=float, nargs=3, default=(0.85, 0.85, 0.85),
                   metavar=("R", "G", "B"),
                   help="RGB (0-1) for vertices with no value (unlabelled, "
                        "medial wall, unmatched regions). Default light grey.")
    p.add_argument("--color-range", type=float, nargs=2, default=None,
                   metavar=("MIN", "MAX"),
                   help="Fixed colour scale shared by all rows (default: each "
                        "row autoscaled to its data).")
    p.add_argument("--zoom", type=float, default=1.55)
    p.add_argument("--size", type=int, nargs=2, default=(1200, 400))
    return p


def main():
    args = _build_arg_parser().parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # VTK-dependent imports stay local so the matching logic above can be
    # imported (and tested) without a VTK stack.
    from brainspace.mesh.mesh_creation import build_polydata
    from brainspace.plotting import plot_hemispheres

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

    plot_hemispheres(
        surf_lh, surf_rh,
        array_name=array_list,
        cmap=args.cmap,
        nan_color=(*args.nan_color, 1.0),
        color_range=tuple(args.color_range) if args.color_range else None,
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
    )
    print("Saved brain plot to:", args.output)


if __name__ == "__main__":
    main()
