# plot_group_gradients_brain_atlas.py

Render region-wise gradients on a cortical surface, for **any atlas**: the
caller supplies the surface and the atlas annotation for that surface.

Atlas-agnostic counterpart of `plot_group_gradients_brain.py`, which only
knows Desikan-Killiany on fsaverage5 (ENIGMA's `aparc_fsa5`). Here the
vertex -> region assignment comes from an annotation file you provide, so
Brainnetome, Schaefer, Glasser, a custom atlas, ... all work the same way, on
whichever surface the annotation was defined on.

## What you provide

| Input | Option | Notes |
| --- | --- | --- |
| Gradient values | `--gradient-csv` | First column = region names, then `Gradient_1` ... (e.g. the CSVs written by `compute_covariance_gradients.py`). |
| Surfaces | `--left-surface`, `--right-surface` | GIFTI `.surf.gii` / `.gii`, or a FreeSurfer surface file. |
| Atlas on the surface | `--left-annot`, `--right-annot` | FreeSurfer `.annot`, or a GIFTI label file. Must have exactly one entry per vertex of the matching surface. |

## How region names are matched

CSV region names and annotation label names are compared after
normalisation: lower-cased, split on `_ - .` and spaces, the `ctx` token
dropped, and a hemisphere token (`lh`, `rh`, `l`, `r`, `left`, `right`)
removed and remembered. So `ctx-lh-bankssts`, `lh_bankssts` and
`bankssts` all match the label `bankssts`; `Sch200_LH_001` matches
`Sch200_LH_001`; `SFG_L_6_1` matches `SFG_6_1`.

- A CSV row with a hemisphere token is drawn only on that hemisphere.
- A CSV row without one (for example region values already averaged over
  left and right) is drawn on both hemispheres.
- Annotation labels named `unknown`, `medialwall`, `corpuscallosum`, ... , or
  with no match, are left uncoloured (`--nan-color`, light grey by default).

For each hemisphere the script prints how many annotation labels were
matched, and lists unmatched labels and CSV regions that were not drawn. If
nothing matches at all it stops with an explicit message. For names the
normalisation cannot reconcile, pass `--name-map` with a JSON
`{"csv region name": "annotation label name"}`.

## Requirements

A VTK stack (`brainspace[plotting]`, `nibabel`, `vtk`; see
`../requirements-brain.txt`); run with a VTK-capable interpreter. The
matching logic itself needs only numpy/pandas and is importable without VTK
(that is what the unit tests exercise); nibabel and VTK are imported only
when reading surfaces/annotations and rendering.

## Options

| Option | Meaning |
| --- | --- |
| `--gradient-csv` | Gradient CSV (required). |
| `--left-surface` / `--right-surface` | Surfaces (required). |
| `--left-annot` / `--right-annot` | Atlas annotations (required). |
| `--output` | PNG to write (required). |
| `--gradient-columns` | Columns to render (default `Gradient_1` .. `Gradient_4`). |
| `--label-text` | Row labels (default `Grad1`, `Grad2`, ...). |
| `--name-map` | JSON of name overrides. |
| `--cmap` | Colormap (default `viridis_r`). |
| `--nan-color R G B` | Colour (0-1) for vertices without a value (default light grey). |
| `--color-range MIN MAX` | One fixed colour scale for all rows (default: each row autoscaled). |
| `--zoom`, `--size` | Camera zoom and image size in pixels. |

## Example

```bash
/path/to/vtk-python plot_group_gradients_brain_atlas.py \
    --gradient-csv group_gradients_rescaled_minus1_1.csv \
    --left-surface lh.surf.gii --right-surface rh.surf.gii \
    --left-annot lh.myatlas.annot --right-annot rh.myatlas.annot \
    --color-range -1 1 \
    --output group_gradients_brain.png
```

## Checked against

- Desikan-Killiany on fsaverage5, annotation rebuilt from ENIGMA's own
  vertex labels: the values painted on all 20484 vertices are identical to
  ENIGMA's `parcel_to_surface` output.
- Schaefer-200 on fsaverage5 with synthetic values (each parcel's mean x/y/z
  coordinate): smooth anterior-posterior and inferior-superior fields in the
  expected directions.

It has not been run against an official Brainnetome annotation; name
matching for that atlas is expected to work but is untested.
