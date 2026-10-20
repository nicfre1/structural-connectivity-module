# plot_group_gradients_brain.py

Render the first group gradients on a cortical surface. One script, three
modes (`dk` is the default when no mode is given):

| Mode | For | You provide |
| --- | --- | --- |
| `dk` | Desikan-Killiany regions on fsaverage5 (ENIGMA's `aparc_fsa5`) | the gradient CSV only |
| `atlas` | region-wise gradients for **any atlas** | the gradient CSV, a surface and an atlas annotation per hemisphere |
| `vertex` | vertex-wise gradients | the gradients `.npz` and the group midthickness surface per hemisphere |

Through the entry point these are `gradients.py plot`, `plot-atlas` and
`plot-vertex`.

## Requirements

A VTK stack (`brainspace[plotting]`, `nibabel`, `vtk`, and `enigmatoolbox` for
`dk`; see `../requirements-brain.txt`); run with a VTK-capable interpreter. On
some interpreters `import vtk` crashes the process, so the covariance scripts
call this one as a subprocess and fall back to a matplotlib figure if it fails.
A dedicated environment (Python 3.10, VTK 9.3.x) is recommended.

The matching logic of the `atlas` mode needs only numpy/pandas and is importable
without VTK (that is what the unit tests exercise); nibabel, ENIGMA and VTK are
imported only when a mode runs.

## Common options

| Option | Meaning |
| --- | --- |
| `--output` | PNG to write (required). |
| `--label-text` | Row labels (default `Grad1`, `Grad2`, ...). |
| `--cmap` | Colormap (default `viridis_r`). |
| `--zoom`, `--size` | Camera zoom and image size in pixels. |

## `dk`

Data-free port of the validated `render_dk68_myrna_style.py`: each of the 68
Desikan-Killiany values is painted on fsaverage5 with ENIGMA Toolbox's
`parcel_to_surface(..., "aparc_fsa5")` and BrainSpace's `plot_hemispheres`.

Input: a gradient CSV with columns `Gradient_1` .. `Gradient_4`, in the
canonical Desikan-Killiany order (lh 34 regions, then rh 34), ideally rescaled
to `[-1, 1]` (`group_gradients_rescaled_minus1_1.csv`).

| Option | Meaning |
| --- | --- |
| `--gradient-csv` | Gradient CSV (required). |
| `--gradient-columns` | Columns to render (default `Gradient_1` .. `Gradient_4`). |

```bash
/path/to/vtk-python plot_group_gradients_brain.py \
    --gradient-csv results/group_covariance/group_gradients_rescaled_minus1_1.csv \
    --output results/group_covariance/group_gradients_brain.png
```

## `atlas`

Atlas-agnostic: the vertex -> region assignment comes from an annotation file
you provide, so Brainnetome, Schaefer, Glasser, a custom atlas, ... all work
the same way, on whichever surface the annotation was defined on.

| Input | Option | Notes |
| --- | --- | --- |
| Gradient values | `--gradient-csv` | First column = region names, then `Gradient_1` ... |
| Surfaces | `--left-surface`, `--right-surface` | GIFTI `.surf.gii` / `.gii`, or a FreeSurfer surface file. |
| Atlas on the surface | `--left-annot`, `--right-annot` | FreeSurfer `.annot`, or a GIFTI label file; one entry per vertex of the matching surface. |

| Option | Meaning |
| --- | --- |
| `--gradient-columns` | Columns to render (default `Gradient_1` .. `Gradient_4`). |
| `--name-map` | JSON `{"csv region name": "annotation label name"}` for names the automatic matching cannot reconcile. |
| `--nan-color R G B` | Colour (0-1) for vertices without a value (default light grey). |
| `--color-range MIN MAX` | One fixed colour scale for all rows (default: each row autoscaled). |

**How region names are matched.** CSV region names and annotation label names
are compared after normalisation: lower-cased, split on `_ - .` and spaces,
the `ctx` token dropped, and a hemisphere token (`lh`, `rh`, `l`, `r`, `left`,
`right`) removed and remembered. So `ctx-lh-bankssts`, `lh_bankssts` and
`bankssts` all match the label `bankssts`; `Sch200_LH_001` matches
`Sch200_LH_001`; `SFG_L_6_1` matches `SFG_6_1`.

- A CSV row with a hemisphere token is drawn only on that hemisphere.
- A CSV row without one (for example values already averaged over left and
  right) is drawn on both hemispheres.
- Labels named `unknown`, `medialwall`, `corpuscallosum`, ... , or with no
  match, are left uncoloured (`--nan-color`).

For each hemisphere the script prints how many annotation labels were matched
and lists unmatched labels and CSV regions that were not drawn. If nothing
matches it stops with an explicit message; if the annotation and the surface
have a different number of vertices it stops too.

```bash
/path/to/vtk-python plot_group_gradients_brain.py atlas \
    --gradient-csv group_gradients_rescaled_minus1_1.csv \
    --left-surface lh.surf.gii --right-surface rh.surf.gii \
    --left-annot lh.myatlas.annot --right-annot rh.myatlas.annot \
    --color-range -1 1 --output group_gradients_brain.png
```

Checked against: Desikan-Killiany on fsaverage5 with an annotation rebuilt from
ENIGMA's own vertex labels (the values painted on all 20484 vertices are
identical to ENIGMA's `parcel_to_surface` output), and Schaefer-200 with
synthetic values (each parcel's mean x/y/z coordinate), which give smooth fields
in the expected directions. It has not been run against an official Brainnetome
annotation; name matching for that atlas is expected to work but is untested.

## `vertex`

Paints vertex-wise gradients (`group_gradients.npz` written by
`compute_covariance_gradients_vertex.py`: `gradients`/`gradients_rescaled` on the
full `2 * n_per_hemi` grid, NaN outside the cortical mask, plus `n_per_hemi`)
on the group midthickness surface you provide. There is no bundled generic
vertex mesh to fall back on: the group midthickness is built from real subject
anatomy by `build_cohort_thickness_vertex.py`.

| Option | Meaning |
| --- | --- |
| `--gradients-npz` | Gradients `.npz` (required). |
| `--left-surface` / `--right-surface` | The group midthickness `.surf.gii` the gradients were computed on (same vertex count and order as `n_per_hemi`). |
| `--n-gradients` | Number of gradients to render (default 4). |

```bash
/path/to/vtk-python plot_group_gradients_brain.py vertex \
    --gradients-npz results/group_covariance/group_gradients.npz \
    --left-surface groupmid_hemi-L.surf.gii --right-surface groupmid_hemi-R.surf.gii \
    --output results/group_covariance/group_gradients_brain.png
```
