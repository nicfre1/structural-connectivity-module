# generate_example_data.py

Generate synthetic example data, so each pipeline can be run end to end
without any real data. One script, three kinds (through the entry point:
`gradients.py example`, `example-vertex`, `example-native`).

In all three kinds, thickness is built from two smooth latent spatial axes (a
dominant one and a weaker one) and one score per axis per subject, so regions
covary across subjects along those axes: the cohort covariance matrix has real
structure and the first group gradient recovers the dominant axis. None of it
is anatomically meaningful.

## `drawem32`

Synthetic DrawEM32 cortical thickness TSVs, so the DrawEM32 -> MCRIBS ->
covariance/gradients chain can be run without subject data. 40 cortical parcels
per hemisphere (`dwm_01` ... `dwm_40`), which is what `mapping_example.json`
expects; each subject's scores are shared by both hemispheres.

| Option | Meaning |
| --- | --- |
| `--out_dir` | Where the synthetic TSV files are written (default `example_data`). |
| `--n_subjects` | Number of synthetic subjects (default 25). |
| `--id_column` | Subject identifier column (default `sample`). |
| `--seed` | Random seed (default 0). |
| `--with_area` | Also write matching surface-area files (for `area_weighted`). |

```bash
python generate_example_data.py drawem32 --out_dir example_data --n_subjects 30
```

## `vertex`

A synthetic cohort vertex-thickness `.npz`, same schema as the real vertex
pipeline's cohort file: `thickness` `(n_subjects, n_valid)`, `mask`/`hemi` on
the full `2 * n_per_hemi` grid, `subject_ids`, `n_per_hemi`. A
`medial_wall_fraction` of low-index vertices per hemisphere is marked invalid
(`mask=False`), standing in for the medial wall / group cortical mask.

| Option | Meaning |
| --- | --- |
| `--out_dir` | Where the cohort `.npz` is written (default `example_data_vertex`). |
| `--n_subjects` | Number of synthetic subjects (default 25). |
| `--n_per_hemi` | Vertices per hemisphere on the synthetic grid (default 162). |
| `--medial_wall_fraction` | Fraction of vertices per hemisphere excluded as non-cortical (default 0.1). |
| `--seed` | Random seed (default 0). |

```bash
python generate_example_data.py vertex --out_dir example_data_vertex --n_subjects 25
```

## `native`

A fully synthetic "template" + cohort "derivatives" tree, in the exact file
layout `build_vertex_grid_assets.py`, `resample_subject_thickness_vertex.py` and
`build_cohort_thickness_vertex.py` expect, so the whole native-surface ->
low-resolution-grid resampling chain can be smoke-tested with Connectome
Workbench and no real anatomy. Surfaces are plain icospheres
(`wb_command -surface-create-sphere`) and metrics are random values with the
right shape: geometrically valid inputs for `wb_command`, not anatomically
meaningful. It proves the resampling *plumbing* end to end; point the same
scripts at a real dHCP symmetric template and real subject derivatives for a
real result.

The "registration" sphere is a copy of the synthetic native mesh itself
(already a valid sphere with the right vertex count), a stand-in that satisfies
`wb_command`'s geometric requirements without claiming any real registration.

Needs Connectome Workbench (`wb_command`) only; no template, no subject data.

```
template/
  week-40_hemi-{left,right}_space-dhcpSym_dens-32k_sphere.surf.gii
  week-40_hemi-LR_space-dhcpSym_dens-32k_surface_area.shape.gii
  week-40_hemi-LR_space-dhcpSym_dens-32k_sulc.shape.gii
derivatives/
  sub-EX{001,002,...}/ses-1/
    anat/*_hemi-{left,right}_{thickness,desc-medialwall_mask}.shape.gii
    anat/*_hemi-{left,right}_midthickness.surf.gii
    xfm/*_hemi-{left,right}_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii
```

| Option | Meaning |
| --- | --- |
| `--out_dir` | Where `template/` and `derivatives/` are written (default `example_native`). |
| `--wb_command` | Path to (or name of) `wb_command` (required). |
| `--n_subjects` | Number of synthetic subjects (default 3). |
| `--n_native` | Vertices per hemisphere on the synthetic native mesh (default 42). |
| `--seed` | Random seed (default 0). |

Full smoke test of the resampling chain:

```bash
WB=/Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command

python generate_example_data.py native --out_dir example_native --wb_command "$WB" --n_subjects 3

python build_vertex_grid_assets.py \
    --template_dir example_native/template --wb_command "$WB" \
    --output_dir example_native/assets --n_vertices 642

python build_cohort_thickness_vertex.py \
    --derivatives_root example_native/derivatives --output_dir example_native/cohort \
    --assets_dir example_native/assets --template_dir example_native/template \
    --wb_command "$WB" --lowres 0.6k --fwhm 5 --n_subjects 3

python compute_covariance_gradients_vertex.py \
    example_native/cohort/cohort_thickness_vertex_dens-0.6k.npz \
    --out_dir example_native/results --no_brain_plot
```
