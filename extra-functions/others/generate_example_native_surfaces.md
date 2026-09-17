# generate_example_native_surfaces.py

Generate a fully synthetic "template" + cohort "derivatives" tree, in the
exact file layout `build_vertex_grid_assets.py`,
`resample_subject_thickness_vertex.py` and `build_cohort_thickness_vertex.py`
expect, so the whole native-surface -> low-resolution-grid resampling chain
can be smoke-tested with Connectome Workbench and no real anatomy at all.

The surfaces are plain icospheres (`wb_command -surface-create-sphere`) and
the metrics are random values with the right shape: geometrically valid
inputs for `wb_command`, but not anatomically meaningful. This proves the
resampling *plumbing* is correct end to end (each `wb_command` call, each
file role, each hand-off between the three scripts); it does not produce a
scientifically meaningful result. Point the same three scripts at a real
dHCP symmetric template and real subject derivatives to get one.

## Requirements

Connectome Workbench (`wb_command`) only -- no template, no subject data.

## What it writes (under `--out_dir`)

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

The "registration" sphere is a copy of the synthetic native mesh itself
(already a valid sphere with the right vertex count) -- a stand-in that
satisfies `wb_command`'s geometric requirements without claiming any real
registration.

## Options

| Option | Meaning |
| --- | --- |
| `--out_dir` | Directory where `template/` and `derivatives/` are written. |
| `--wb_command` | Path to (or name of) `wb_command` (required). |
| `--n_subjects` | Number of synthetic subjects (default 3). |
| `--n_native` | Vertices per hemisphere on the synthetic "native" mesh (default 42). |
| `--seed` | Random seed for reproducibility. |

## Example (full smoke test of the resampling chain)

```bash
WB=/Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command

python generate_example_native_surfaces.py --out_dir example_native --wb_command "$WB" --n_subjects 3

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
