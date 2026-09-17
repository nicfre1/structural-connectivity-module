# compute_covariance_gradients_vertex.py

Vertex-wise counterpart of `compute_covariance_gradients.py`: individualized
(cohort-referenced) covariance, group covariance matrix, group gradients and
Procrustes-aligned individual gradients, one row per cortical vertex of a
common surface grid instead of one row per Desikan-Killiany region.

Data-free port of the vertex branch of the validated reference pipeline
(`etape3_covariance_vertex_individualisee.py` +
`etape4_gradients_vertex_groupe_et_individuels.py`). Every equation is
unchanged from the region-level script and is in fact reused from it
(`individualized_covariance_from_z`, `rescale_to_minus1_1`):

```
cohort z-score per vertex:   X_z = (X - X.mean(0)) / X.std(0)
individual matrix:           M[i, j] = 1 / exp((z_i - z_j) ** 2)
group covariance matrix:     np.corrcoef(X_z, rowvar=False)
group gradients:             GradientMaps(n_components=10, random_state=0,
                                          approach="dm", kernel=None)
                             .fit(group_matrix, sparsity=0.9)
individual gradients:        GradientMaps(..., alignment="procrustes")
                             .fit(individual_matrix, reference=group_grad,
                                  sparsity=0.9)   # one subject at a time
gradient rescaling:          2 * (x - min) / (max - min) - 1   (first 4)
```

Individual gradients are fit one subject at a time (instead of passing the
whole list to a single `.fit()`, as the region-level script does): at
10k-20k vertices per hemisphere, keeping every subject's `V x V` matrix in
memory at once does not fit. Procrustes alignment is independent per
subject, so the result is identical either way.

## What is out of scope here

Getting from a subject's native cortical surface to a shared low-resolution
vertex grid (resampling with Connectome Workbench, area correction, a
template registration sphere, surface smoothing, building the group cortical
mask) needs real per-subject anatomical surfaces and `wb_command`. That
resampling is not data-free and is not ported here — this script starts
where the region-level one starts: from an already-built cohort file.

## Input

A cohort `.npz` with:

| Key | Shape | Meaning |
| --- | --- | --- |
| `thickness` | `(n_subjects, n_valid)` | thickness at every valid (masked) vertex |
| `mask` | `(2 * n_per_hemi,)` bool | which vertices of the full L+R grid are valid cortex |
| `hemi` | `(2 * n_per_hemi,)` | `'L'`/`'R'` per vertex of the full grid |
| `subject_ids` | `(n_subjects,)` | subject identifiers |
| `n_per_hemi` | scalar | vertices per hemisphere on the full grid |

This is the schema written by the real vertex pipeline's cohort-assembly
step; `generate_example_thickness_vertex.py` writes a synthetic one with the
same keys.

## Options

| Option | Meaning |
| --- | --- |
| `--out_dir` | Directory where CSV/NPZ files and figures are written. |
| `--n_components` | Number of gradient components (default 10). |
| `--sparsity` | Sparsity passed to `GradientMaps` (default 0.9). |
| `--save_matrices` | Also write each subject's `V x V` matrix as `.npy` (large). |
| `--figure_max_vertices` | Subsampling of the individualized-covariance heatmaps. |
| `--left_surface` / `--right_surface` | Group midthickness `.surf.gii`, for the cortical-surface render. |
| `--no_brain_plot` | Skip the render even if surfaces are given. |
| `-v` | Verbosity level. Default `WARNING`, or `INFO` when `-v` is given. |

## Outputs (in `--out_dir`)

```
individualized_covariance/png/<subject>_individualized_covariance.png
individualized_covariance/npy/<subject>_individualized_covariance.npy   (--save_matrices only)
group_covariance/group_gradients.npz                 (full grid, NaN outside mask)
group_covariance/group_gradients_valid_vertices.csv  (valid vertices only)
group_covariance/group_lambdas.csv
group_covariance/group_lambdas.png
group_covariance/group_gradients_brain.png           (only with --left_surface/--right_surface)
individual_gradients/<subject>_gradients_aligned.npz
individual_gradients/<subject>_gradients_aligned_valid_vertices.csv
```

No `.csv` is written for the full grid (it would carry a `NaN` row per
excluded vertex); use the `.npz` (`gradients`/`gradients_rescaled` + `mask`
+ `n_per_hemi`) whenever the full-grid, surface-alignable layout is needed,
e.g. for `plot_group_gradients_brain_vertex.py`.

## Example

```bash
python generate_example_thickness_vertex.py --out_dir example_data_vertex --n_subjects 25
python compute_covariance_gradients_vertex.py \
    example_data_vertex/cohort_thickness_vertex_dens-demo.npz \
    --out_dir example_data_vertex/results
```
