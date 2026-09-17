# generate_example_thickness_vertex.py

Generate a synthetic cohort vertex-thickness `.npz`, so
`compute_covariance_gradients_vertex.py` can be run end to end without any
real surface data.

Same schema as the real vertex pipeline's cohort file: `thickness`
`(n_subjects, n_valid)`, `mask`/`hemi` on the full `2 * n_per_hemi` grid,
`subject_ids`, `n_per_hemi`. A `medial_wall_fraction` of low-index vertices
per hemisphere is marked invalid (`mask=False`), standing in for the medial
wall / group cortical mask of the real pipeline.

As in `generate_example_thickness.py`, thickness is built from two smooth
latent spatial axes plus a per-subject score, so the cohort covariance has
real structure and the first group gradient recovers the dominant axis.

## Options

| Option | Meaning |
| --- | --- |
| `--out_dir` | Directory where the synthetic cohort `.npz` is written. |
| `--n_subjects` | Number of synthetic subjects (default 25). |
| `--n_per_hemi` | Vertices per hemisphere on the synthetic grid (default 162). |
| `--medial_wall_fraction` | Fraction of vertices per hemisphere excluded as non-cortical (default 0.1). |
| `--seed` | Random seed for reproducibility (default 0). |

## Example

```bash
python generate_example_thickness_vertex.py --out_dir example_data_vertex --n_subjects 25
```
