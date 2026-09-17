# others

Detailed notes for the scripts in `extra-functions/`. The scripts themselves
keep a short header (description + example), same style as the top-level
`compute_covariance.py`; the full explanation of each one lives here.

| Script | Notes |
| --- | --- |
| `drawem32_to_mcribs.py` | [drawem32_to_mcribs.md](drawem32_to_mcribs.md) |
| `compute_covariance_gradients.py` | [compute_covariance_gradients.md](compute_covariance_gradients.md) |
| `plot_group_gradients_brain.py` | [plot_group_gradients_brain.md](plot_group_gradients_brain.md) |
| `generate_example_thickness.py` | [generate_example_thickness.md](generate_example_thickness.md) |
| `compute_covariance_gradients_vertex.py` | [compute_covariance_gradients_vertex.md](compute_covariance_gradients_vertex.md) |
| `plot_group_gradients_brain_vertex.py` | [plot_group_gradients_brain_vertex.md](plot_group_gradients_brain_vertex.md) |
| `generate_example_thickness_vertex.py` | [generate_example_thickness_vertex.md](generate_example_thickness_vertex.md) |

## Method summary

For a cohort of `S` subjects and `R = 68` Desikan-Killiany regions (lh then rh,
no left/right averaging):

1. `X` is the `S x R` matrix of regional mean thickness.
2. Cohort z-score, per region across subjects: `X_z = (X - X.mean(0)) / X.std(0)`.
3. Individualized (cohort-referenced) covariance, one per subject:
   `M[i, j] = 1 / exp((z_i - z_j) ** 2)`.
4. Group covariance matrix: `np.corrcoef(X_z, rowvar=False)`.
5. Group gradients: `GradientMaps(n_components=10, random_state=0).fit(group_matrix)`.
6. Individual gradients: `GradientMaps(..., alignment="procrustes").fit(individual_covariances, reference=group_gradients)`.
7. First four gradients rescaled to `[-1, 1]` via `2*(x-min)/(max-min) - 1`.

## Vertex-wise branch

`compute_covariance_gradients_vertex.py` runs the exact same seven steps at
the vertex level (a common surface grid, ~10k-20k vertices, instead of the
68 regions above); see
[compute_covariance_gradients_vertex.md](compute_covariance_gradients_vertex.md).
Getting from native cortical surfaces to that common grid needs real
per-subject anatomy and Connectome Workbench, and is out of scope here — the
vertex script starts from an already-built cohort file, real or synthetic
(`generate_example_thickness_vertex.py`).
