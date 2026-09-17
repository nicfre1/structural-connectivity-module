# plot_group_gradients_brain_vertex.py

Render the first four vertex-wise group gradients on a group cortical
surface.

Vertex-wise counterpart of `plot_group_gradients_brain.py`: instead of
ENIGMA Toolbox's bundled fsa5 (parcel-based) surface and
`parcel_to_surface`, values are painted directly on a group midthickness
mesh supplied by the caller (`--left-surface`/`--right-surface`), because
there is no bundled generic vertex mesh to fall back on -- the group
midthickness is itself built from real subject anatomy by the vertex
pipeline's cohort-assembly step.

Data-free port of `etape5_rendu_vertex_dhcpsym.py`: same BrainSpace call
(`plot_hemispheres`, colormap `viridis_r`, colorbar, four rows `Grad1` ..
`Grad4`).

## Requirements

Needs a working VTK stack: `brainspace[plotting]`, `nibabel`, `vtk` (see
`../requirements-brain.txt`). On some interpreters `import vtk` crashes the
process, so this lives in its own script and
`compute_covariance_gradients_vertex.py` calls it as a subprocess. A
dedicated environment (Python 3.10, VTK 9.3.x) is recommended.

## Input

- `--gradients-npz`: a gradients `.npz` as written by
  `compute_covariance_gradients_vertex.py` (`gradients`/`gradients_rescaled`
  on the full `2 * n_per_hemi` grid, NaN outside the cortical mask, +
  `n_per_hemi`).
- `--left-surface` / `--right-surface`: the group midthickness
  `.surf.gii` the gradients were computed on (same vertex count and order
  as `n_per_hemi`).

## Options

| Option | Meaning |
| --- | --- |
| `--gradients-npz` | Path to the gradients `.npz` (required). |
| `--left-surface` / `--right-surface` | Group midthickness surfaces (required). |
| `--output` | Path of the PNG to write (required). |
| `--n-gradients` | Number of gradients to render (default 4). |
| `--label-text` | Row labels (default `Grad1 Grad2 Grad3 Grad4`). |
| `--zoom` | BrainSpace camera zoom (default 1.55). |
| `--size` | Render size in pixels (default `1200 400`). |
| `--cmap` | Colormap (default `viridis_r`). |

## Example

```bash
/path/to/vtk-python plot_group_gradients_brain_vertex.py \
    --gradients-npz results/group_covariance/group_gradients.npz \
    --left-surface groupmid_hemi-L.surf.gii \
    --right-surface groupmid_hemi-R.surf.gii \
    --output results/group_covariance/group_gradients_brain.png
```
