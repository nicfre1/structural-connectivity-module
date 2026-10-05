# extra-functions

Helper scripts to move cortical thickness statistics from the **DrawEM32**
parcellation to the **MCRIBS** (M-CRIB-S) Desikan-Killiany parcellation with
**68 cortical regions** (34 per hemisphere), then run the validated structural
covariance and gradient analysis on them.

No study data is shipped here. Every script is data-agnostic: anyone can run the
chain on their own thickness files and obtain their own results. Synthetic files
are provided only so the pipeline can be exercised end to end.

Tested with Python 3.11 and the packages in `../requirements.txt`. The cortical
surface render of the group gradients additionally needs a working VTK stack
(`enigmatoolbox`, `brainspace[plotting]`, `nibabel`, `vtk`); without it the
pipeline still completes and writes a matplotlib fallback figure.

## Contents

| File | Purpose |
| --- | --- |
| `drawem32_to_mcribs.py` | Convert DrawEM32 thickness TSVs to the MCRIBS 68-region convention, driven by an external mapping file. |
| `mapping_drawem32_to_mcribs.json` | Mapping **template**: the 34 Desikan-Killiany target regions with empty source lists to fill in. |
| `compute_covariance_gradients.py` | Individualized covariance matrices, group covariance matrix, group + Procrustes-aligned individual gradients, figures. Data-free port of the validated reference scripts. |
| `plot_group_gradients_brain.py` | Standalone cortical-surface render of the 4 group gradients (aparc_fsa5). Run with a VTK-capable interpreter. |
| `generate_example_thickness.py` | Write synthetic DrawEM32 thickness files so the chain can be run without real data. |
| `mapping_example.json` | Example mapping matching the synthetic files (demo only, not anatomically valid). |
| `compute_covariance_gradients_vertex.py` | Vertex-wise counterpart of `compute_covariance_gradients.py`: same equations, one row per cortical vertex of a common surface grid instead of one row per DK region. |
| `plot_group_gradients_brain_vertex.py` | Cortical-surface render of the vertex-wise group gradients, on a caller-supplied group midthickness mesh. Run with a VTK-capable interpreter. |
| `generate_example_thickness_vertex.py` | Write a synthetic cohort vertex-thickness file so the vertex-wise chain can be run without real surface data. |
| `vertex_resample_utils.py` | Shared geometry helpers: native surface -> common template mesh -> low-resolution grid, via Connectome Workbench. |
| `build_vertex_grid_assets.py` | Build the low-resolution icosphere grid from a dHCP symmetric template. Needs `wb_command` + your own copy of the template. |
| `resample_subject_thickness_vertex.py` | Resample one subject's native thickness onto the low-resolution grid. Needs `wb_command` + your own subject derivatives. |
| `build_cohort_thickness_vertex.py` | Loop the above over a cohort, build the group midthickness surface, smooth, define the group mask, write the cohort file `compute_covariance_gradients_vertex.py` consumes. |
| `generate_example_native_surfaces.py` | Write a fully synthetic template + subject derivatives tree so the whole resampling chain can be smoke-tested with `wb_command` and no real anatomy. |
| `others/` | Detailed notes for each script (kept out of the code headers). |

The method (cohort z-score, `np.corrcoef` group matrix, `GradientMaps` +
Procrustes, rescale to `[-1, 1]`) is summarised in [`others/README.md`](others/README.md)
and unchanged from the validated reference scripts.

## 1. Provide the DrawEM32 -> MCRIBS correspondence

Edit `mapping_drawem32_to_mcribs.json`. For every MCRIBS/Desikan-Killiany region,
list the DrawEM32 source labels that overlap it (column name minus the
`lh_`/`rh_` prefix and the `_thickness` suffix). Same list for both hemispheres.
Validate it against your atlas documentation before trusting the results.

## 2. Convert

```bash
python drawem32_to_mcribs.py \
    lh_drawem32_thickness.tsv rh_drawem32_thickness.tsv \
    --mapping mapping_drawem32_to_mcribs.json \
    --out_dir mcribs_thickness_68
```

Produces `thickness_68_lh_stats.tsv`, `thickness_68_rh_stats.tsv` and a
`conversion_report_{lh,rh}.tsv`.

## 3. Covariance and gradients

```bash
python compute_covariance_gradients.py \
    mcribs_thickness_68/thickness_68_lh_stats.tsv \
    mcribs_thickness_68/thickness_68_rh_stats.tsv \
    --out_dir mcribs_thickness_68/results
```

Or, if you already have one `<subject>_thickness.csv` per subject (columns
`Region`, `ThickAvg`, 68 Desikan-Killiany rows):

```bash
python compute_covariance_gradients.py --thickness_dir dk68_stats --out_dir results
```

Outputs:

```
individualized_covariance/csv/<subject>_individualized_covariance.csv
individualized_covariance/png/<subject>_individualized_covariance.png
group_covariance/group_covariance_matrix.csv
group_covariance/group_gradients.csv
group_covariance/group_gradients_rescaled_minus1_1.csv
group_covariance/group_lambdas.csv
group_covariance/group_lambdas.png
group_covariance/group_gradients_brain.png
individual_gradients4/<subject>_gradients_aligned.csv
individual_gradients4/<subject>_gradients_aligned_rescaled_minus1_1.csv
```

`group_gradients_brain.png` is the cortical-surface render when a VTK stack is
available, otherwise a matplotlib fallback. To produce (or re-produce) the real
render with a dedicated interpreter:

```bash
/path/to/vtk-python plot_group_gradients_brain.py \
    --gradient-csv results/group_covariance/group_gradients_rescaled_minus1_1.csv \
    --output results/group_covariance/group_gradients_brain.png
```

## Run the full chain on synthetic data

```bash
python generate_example_thickness.py --out_dir example_data --n_subjects 30
python drawem32_to_mcribs.py \
    example_data/drawem32_thickness_lh_stats.tsv \
    example_data/drawem32_thickness_rh_stats.tsv \
    --mapping mapping_example.json --out_dir example_out
python compute_covariance_gradients.py \
    example_out/thickness_68_lh_stats.tsv \
    example_out/thickness_68_rh_stats.tsv \
    --out_dir example_out/results
```

## Vertex-wise branch

For a vertex-wise analysis (a common surface grid instead of the 68 DK
regions), `compute_covariance_gradients_vertex.py` runs the exact same
covariance/gradient equations one vertex at a time. It starts from an
already-built cohort file (`thickness`, `mask`, `hemi`, `subject_ids`,
`n_per_hemi` — see
[`others/compute_covariance_gradients_vertex.md`](others/compute_covariance_gradients_vertex.md)):

```bash
python generate_example_thickness_vertex.py --out_dir example_data_vertex --n_subjects 25
python compute_covariance_gradients_vertex.py \
    example_data_vertex/cohort_thickness_vertex_dens-demo.npz \
    --out_dir example_data_vertex/results
```

The cortical-surface render (`plot_group_gradients_brain_vertex.py`) needs a
group midthickness mesh (`--left_surface`/`--right_surface`) supplied by the
caller; there is no bundled generic vertex mesh to fall back on, so it is
skipped unless those are given.

### Building the cohort file from native surfaces

Going from native cortical surfaces to that common grid is itself ported,
data-free: `build_vertex_grid_assets.py` + `resample_subject_thickness_vertex.py`
+ `build_cohort_thickness_vertex.py` (see
[`others/README.md`](others/README.md#vertex-wise-branch)). It needs two
things this repo does not ship: Connectome Workbench (`wb_command`, a free
external tool) and a local copy of the dHCP symmetric surface template (a
public, group-average atlas, not subject data). With those and your own
subject derivatives (standard `anat/`+`xfm/` layout):

```bash
WB=/path/to/wb_command

python build_vertex_grid_assets.py \
    --template_dir /path/to/dhcp-symmetric-template --wb_command "$WB" \
    --output_dir assets --n_vertices 10242

python build_cohort_thickness_vertex.py \
    --derivatives_root /path/to/derivatives --output_dir cohort \
    --assets_dir assets --template_dir /path/to/dhcp-symmetric-template \
    --wb_command "$WB" --lowres 10k --fwhm 10 --manifest my_cohort.csv

python compute_covariance_gradients_vertex.py \
    cohort/cohort_thickness_vertex_dens-10k.npz --out_dir cohort/results \
    --left_surface assets/groupmid_hemi-L_dhcpSym_dens-10k.surf.gii \
    --right_surface assets/groupmid_hemi-R_dhcpSym_dens-10k.surf.gii
```

Without either the template or real subject data,
`generate_example_native_surfaces.py` fabricates a stand-in for both so the
whole chain can still be run and smoke-tested end to end:

```bash
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
