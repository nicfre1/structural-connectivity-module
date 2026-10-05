# build_cohort_thickness_vertex.py

Vertex-wise cortical thickness for a whole COHORT: loops
`resample_subject_thickness_vertex.py` over a batch of subjects, builds a
group-average midthickness surface, smooths, defines the group cortical
mask, and writes the cohort `.npz` that
`compute_covariance_gradients_vertex.py` consumes directly.

Data-free port of `etape2_thickness_vertex_cohorte.py`.

## On top of the single-subject step

1. **Group-average midthickness surface**, in the common template space:
   each subject's native midthickness is resampled onto the 32k mesh, then
   averaged, then downsampled to the low-resolution grid. This is the
   geometric support for smoothing and for cortical-surface renders
   (`plot_group_gradients_brain_vertex.py`).
2. **Surface smoothing** of each subject's thickness on that group surface,
   Gaussian kernel, `--fwhm` mm (default 10, standard in morphometry;
   needed at the vertex scale to control single-vertex noise -- there is no
   equivalent at the region level, where the parcel average already plays
   that role).
3. **Group cortical mask** = vertices valid in at least `--min_coverage`
   (default 1.0 = all) of the subjects, excluding the medial wall and any
   vertex invalid in even one subject.

## Requirements

Connectome Workbench (`wb_command`) and subject session directories in the
standard `anat/` + `xfm/` layout (see
`others/resample_subject_thickness_vertex.md`). See
`generate_example_native_surfaces.py` for a synthetic stand-in that
exercises this whole script without real data.

## Subject selection

Either `--manifest` (a CSV with columns `participant_id[,session_id]` --
your own cohort list; nothing here ships one) or the first `--n_subjects`
sessions found under `--derivatives_root`, plus any `--must_include`.

## Options

| Option | Meaning |
| --- | --- |
| `--derivatives_root` | Root directory of subject derivatives (`sub-*/ses-*/{anat,xfm}`). |
| `--output_dir` | Where the cohort file and working files are written. |
| `--assets_dir` | Directory from `build_vertex_grid_assets.py`. |
| `--template_dir` | dHCP symmetric template directory. |
| `--wb_command` | Path to (or name of) `wb_command`. |
| `--lowres` | Grid label, matching `build_vertex_grid_assets.py`. |
| `--fwhm` | Smoothing kernel, mm (default 10). |
| `--manifest` | CSV of `participant_id[,session_id]` to use exactly. |
| `--n_subjects` / `--must_include` | Selection when no manifest is given. |
| `--min_coverage` | Fraction of subjects a vertex must be valid in (default 1.0). |
| `--thickness_desc` | `thickness` (default) or e.g. `desc-corr_thickness`. |

## Output

`<output_dir>/cohort_thickness_vertex_dens-{lowres}.npz`: `thickness`
`(n_subjects, V)`, `mask`/`hemi` `(2*n_per_hemi,)`, `subject_ids`,
`fwhm`, `lowres`, `min_coverage`, `n_per_hemi` -- the schema
`compute_covariance_gradients_vertex.py` expects (see
`others/compute_covariance_gradients_vertex.md`).

## Example

```bash
python build_cohort_thickness_vertex.py \
    --derivatives_root /path/to/derivatives --output_dir ./cohort \
    --assets_dir ./assets --template_dir /path/to/dhcp-symmetric-template \
    --wb_command /Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command \
    --lowres 10k --fwhm 10 --manifest my_cohort.csv
```
