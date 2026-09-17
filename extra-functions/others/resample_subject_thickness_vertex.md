# resample_subject_thickness_vertex.py

Vertex-wise cortical thickness for ONE subject, native surface -> the common
low-resolution grid built by `build_vertex_grid_assets.py`.

Data-free port of `etape1_thickness_vertex_par_sujet.py`. Uses
`vertex_resample_utils.resample_subject` for the actual chain (see
`others/vertex_resample_utils.md`); no smoothing here -- that happens once,
on the group-average surface, in `build_cohort_thickness_vertex.py`.

## Requirements

Connectome Workbench (`wb_command`) and a subject session directory in the
standard `anat/` + `xfm/` layout (native thickness, medial-wall mask,
midthickness, and the native-to-common-space registration sphere). See
`generate_example_native_surfaces.py` for a synthetic stand-in.

## Options

| Option | Meaning |
| --- | --- |
| `session_dir` | `.../sub-XXX/ses-YYY` (must contain `anat/` and `xfm/`). |
| `--assets_dir` | Directory from `build_vertex_grid_assets.py` (default `assets`). |
| `--template_dir` | dHCP symmetric template directory. |
| `--wb_command` | Path to (or name of) `wb_command`. |
| `--lowres` | Grid label, matching `build_vertex_grid_assets.py` (default `10k`). |
| `--thickness_desc` | `thickness` (default) or e.g. `desc-corr_thickness`. |
| `--output` | Output `.npz` path. |

## Output

`<output>.npz`: `thickness` (2N,), `roi` (2N,) bool, `n_per_hemi`,
`subject`, `lowres`; L then R. The intermediate `.shape.gii` files are kept
under `<output's parent>/_work/<subject>/` for reuse by
`build_cohort_thickness_vertex.py`'s smoothing step.

## Example

```bash
python resample_subject_thickness_vertex.py \
    /path/to/derivatives/sub-XXX/ses-YYY \
    --assets_dir ./assets --template_dir /path/to/dhcp-symmetric-template \
    --wb_command /Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command \
    --lowres 10k --output vertex_thickness/sub-XXX_ses-YYY_thickness_vertex.npz
```
