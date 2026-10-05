# vertex_resample_utils.py

Shared geometry helpers for the native-surface -> low-resolution-grid
resampling chain, used by `build_vertex_grid_assets.py`,
`resample_subject_thickness_vertex.py` and `build_cohort_thickness_vertex.py`
so a single subject and a whole cohort go through exactly the same steps.

No gradient equations and no study data here -- purely geometry, via
Connectome Workbench (`wb_command`).

## Resampling chain (per hemisphere)

```
native thickness (subject's own native mesh)
  |  wb_command -metric-resample  (ADAP_BARY_AREA, area-corrected)
  |  sphere: subject's own native-to-common-space registration
  v
thickness on the common template mesh (inter-subject correspondence)
  |  wb_command -metric-resample  (BARYCENTRIC)
  |  sphere: common template mesh -> low-resolution icosphere
  v
thickness on the low-resolution grid (what gets stacked for the gradients)
```

## What it depends on

Only on file *roles* and a naming convention (native thickness, native
medial-wall mask, native midthickness, native-to-template registration
sphere, template sphere/area), matching common structural-pipeline
(BIDS-derivatives-style) outputs -- nothing here is specific to any one
study or cohort. Validated on dHCP neonatal data; reusable with any dataset
that shares that layout and has been registered to the same common template
space.

## Key functions

| Function | Role |
| --- | --- |
| `resample_hemisphere` | One subject, one hemisphere -> thickness + ROI on the low-res grid. |
| `resample_subject` | Both hemispheres for one subject. |
| `stack_lowres` | Concatenate L then R into flat `(thickness, roi)` arrays. |
| `session_dirs` | `(anat_dir, xfm_dir)` for a `sub-XXX/ses-YYY` directory. |
| `validate_sessions` | Drop sessions with a missing/unreadable required file (e.g. a truncated download). |
| `find_one` / `first_array` / `save_shape` / `run` | Small file/subprocess utilities. |
