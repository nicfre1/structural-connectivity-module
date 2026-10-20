# gradients.py

Single entry point for the `extra-functions/` scripts. Each step stays its own
script and still runs on its own; `gradients.py` only picks the right one from
the first word and hands it the remaining options, so you do not need to
remember thirteen file names.

```bash
python gradients.py --help                # list the commands
python gradients.py <command> --help      # options of one command
python gradients.py <command> [options]
```

## Commands

| Command | Runs | Needs |
| --- | --- | --- |
| `convert` | `drawem32_to_mcribs.py` | |
| `compute` | `compute_covariance_gradients.py` | |
| `plot` | `plot_group_gradients_brain.py` | VTK |
| `plot-atlas` | `plot_group_gradients_brain_atlas.py` | VTK |
| `compute-vertex` | `compute_covariance_gradients_vertex.py` | |
| `plot-vertex` | `plot_group_gradients_brain_vertex.py` | VTK |
| `build-grid` | `build_vertex_grid_assets.py` | Connectome Workbench |
| `resample-subject` | `resample_subject_thickness_vertex.py` | Connectome Workbench |
| `build-cohort` | `build_cohort_thickness_vertex.py` | Connectome Workbench |
| `example` | `generate_example_thickness.py` | |
| `example-vertex` | `generate_example_thickness_vertex.py` | |
| `example-native` | `generate_example_native_surfaces.py` | Connectome Workbench |

Everything after the command is passed to the script unchanged, including
`--help`, so the options of each step are documented in its own page (see
[README.md](README.md)).

## Choosing the Python interpreter

By default a command runs with the same Python that runs `gradients.py`. The
brain renders (`plot`, `plot-atlas`, `plot-vertex`) need a Python with a VTK
stack (see `../requirements-brain.txt`), often a different environment. Pass
it with `--python`; this option is consumed by `gradients.py` and is not
forwarded to the script:

```bash
python gradients.py plot --python /path/to/vtk-python \
    --gradient-csv group_gradients_rescaled_minus1_1.csv --output brain.png
```

## What it does not do

It adds no logic and no dependency. Every command runs as a separate process
(so a VTK crash cannot take anything else down, and `convert`/`compute` do
not import nibabel or VTK). Anything the scripts accept is accepted here, and
the scripts keep working when called directly.

When adding a script to `extra-functions/`, register it in `COMMANDS` in
`gradients.py`; a test fails if a script has no command.

## Example: DrawEM32 to MCRIBS, then covariance and gradients

```bash
python gradients.py example --out_dir example_data --n_subjects 30
python gradients.py convert \
    example_data/drawem32_thickness_lh_stats.tsv \
    example_data/drawem32_thickness_rh_stats.tsv \
    --mapping mapping_example.json --out_dir example_out
python gradients.py compute \
    example_out/thickness_68_lh_stats.tsv \
    example_out/thickness_68_rh_stats.tsv \
    --out_dir example_out/results
```
