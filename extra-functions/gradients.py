#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single entry point for the extra-functions scripts.

Each step stays its own script (and still runs on its own); this only picks
the right one from the first word and hands it the remaining options, so you
do not need to remember thirteen file names.

See others/gradients.md for details.

Examples:
gradients.py --help
gradients.py convert lh.tsv rh.tsv --mapping mapping.json --out_dir out
gradients.py compute out/thickness_68_lh_stats.tsv out/thickness_68_rh_stats.tsv --out_dir res
gradients.py compute --help
gradients.py plot --python /path/to/vtk-python --gradient-csv g.csv --output brain.png
"""

import difflib
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# command -> (script, one-line description, extra requirement shown in --help)
COMMANDS = {
    "convert": (
        "drawem32_to_mcribs.py",
        "DrawEM32 thickness TSVs -> MCRIBS Desikan-Killiany (68 regions)",
        "",
    ),
    "compute": (
        "compute_covariance_gradients.py",
        "covariance matrices + group/individual gradients (regions)",
        "",
    ),
    "plot": (
        "plot_group_gradients_brain.py",
        "brain render of gradients, Desikan-Killiany (fsaverage5)",
        "VTK",
    ),
    "plot-atlas": (
        "plot_group_gradients_brain_atlas.py",
        "brain render for ANY atlas (your surface + annotation)",
        "VTK",
    ),
    "compute-vertex": (
        "compute_covariance_gradients_vertex.py",
        "covariance matrices + gradients, one row per cortical vertex",
        "",
    ),
    "plot-vertex": (
        "plot_group_gradients_brain_vertex.py",
        "brain render of vertex-wise gradients on your group surface",
        "VTK",
    ),
    "build-grid": (
        "build_vertex_grid_assets.py",
        "low-resolution icosphere grid from a dHCP symmetric template",
        "Workbench",
    ),
    "resample-subject": (
        "resample_subject_thickness_vertex.py",
        "native surface thickness of ONE subject -> common grid",
        "Workbench",
    ),
    "build-cohort": (
        "build_cohort_thickness_vertex.py",
        "whole cohort -> the cohort file compute-vertex reads",
        "Workbench",
    ),
    "example": (
        "generate_example_thickness.py",
        "synthetic DrawEM32 thickness files (no real data needed)",
        "",
    ),
    "example-vertex": (
        "generate_example_thickness_vertex.py",
        "synthetic vertex-wise cohort file (no real data needed)",
        "",
    ),
    "example-native": (
        "generate_example_native_surfaces.py",
        "synthetic template + subject surfaces for the resampling chain",
        "Workbench",
    ),
}

REQUIREMENT_NOTES = {
    "VTK": "needs a VTK-capable Python: pass --python /path/to/that/python",
    "Workbench": "needs Connectome Workbench (wb_command)",
}


def help_text():
    width = max(len(name) for name in COMMANDS)
    lines = [
        "usage: gradients.py <command> [options]",
        "",
        "Single entry point for the extra-functions scripts. Each command runs",
        "the matching script unchanged; options after the command go to it.",
        "",
        "commands:",
    ]
    for name, (_script, description, needs) in COMMANDS.items():
        suffix = f"  [{needs}]" if needs else ""
        lines.append(f"  {name:<{width}}  {description}{suffix}")
    lines += [
        "",
        "Run `gradients.py <command> --help` for the options of one command.",
    ]
    for key, note in REQUIREMENT_NOTES.items():
        lines.append(f"[{key}] {note}")
    return "\n".join(lines)


def split_python_option(args):
    """Remove `--python PATH` / `--python=PATH` from args -> (python, rest)."""
    rest, python, i = [], None, 0
    while i < len(args):
        arg = args[i]
        if arg == "--python":
            if i + 1 >= len(args):
                raise SystemExit("gradients.py: --python needs a path")
            python = args[i + 1]
            i += 2
        elif arg.startswith("--python="):
            python = arg.split("=", 1)[1]
            i += 1
        else:
            rest.append(arg)
            i += 1
    return python, rest


def build_command(name, args):
    """-> argv list that runs the script behind `name` with `args`."""
    if name not in COMMANDS:
        close = difflib.get_close_matches(name, COMMANDS, n=1)
        hint = f" Did you mean '{close[0]}'?" if close else ""
        raise SystemExit(
            f"gradients.py: unknown command '{name}'.{hint} "
            "Run `gradients.py --help` to list the commands."
        )
    python, rest = split_python_option(list(args))
    script = HERE / COMMANDS[name][0]
    return [python or sys.executable, str(script), *rest]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(help_text())
        return 0
    name, args = argv[0], argv[1:]
    return subprocess.call(build_command(name, args))


if __name__ == "__main__":
    sys.exit(main())
