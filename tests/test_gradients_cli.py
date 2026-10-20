import sys
from pathlib import Path

import pytest

EXTRA = Path(__file__).resolve().parents[1] / "extra-functions"
sys.path.insert(0, str(EXTRA))

import gradients as cli  # noqa: E402

NOT_COMMANDS = {"gradients.py", "vertex_resample_utils.py"}


def test_every_command_points_to_an_existing_script():
    for name, (script, _preset, _description, _needs) in cli.COMMANDS.items():
        assert (EXTRA / script).exists(), f"{name} -> missing {script}"


def test_every_script_in_extra_functions_is_reachable_from_a_command():
    registered = {script for script, _p, _d, _n in cli.COMMANDS.values()}
    scripts = {p.name for p in EXTRA.glob("*.py")} - NOT_COMMANDS
    assert scripts == registered, (
        "scripts without a command: "
        f"{sorted(scripts - registered)}; commands without a script: "
        f"{sorted(registered - scripts)}"
    )


def test_build_command_forwards_options_and_defaults_to_current_python():
    argv = cli.build_command("compute", ["a.tsv", "b.tsv", "--out_dir", "res"])
    assert argv[0] == sys.executable
    assert Path(argv[1]).name == "compute_covariance_gradients.py"
    assert argv[2:] == ["a.tsv", "b.tsv", "--out_dir", "res"]


@pytest.mark.parametrize("given", [["--python", "/vtk/py"], ["--python=/vtk/py"]])
def test_python_option_selects_interpreter_and_is_not_forwarded(given):
    argv = cli.build_command("plot", [*given, "--output", "x.png"])
    assert argv[0] == "/vtk/py"
    assert "--python" not in argv and not any(a.startswith("--python=") for a in argv)
    assert argv[2:] == ["dk", "--output", "x.png"]


def test_commands_sharing_a_script_select_its_mode_with_a_preset():
    modes = {
        "plot": "dk", "plot-atlas": "atlas", "plot-vertex": "vertex",
        "example": "drawem32", "example-vertex": "vertex", "example-native": "native",
    }
    for command, mode in modes.items():
        argv = cli.build_command(command, ["--out_dir", "x"])
        assert argv[2] == mode and argv[3:] == ["--out_dir", "x"]


def test_python_option_without_value_is_an_error():
    with pytest.raises(SystemExit):
        cli.build_command("plot", ["--python"])


def test_unknown_command_suggests_the_closest_one():
    with pytest.raises(SystemExit) as error:
        cli.build_command("comput", [])
    assert "compute" in str(error.value)


def test_help_lists_every_command(capsys):
    assert cli.main([]) == 0
    out = capsys.readouterr().out
    for name in cli.COMMANDS:
        assert name in out


def test_chain_runs_through_the_single_entry_point(tmp_path):
    example, converted, results = tmp_path / "ex", tmp_path / "mc", tmp_path / "res"
    assert cli.main(["example", "--out_dir", str(example), "--n_subjects", "12"]) == 0
    assert cli.main([
        "convert",
        str(example / "drawem32_thickness_lh_stats.tsv"),
        str(example / "drawem32_thickness_rh_stats.tsv"),
        "--mapping", str(EXTRA / "mapping_example.json"),
        "--out_dir", str(converted), "--strict",
    ]) == 0
    assert cli.main([
        "compute",
        str(converted / "thickness_68_lh_stats.tsv"),
        str(converted / "thickness_68_rh_stats.tsv"),
        "--out_dir", str(results), "--no_brain_plot",
    ]) == 0
    assert (results / "group_covariance" / "group_gradients.csv").exists()
