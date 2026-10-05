#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Vertex-wise cortical thickness for a whole COHORT: loop
resample_subject_thickness_vertex.py over a batch of subjects, build a
group-average midthickness surface, smooth, define the group cortical mask,
and stack the result into the cohort .npz that
compute_covariance_gradients_vertex.py consumes.

See others/build_cohort_thickness_vertex.md for details.

On top of the single-subject step, this also:

  1. Builds the group-average midthickness surface in the common template
     space (each subject's native midthickness resampled onto the 32k mesh,
     then averaged, then downsampled to the ~N k grid). This is the
     geometric support for smoothing and for cortical-surface renders
     (plot_group_gradients_brain_vertex.py).
  2. Smooths each subject's thickness on that group surface, Gaussian
     kernel, --fwhm mm (default 10, standard in morphometry; needed at the
     vertex scale to control single-vertex noise).
  3. Defines the GROUP cortical mask = vertices valid in at least
     --min_coverage (default 1.0 = all) of the subjects.

Output: <output_dir>/cohort_thickness_vertex_dens-{N}k.npz
  thickness   : (n_subjects, V)  smoothed thickness, restricted to the group mask
  mask        : (2N,) bool   group cortical mask
  hemi        : (2N,) 'L'/'R'
  subject_ids : (n_subjects,)
  fwhm, lowres, min_coverage
  + assets/groupmid_hemi-{L,R}_dhcpSym_dens-{N}k.surf.gii

Subject selection: either --manifest (a CSV with columns
participant_id[,session_id], your own cohort list -- nothing here ships any
subject list) or the first --n_subjects sessions found under
--derivatives_root, plus any --must_include.

Example:
build_cohort_thickness_vertex.py \\
    --derivatives_root /path/to/derivatives --output_dir ./cohort \\
    --assets_dir ./assets --template_dir /path/to/dhcp-symmetric-template \\
    --wb_command /Applications/ConnectomeWorkbench.app/Contents/usr/bin/wb_command \\
    --lowres 10k --fwhm 10 --manifest my_cohort.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np

from vertex_resample_utils import (
    HEMISPHERES,
    find_one,
    first_array,
    run,
    session_dirs,
    validate_sessions,
)
from resample_subject_thickness_vertex import build_subject_vertex_npz


def load_surface(path: Path) -> tuple[np.ndarray, np.ndarray]:
    image = nib.load(str(path))
    return (
        np.asarray(image.darrays[0].data, dtype=np.float64),
        np.asarray(image.darrays[1].data, dtype=np.int32),
    )


def save_surface(points: np.ndarray, faces: np.ndarray, path: Path) -> None:
    pts = nib.gifti.GiftiDataArray(points.astype(np.float32), intent="NIFTI_INTENT_POINTSET")
    tri = nib.gifti.GiftiDataArray(faces.astype(np.int32), intent="NIFTI_INTENT_TRIANGLE")
    nib.save(nib.gifti.GiftiImage(darrays=[pts, tri]), str(path))


def build_group_midthickness(
    *, sessions: list[Path], template_dir: Path, assets_dir: Path,
    work_dir: Path, wb_command: str, lowres: str,
) -> None:
    """Group-average midthickness, 32k then ~N k grid, per hemisphere."""
    for hemisphere, short in HEMISPHERES:
        structure = "CORTEX_LEFT" if hemisphere == "left" else "CORTEX_RIGHT"
        dhcp32k_sphere = (
            template_dir
            / f"week-40_hemi-{hemisphere}_space-dhcpSym_dens-32k_sphere.surf.gii"
        )
        accum: np.ndarray | None = None
        faces: np.ndarray | None = None
        used = 0
        for session_dir in sessions:
            anat_dir, xfm_dir = session_dirs(session_dir)
            try:
                native_mid = find_one(anat_dir, f"*_hemi-{hemisphere}_midthickness.surf.gii")
                reg_sphere = find_one(
                    xfm_dir,
                    f"*_hemi-{hemisphere}_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii",
                )
            except (FileNotFoundError, RuntimeError):
                continue
            out_32k = work_dir / f"mid32k_{short}_{session_dir.parent.name}.surf.gii"
            try:
                run(
                    wb_command, "-surface-resample",
                    native_mid, reg_sphere, dhcp32k_sphere, "BARYCENTRIC", out_32k,
                )
                pts, tri = load_surface(out_32k)
            except Exception as error:  # noqa: BLE001
                print(f"    [SKIP midthickness] {session_dir.parent.name} hemi-{short}: {error}")
                out_32k.unlink(missing_ok=True)
                continue
            out_32k.unlink(missing_ok=True)
            accum = pts if accum is None else accum + pts
            faces = tri if faces is None else faces
            used += 1
        if accum is None:
            raise RuntimeError(f"hemi-{hemisphere}: no usable midthickness")
        group_32k = accum / used
        group_32k_path = work_dir / f"groupmid_{short}_32k.surf.gii"
        save_surface(group_32k, faces, group_32k_path)

        target_sphere = assets_dir / f"sphere_hemi-{short}_dhcpSym_dens-{lowres}.surf.gii"
        group_low_path = assets_dir / f"groupmid_hemi-{short}_dhcpSym_dens-{lowres}.surf.gii"
        run(
            wb_command, "-surface-resample",
            group_32k_path, dhcp32k_sphere, target_sphere, "BARYCENTRIC", group_low_path,
        )
        run(wb_command, "-set-structure", group_low_path, structure)
        print(f"  group midthickness hemi-{short}: average of {used} subjects -> {group_low_path.name}")


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--derivatives_root", required=True, type=Path)
    p.add_argument("--output_dir", required=True, type=Path)
    p.add_argument("--assets_dir", type=Path, default=Path("assets"))
    p.add_argument("--template_dir", required=True, type=Path)
    p.add_argument("--wb_command", required=True, type=str)
    p.add_argument("--lowres", default="10k")
    p.add_argument("--fwhm", type=float, default=10.0)
    p.add_argument("--manifest", type=Path, default=None,
                    help="CSV with columns participant_id[,session_id]: uses "
                         "exactly these sessions (takes priority over "
                         "--n_subjects/--must_include). Your own cohort list; "
                         "nothing here ships one.")
    p.add_argument("--n_subjects", type=int, default=25)
    p.add_argument("--must_include", nargs="*", default=[])
    p.add_argument("--min_coverage", type=float, default=1.0,
                    help="Fraction of subjects a vertex must be valid in to "
                         "enter the group mask (1.0 = all).")
    p.add_argument("--thickness_desc", default="thickness")
    return p


def main():
    args = _build_arg_parser().parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.output_dir / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    # --- subject selection ---------------------------------------------
    sessions: list[Path] = []
    if args.manifest is not None:
        import pandas as pd

        manifest = pd.read_csv(args.manifest)
        has_session = "session_id" in manifest.columns
        for _, row in manifest.iterrows():
            sub = str(row["participant_id"])
            if not sub.startswith("sub-"):
                sub = f"sub-{sub}"
            subject_dir = args.derivatives_root / sub
            if has_session and pd.notna(row["session_id"]):
                cand = subject_dir / f"ses-{int(row['session_id'])}"
                ses_dir = cand if cand.is_dir() else None
            else:
                ses_dir = None
            if ses_dir is None:
                found = sorted(subject_dir.glob("ses-*"))
                ses_dir = found[0] if found else None
            if ses_dir is not None:
                sessions.append(ses_dir)
            else:
                print(f"  [SKIP manifest] {sub}: session not found on disk")
        print(f"Manifest {args.manifest.name}: {len(sessions)}/{len(manifest)} sessions found")
    else:
        subject_dirs = sorted(d for d in args.derivatives_root.glob("sub-*") if d.is_dir())
        selected = subject_dirs[: args.n_subjects]
        selected_ids = {d.name for d in selected}
        for forced in args.must_include:
            if forced not in selected_ids and (args.derivatives_root / forced).exists():
                selected.append(args.derivatives_root / forced)
                selected_ids.add(forced)
        for subject_dir in selected:
            ses = sorted(subject_dir.glob("ses-*"))
            if ses:
                sessions.append(ses[0])

    # --- validation: drop incomplete downloads ---------------------------
    sessions, unreadable = validate_sessions(sessions, args.thickness_desc)
    if unreadable:
        print(f"\n[{len(unreadable)} session(s) dropped -- missing or unreadable file]")
        for subject, reason in unreadable:
            print(f"  [EXCLUDED] {subject}: {reason}")
    print(f"Retained cohort: {len(sessions)} sessions\n")

    # --- 1. group-average midthickness surface ----------------------------
    print("=== group-average midthickness (dHCPsym40) ===")
    build_group_midthickness(
        sessions=sessions, template_dir=args.template_dir, assets_dir=args.assets_dir,
        work_dir=work_dir, wb_command=args.wb_command, lowres=args.lowres,
    )
    group_mid = {
        short: args.assets_dir / f"groupmid_hemi-{short}_dhcpSym_dens-{args.lowres}.surf.gii"
        for _, short in HEMISPHERES
    }

    # --- 2. vertex thickness + smoothing, subject by subject ---------------
    print("\n=== vertex thickness + FWHM %.1f mm smoothing ===" % args.fwhm)
    rows_thickness: list[np.ndarray] = []
    rows_roi: list[np.ndarray] = []
    subject_ids: list[str] = []
    skipped: list[tuple[str, str]] = []

    for session_dir in sessions:
        subject = f"{session_dir.parent.name}_{session_dir.name}"
        npz_path = args.output_dir / "vertex_thickness" / f"{subject}_thickness_vertex.npz"
        subj_work = work_dir / subject
        try:
            build_subject_vertex_npz(
                session_dir=session_dir,
                assets_dir=args.assets_dir,
                template_dir=args.template_dir,
                wb_command=args.wb_command,
                lowres=args.lowres,
                output_npz=npz_path,
                work_dir=subj_work,
                thickness_desc=args.thickness_desc,
            )
        except Exception as error:  # noqa: BLE001
            skipped.append((subject, str(error)))
            continue

        # surface smoothing per hemisphere, on the group midthickness
        try:
            smoothed_parts: list[np.ndarray] = []
            roi_parts: list[np.ndarray] = []
            for _, short in HEMISPHERES:
                raw = subj_work / f"thick_dhcpSym{args.lowres}_{short}.shape.gii"
                roi = subj_work / f"roi_dhcpSym{args.lowres}_{short}.shape.gii"
                sm = subj_work / f"thick_dhcpSym{args.lowres}_{short}_sm.shape.gii"
                run(
                    args.wb_command, "-metric-smoothing",
                    group_mid[short], raw, args.fwhm, sm, "-fwhm", "-roi", roi,
                )
                smoothed_parts.append(first_array(sm).astype(np.float32))
                roi_parts.append(first_array(roi) > 0.5)
        except Exception as error:  # noqa: BLE001
            skipped.append((subject, f"smoothing: {error}"))
            continue

        rows_thickness.append(np.concatenate(smoothed_parts))
        rows_roi.append(np.concatenate(roi_parts))
        subject_ids.append(subject)

    if not rows_thickness:
        raise SystemExit("No usable subject.")

    thickness_full = np.vstack(rows_thickness)          # (n_subjects, 2N)
    roi_stack = np.vstack(rows_roi)                     # (n_subjects, 2N)
    n_per_hemi = thickness_full.shape[1] // 2

    coverage = roi_stack.mean(axis=0)
    mask = coverage >= args.min_coverage
    # safety: never an infinite or non-positive thickness vertex
    mask &= np.isfinite(thickness_full).all(axis=0)
    mask &= (thickness_full > 0).all(axis=0)

    hemi = np.array(["L"] * n_per_hemi + ["R"] * n_per_hemi)
    thickness_masked = thickness_full[:, mask].astype(np.float32)

    out_npz = args.output_dir / f"cohort_thickness_vertex_dens-{args.lowres}.npz"
    np.savez_compressed(
        out_npz,
        thickness=thickness_masked,
        mask=mask,
        hemi=hemi,
        subject_ids=np.array(subject_ids),
        fwhm=np.float64(args.fwhm),
        lowres=args.lowres,
        min_coverage=np.float64(args.min_coverage),
        n_per_hemi=np.int64(n_per_hemi),
    )

    print(f"\n=== done: {len(subject_ids)} subjects, {mask.sum()} group cortical vertices "
          f"(of {mask.size}) ===")
    print(f"cohort matrix: {thickness_masked.shape} -> {out_npz}")
    for name, reason in skipped:
        print(f"  [SKIP] {name}: {reason}")


if __name__ == "__main__":
    main()
