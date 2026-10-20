"""
Annotation/surface readers of plot_group_gradients_brain_atlas.py.

Own module because it needs nibabel (not in ../requirements.txt, only in
extra-functions/requirements-brain.txt): a module-level importorskip skips
this whole file where nibabel is absent (e.g. the standard CI runner)
without taking the other tests down with it.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra-functions"))

nib = pytest.importorskip("nibabel")

import plot_group_gradients_brain_atlas as atlas  # noqa: E402


def test_load_annotation_reads_freesurfer_annot(tmp_path):
    labels = np.array([0, 1, 1, 2, 0])
    names = ["bankssts", "insula", "unknown"]
    ctab = np.array([[10, 20, 30, 0], [40, 50, 60, 0], [70, 80, 90, 0]])
    path = tmp_path / "lh.test.annot"
    nib.freesurfer.write_annot(str(path), labels, ctab, names, fill_ctab=True)

    read_labels, read_names = atlas.load_annotation(path)
    assert read_names == names
    assert np.array_equal(read_labels, labels)


def test_load_surface_reads_gifti_points_and_faces(tmp_path):
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float32)
    faces = np.array([[0, 1, 2], [0, 1, 3]], dtype=np.int32)
    image = nib.gifti.GiftiImage(darrays=[
        nib.gifti.GiftiDataArray(points, intent="NIFTI_INTENT_POINTSET"),
        nib.gifti.GiftiDataArray(faces, intent="NIFTI_INTENT_TRIANGLE"),
    ])
    path = tmp_path / "lh.surf.gii"
    nib.save(image, str(path))

    read_points, read_faces = atlas.load_surface(path)
    assert np.allclose(read_points, points)
    assert np.array_equal(read_faces, faces)
