"""Axis-angle public API: principal angles, limiting cases and transformations."""

import math

import numpy as np
import pytest
from _rotation_cases import (
    CASES,
    LAYOUTS,
    MIXED_AXIS,
    MODES,
    axis_matrix,
    check_axis_angle,
    check_canonical,
    close,
    euler_matrix,
    orientation_payload,
    quaternion_matrix,
    tensor_reference,
)
from numpy.testing import assert_allclose

import rotlib


@pytest.mark.parametrize("angle", [-7.0, -0.3, 0.0, 0.7, np.pi, 4.0, 1e20])
@pytest.mark.parametrize("normalize", [False, True])
def test_aa_fix(angle, normalize):
    axis = np.array([1.0, 2.0, 3.0])
    original = axis.copy()
    ax, ang = rotlib.aa_fix(axis, angle, normalize=normalize)
    close(axis, original, atol=0.0)
    close(math.hypot(*ax), 1.0 if normalize else math.hypot(*original))
    check_axis_angle(
        ax / math.hypot(*ax),
        ang,
        axis_matrix(original / math.hypot(*original), angle),
    )


@pytest.mark.parametrize("scale", [1e-300, 1e308, np.nextafter(0.0, 1.0)])
def test_aa_fix_scale_safe_normalization(scale):
    axis, angle = rotlib.aa_fix(np.full(3, scale), 0.7)
    close(axis, np.full(3, 1 / np.sqrt(3.0)))
    close(angle, 0.7)


def test_aa_fix_zero_axis_rejected():
    with pytest.raises(ValueError, match="zero"):
        rotlib.aa_fix(np.zeros(3), 0.7)


def test_aa_rand():
    first, second = np.random.default_rng(83), np.random.default_rng(83)
    for _ in range(64):
        axis, angle = rotlib.aa_rand(first)
        other_axis, other_angle = rotlib.aa_rand(second)
        close(axis, other_axis, atol=0.0)
        assert angle == other_angle
        close(math.hypot(*axis), 1.0)
        assert 0.0 <= angle <= np.pi


@pytest.mark.parametrize("axis, angle", CASES)
def test_aa_rotmat(axis, angle):
    r = rotlib.aa_rotmat(axis, angle)
    close(r, axis_matrix(axis, angle))
    close(r @ r.T, np.eye(3))
    close(np.linalg.det(r), 1.0)


def test_aa_rotmat_small_second_order_term():
    # With z-axis component zero, R01 is purely second order in the angle.
    axis = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    angle = 1e-8
    expected = 2 * math.sin(angle / 2) ** 2 * axis[0] * axis[1]
    assert_allclose(
        rotlib.aa_rotmat(axis, angle)[0, 1], expected, rtol=2e-14, atol=0.0
    )


@pytest.mark.parametrize(
    "axis",
    [
        np.array([3.0, -2.0, 1.0]) / np.sqrt(14.0),
        np.array([1.0, -3.0, 2.0]) / np.sqrt(14.0),
        np.array([1.0, 2.0, -3.0]) / np.sqrt(14.0),
        np.ones(3) / np.sqrt(3.0),
    ],
)
@pytest.mark.parametrize(
    "angle",
    [
        0.0,
        1e-200,
        1e-8,
        np.pi / 2 - 1e-8,
        np.pi / 2 + 1e-8,
        np.pi - 1e-8,
        np.pi,
    ],
)
def test_aa_from_rotmat(axis, angle):
    r = axis_matrix(axis, angle)
    ax, ang = rotlib.aa_from_rotmat(r)
    check_axis_angle(ax, ang, r)
    if angle == 0.0:
        close(ax, [1.0, 0.0, 0.0], atol=0.0)
    elif angle < 1e-7:
        assert_allclose(ang, angle, rtol=2e-14, atol=0.0)


@pytest.mark.parametrize("v", LAYOUTS)
def test_aa_rotate_vectors(v):
    out = rotlib.aa_rotate_vectors(v, MIXED_AXIS, 0.7)
    assert out.shape == v.shape
    close(out, v @ axis_matrix(MIXED_AXIS, 0.7).T)


@pytest.mark.parametrize("angle", [-0.7, 0.0, 0.7, np.pi, 4.0])
def test_aa_shiftmat_and_to_dcm(angle):
    r = axis_matrix(MIXED_AXIS, angle)
    close(rotlib.aa_shiftmat(MIXED_AXIS, angle, forward=True), r.T)
    close(rotlib.aa_shiftmat(MIXED_AXIS, angle), r)
    close(rotlib.aa_to_dcm(MIXED_AXIS, angle), r.T)


@pytest.mark.parametrize("v", LAYOUTS)
@pytest.mark.parametrize("forward", [False, True])
def test_aa_shift_vectors(v, forward):
    r = axis_matrix(MIXED_AXIS, 0.7)
    s = r.T if forward else r
    out = rotlib.aa_shift_vectors(v, MIXED_AXIS, 0.7, forward=forward)
    assert out.shape == v.shape
    close(out, v @ s.T)


@pytest.mark.parametrize(
    "rank, function",
    [
        (2, rotlib.aa_shift_tensor2),
        (3, rotlib.aa_shift_tensor3),
    ],
)
@pytest.mark.parametrize("forward", [False, True])
def test_aa_shift_tensors(rank, function, forward):
    a = np.arange(3**rank, dtype=np.float64).reshape((3,) * rank) - 4.0
    r = axis_matrix(MIXED_AXIS, 0.7)
    close(
        function(a, MIXED_AXIS, 0.7, forward=forward),
        tensor_reference(a, r.T if forward else r),
    )


@pytest.mark.parametrize("angle", [-7.0, 0.0, 0.7, np.pi, 4.0])
def test_aa_to_quat(angle):
    q = rotlib.aa_to_quat(MIXED_AXIS, angle)
    close(math.hypot(*q), 1.0)
    check_canonical(q)
    close(quaternion_matrix(q), axis_matrix(MIXED_AXIS, angle))


@pytest.mark.parametrize("seq, world", MODES)
def test_aa_to_euler(seq, world):
    angles = rotlib.aa_to_euler(MIXED_AXIS, 0.7, seq, world)
    close(euler_matrix(angles, seq, world), axis_matrix(MIXED_AXIS, 0.7))


@pytest.mark.parametrize(
    "representation", ["quat", "euler", "axis_angle", "dcm"]
)
def test_aa_from_any(representation):
    payload, matrix = orientation_payload(representation)
    check_axis_angle(*rotlib.aa_from_any(payload), matrix)


def test_aa_from_any_canonicalizes_existing_axis_angle():
    axis, angle = rotlib.aa_from_any(
        {
            "repr": "axis_angle",
            "axis": 3 * MIXED_AXIS,
            "angle": 4.0,
        }
    )
    check_axis_angle(axis, angle, axis_matrix(MIXED_AXIS, 4.0))


def test_aa_from_any_unknown_representation():
    with pytest.raises(ValueError, match="Unrecognized orientation repr"):
        rotlib.aa_from_any({"repr": "unknown"})
