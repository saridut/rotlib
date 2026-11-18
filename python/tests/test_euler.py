"""Public Euler API: conventions, all sequences, and singularity boundaries.

Angles are stored in X/Y/Z order even when the rotation sequence differs.
Expected matrices come from products of elementary rotations, not Rotlib.
"""

import math

import numpy as np
import pytest
from _rotation_cases import (
    LAYOUTS,
    MODES,
    VECTORS,
    check_axis_angle,
    check_canonical,
    close,
    euler_matrix,
    euler_quaternion,
    orientation_payload,
    quaternion_matrix,
    tensor_reference,
)

import rotlib

ANGLES = np.array([0.3, -0.4, 0.7])


@pytest.mark.parametrize("seq, world", MODES)
def test_euler_rotmat(seq, world):
    expected = euler_matrix(ANGLES, seq, world)
    close(rotlib.euler_rotmat(ANGLES, seq, world), expected)
    # A noncommuting example distinguishes fixed and moving axes.
    assert not np.allclose(expected, euler_matrix(ANGLES, seq, not world))


@pytest.mark.parametrize("seq, world", MODES)
@pytest.mark.parametrize(
    "angles", [np.zeros(3), ANGLES, np.array([2.8, -2.4, 2.1])]
)
def test_euler_factor_rotmat(seq, world, angles):
    expected = euler_matrix(angles, seq, world)
    factors = rotlib.euler_factor_rotmat(expected, seq, world)
    middle = "XYZ".index(seq[1])
    assert abs(factors[middle]) <= np.pi / 2
    assert np.all(np.abs(factors) <= np.pi)
    close(euler_matrix(factors, seq, world), expected)


@pytest.mark.parametrize("seq, world", MODES)
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("offset", [-1e-8, -5e-15, 0.0, 5e-15, 1e-8])
def test_euler_factor_rotmat_near_gimbal_lock(seq, world, sign, offset):
    angles = ANGLES.copy()
    middle = "XYZ".index(seq[1])
    angles[middle] = sign * (np.pi / 2 + offset)
    expected = euler_matrix(angles, seq, world)
    factors = rotlib.euler_factor_rotmat(expected, seq, world)
    assert abs(factors[middle]) <= np.pi / 2
    assert np.all(np.abs(factors) <= np.pi)
    close(euler_matrix(factors, seq, world), expected)
    if abs(offset) < 1e-14:
        assert factors["XYZ".index(seq[0])] == 0.0


@pytest.mark.parametrize("seq, world", MODES)
def test_euler_factor_rotmat_custom_tolerance(seq, world):
    angles = ANGLES.copy()
    angles["XYZ".index(seq[1])] = np.pi / 2 - 5e-5
    matrix = euler_matrix(angles, seq, world)
    ordinary = rotlib.euler_factor_rotmat(matrix, seq, world)
    loose = rotlib.euler_factor_rotmat(matrix, seq, world, tol=1e-4)
    close(euler_matrix(ordinary, seq, world), matrix)
    assert loose["XYZ".index(seq[0])] == 0.0
    # A loose cutoff deliberately discards the small outer-angle contribution.
    close(euler_matrix(loose, seq, world), matrix, atol=1e-4)


@pytest.mark.parametrize("seq, world", MODES)
def test_euler_shiftmat_and_rotate_vectors(seq, world):
    expected = euler_matrix(ANGLES, seq, world)
    close(rotlib.euler_shiftmat(ANGLES, seq, world), expected)
    close(rotlib.euler_shiftmat(ANGLES, seq, world, forward=True), expected.T)
    close(
        rotlib.euler_rotate_vectors(VECTORS, ANGLES, seq, world),
        VECTORS @ expected.T,
    )


@pytest.mark.parametrize("v", LAYOUTS)
def test_euler_rotate_vectors_layout(v):
    result = rotlib.euler_rotate_vectors(v, ANGLES)
    assert result.shape == v.shape
    close(result, v @ euler_matrix(ANGLES, "XYZ", True).T)


@pytest.mark.parametrize("seq, world", MODES)
@pytest.mark.parametrize("forward", [False, True])
def test_euler_shift_vectors(seq, world, forward):
    matrix = euler_matrix(ANGLES, seq, world)
    s = matrix.T if forward else matrix
    close(
        rotlib.euler_shift_vectors(VECTORS, ANGLES, seq, world, forward),
        VECTORS @ s.T,
    )


@pytest.mark.parametrize("v", LAYOUTS)
def test_euler_shift_vectors_layout_and_defaults(v):
    result = rotlib.euler_shift_vectors(v, ANGLES)
    assert result.shape == v.shape
    close(result, v @ euler_matrix(ANGLES, "XYZ", True).T)


@pytest.mark.parametrize(
    "rank, function",
    [
        (2, rotlib.euler_shift_tensor2),
        (3, rotlib.euler_shift_tensor3),
    ],
)
@pytest.mark.parametrize("seq, world", MODES)
@pytest.mark.parametrize("forward", [False, True])
def test_euler_shift_tensors(rank, function, seq, world, forward):
    a = np.arange(3**rank, dtype=np.float64).reshape((3,) * rank) - 4.0
    matrix = euler_matrix(ANGLES, seq, world)
    s = matrix.T if forward else matrix
    close(function(a, ANGLES, seq, world, forward), tensor_reference(a, s))


@pytest.mark.parametrize("seq, world", MODES)
@pytest.mark.parametrize(
    "angles", [ANGLES, np.zeros(3), np.array([np.pi, 0.0, 0.0])]
)
def test_euler_to_quat_aa_dcm(seq, world, angles):
    expected = euler_matrix(angles, seq, world)
    q = rotlib.euler_to_quat(angles, seq, world)
    check_canonical(q)
    close(math.hypot(*q), 1.0)
    close(quaternion_matrix(q), expected)
    check_axis_angle(*rotlib.euler_to_aa(angles, seq, world), expected)
    close(rotlib.euler_to_dcm(angles, seq, world), expected.T)


@pytest.mark.parametrize("seq, world", MODES)
@pytest.mark.parametrize("to_seq, to_world", MODES)
def test_euler_to_euler(seq, world, to_seq, to_world):
    out = rotlib.euler_to_euler(ANGLES, seq, world, to_seq, to_world)
    close(
        euler_matrix(out, to_seq, to_world), euler_matrix(ANGLES, seq, world)
    )


@pytest.mark.parametrize(
    "representation", ["quat", "euler", "axis_angle", "dcm"]
)
@pytest.mark.parametrize("seq, world", MODES)
def test_euler_from_any(representation, seq, world):
    payload, expected = orientation_payload(representation)
    out = rotlib.euler_from_any(payload, seq, world)
    close(euler_matrix(out, seq, world), expected)


@pytest.mark.parametrize(
    "representation", ["quat", "euler", "axis_angle", "dcm"]
)
def test_euler_from_any_forwards_tolerance(representation):
    angles = np.array([0.7, np.pi / 2 - 5e-5, -0.4])
    q = euler_quaternion(angles, "XYZ", True)
    axis = q[1:] / math.hypot(*q[1:])
    angle = 2 * math.atan2(math.hypot(*q[1:]), q[0])
    payloads = {
        "quat": {"quat": q},
        "euler": {"euler": angles, "seq": "XYZ", "world": True},
        "axis_angle": {"axis": axis, "angle": angle},
        "dcm": {"dcm": euler_matrix(angles, "XYZ", True).T},
    }
    payload = {"repr": representation, **payloads[representation]}
    ordinary = rotlib.euler_from_any(payload, "XYZ", True)
    loose = rotlib.euler_from_any(payload, "XYZ", True, tol=1e-4)
    close(
        euler_matrix(ordinary, "XYZ", True), euler_matrix(angles, "XYZ", True)
    )
    assert loose[0] == 0.0


def test_euler_from_any_unknown_representation():
    with pytest.raises(ValueError, match="Unrecognized orientation repr"):
        rotlib.euler_from_any({"repr": "unknown"}, "XYZ", True)


def test_euler_to_euler_forwards_tolerance():
    angles = np.array([0.7, np.pi / 2 - 5e-5, -0.4])
    ordinary = rotlib.euler_to_euler(angles, "XYZ", True, "XYZ", True)
    loose = rotlib.euler_to_euler(angles, "XYZ", True, "XYZ", True, tol=1e-4)
    close(
        euler_matrix(ordinary, "XYZ", True), euler_matrix(angles, "XYZ", True)
    )
    assert loose[0] == 0.0
