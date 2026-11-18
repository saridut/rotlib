"""DCM public API: components map A to B; active matrices act in A."""

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

import rotlib


@pytest.mark.parametrize("axis, angle", CASES)
def test_dcm_from_frames(axis, angle):
    # Neither frame coincides with the underlying coordinate frame.
    a = np.array([[1.0, 2.0, 2.0], [2.0, 1.0, -2.0], [-2.0, 2.0, -1.0]]) / 3.0
    r_f = axis_matrix(axis, angle)
    b = a @ r_f.T
    d = rotlib.dcm_from_frames(a, b)
    # Check each direction cosine independently of matrix-product ordering.
    expected = np.array([[np.dot(bi, aj) for aj in a] for bi in b])
    close(d, expected)
    x_f = np.array([2.0, -3.0, 4.0])
    close(d @ (a @ x_f), b @ x_f)
    # R acts in A, rather than in the underlying F coordinates.
    close(rotlib.dcm_rotmat(d), a @ r_f @ a.T)


@pytest.mark.parametrize("axis, angle", CASES)
def test_dcm_rotmat_and_shiftmat(axis, angle):
    r = axis_matrix(axis, angle)
    d = r.T
    out = rotlib.dcm_rotmat(d)
    close(out, r)
    assert not np.shares_memory(out, d)
    for forward in (False, True):
        s = rotlib.dcm_shiftmat(d, forward=forward)
        close(s, d if forward else r)
        assert not np.shares_memory(s, d)


@pytest.mark.parametrize("v", LAYOUTS)
def test_dcm_rotate_vectors(v):
    r = axis_matrix(MIXED_AXIS, 0.7)
    out = rotlib.dcm_rotate_vectors(v, r.T)
    assert out.shape == v.shape
    close(out, v @ r.T)


@pytest.mark.parametrize("v", LAYOUTS)
@pytest.mark.parametrize("forward", [False, True])
def test_dcm_shift_vectors(v, forward):
    r = axis_matrix(MIXED_AXIS, 0.7)
    out = rotlib.dcm_shift_vectors(v, r.T, forward=forward)
    assert out.shape == v.shape
    close(out, v @ (r if forward else r.T))


@pytest.mark.parametrize(
    "rank, function",
    [
        (2, rotlib.dcm_shift_tensor2),
        (3, rotlib.dcm_shift_tensor3),
    ],
)
@pytest.mark.parametrize("forward", [False, True])
def test_dcm_shift_tensors(rank, function, forward):
    a = np.arange(3**rank, dtype=np.float64).reshape((3,) * rank) - 4.0
    r = axis_matrix(MIXED_AXIS, 0.7)
    close(
        function(a, r.T, forward=forward),
        tensor_reference(a, r.T if forward else r),
    )


@pytest.mark.parametrize("axis, angle", CASES)
def test_dcm_to_quat_and_aa(axis, angle):
    r = axis_matrix(axis, angle)
    q = rotlib.dcm_to_quat(r.T)
    close(math.hypot(*q), 1.0)
    check_canonical(q)
    close(quaternion_matrix(q), r)
    check_axis_angle(*rotlib.dcm_to_aa(r.T), r)


@pytest.mark.parametrize("seq, world", MODES)
def test_dcm_to_euler(seq, world):
    r = axis_matrix(MIXED_AXIS, 0.7)
    close(euler_matrix(rotlib.dcm_to_euler(r.T, seq, world), seq, world), r)


@pytest.mark.parametrize(
    "representation", ["quat", "euler", "axis_angle", "dcm"]
)
def test_dcm_from_any(representation):
    payload, matrix = orientation_payload(representation)
    out = rotlib.dcm_from_any(payload)
    assert out.shape == (3, 3)
    close(out, matrix.T)
    if representation == "dcm":
        assert not np.shares_memory(out, payload["dcm"])


def test_dcm_from_any_unknown_representation():
    with pytest.raises(ValueError, match="Unrecognized orientation repr"):
        rotlib.dcm_from_any({"repr": "unknown"})
