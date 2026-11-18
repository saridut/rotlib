"""Tests for every public quaternion functions; scalar-first Hamilton convention."""

import math

import numpy as np
import pytest
from _rotation_cases import (
    CASES,
    LAYOUTS,
    MIXED_AXIS,
    MODES,
    axis_quaternion,
    check_axis_angle,
    check_canonical,
    close,
    euler_matrix,
    euler_quaternion,
    hamilton,
    orientation_payload,
    quaternion_matrix,
    tensor_reference,
)
from numpy.testing import assert_allclose

import rotlib


@pytest.mark.parametrize(
    "values",
    [
        [-1.0, 0.0, 0.0, 0.0],
        [0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0, 0.0],
    ],
)
def test_quat_canonicalized(values):
    q = np.array(values)
    before = q.copy()
    assert rotlib.quat_canonicalized(q) is None
    check_canonical(q)
    close(quaternion_matrix(q), quaternion_matrix(before))


def test_quat_rand():
    # Reproducibility and output invariants, not a statistical uniformity test.
    first, second = np.random.default_rng(82), np.random.default_rng(82)
    for _ in range(64):
        q = rotlib.quat_rand(first)
        close(q, rotlib.quat_rand(second), atol=0.0)
        close(math.hypot(*q), 1.0)
        check_canonical(q)


def test_quat_identity():
    close(rotlib.quat_identity(), [1.0, 0.0, 0.0, 0.0], atol=0.0)


def test_quat_conjugated():
    q = np.array([1.0, 2.0, -3.0, 4.0])
    assert rotlib.quat_conjugated(q) is q
    close(q, [1.0, -2.0, 3.0, -4.0], atol=0.0)


@pytest.mark.parametrize("scale", [1.0, 1e200, 1e-200, 1e308, 1e-300])
def test_quat_inverted(scale):
    original = np.full(4, scale)
    q = original.copy()
    assert rotlib.quat_inverted(q) is q
    expected = np.array([0.25, -0.25, -0.25, -0.25]) / scale
    assert_allclose(q, expected, rtol=2e-14, atol=0.0)
    close(hamilton(original, q), [1.0, 0.0, 0.0, 0.0])


@pytest.mark.parametrize(
    "function", [rotlib.quat_inverted, rotlib.quat_normalized]
)
def test_zero_quaternion_rejected(function):
    with pytest.raises(ValueError, match="zero"):
        function(np.zeros(4))


@pytest.mark.parametrize(
    "scale",
    [
        3.0,
        1e-300,
        np.finfo(np.float64).max,
        np.nextafter(0.0, 1.0),
    ],
)
@pytest.mark.parametrize("canonicalize", [False, True])
def test_quat_normalized(scale, canonicalize):
    q = np.full(4, -scale)
    assert rotlib.quat_normalized(q, canonicalize=canonicalize) is q
    close(q, np.full(4, 0.5 if canonicalize else -0.5))


@pytest.mark.parametrize(
    "norm, expected",
    [
        (1.0, True),
        (0.0, False),
        (0.5, False),
        (2.0, False),
        (1e308, False),
        (1e-300, False),
    ],
)
def test_quat_is_normalized(norm, expected):
    assert (
        rotlib.quat_is_normalized(np.array([norm, 0.0, 0.0, 0.0])) == expected
    )


@pytest.mark.parametrize("delta, expected", [(0.5e-12, True), (2e-12, False)])
def test_quat_is_normalized_tolerance(delta, expected):
    q = np.array([1.0 + delta, 0.0, 0.0, 0.0])
    assert rotlib.quat_is_normalized(q, tol=1e-12) == expected
    q[0] = 1.0 - delta
    assert rotlib.quat_is_normalized(q, tol=1e-12) == expected


def test_quat_prod_hamilton_order_and_nonunit_values():
    i, j, k = np.eye(4)[1:]
    close(rotlib.quat_prod(i, j), k)
    close(rotlib.quat_prod(j, i), -k)
    p = np.array([1.0, 2.0, 3.0, 4.0])
    q = np.array([5.0, 6.0, 7.0, 8.0])
    close(rotlib.quat_prod(p, q), [-60.0, 12.0, 30.0, 24.0])
    close(rotlib.quat_prod(np.zeros(4), q), np.zeros(4))


@pytest.mark.parametrize("angle", [0.0, 1e-200, 0.7, np.pi, 1.5 * np.pi])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_quat_angle_between(angle, sign):
    q = sign * axis_quaternion(np.array([0.0, 0.0, 1.0]), angle)
    actual = rotlib.quat_angle_between(np.array([1.0, 0.0, 0.0, 0.0]), q)
    expected = min(angle, 2 * np.pi - angle)
    assert_allclose(actual, expected, rtol=2e-14, atol=0.0)
    assert rotlib.quat_angle_between(q, -q) == 0.0


@pytest.mark.parametrize("angle", [0.0, 1e-200, 1e-8, 0.7, np.pi])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
@pytest.mark.parametrize("t", [0.0, 0.2, 0.5, 1.0])
def test_quat_interpolate(angle, sign, t):
    axis = np.array([0.0, 0.0, 1.0])
    q1 = np.array([1.0, 0.0, 0.0, 0.0])
    # Exact scalar zero makes the half-turn path tie deterministic.
    q2 = axis_quaternion(axis, angle)
    if angle == np.pi:
        q2[0] = 0.0
    expected_angle = sign * t * angle if angle == np.pi else t * angle
    out = rotlib.quat_interpolate(q1, sign * q2, t)
    assert_allclose(
        out, axis_quaternion(axis, expected_angle), rtol=2e-14, atol=1e-16
    )
    close(math.hypot(*out), 1.0)
    if angle == 1e-200 and t > 0.0:
        assert_allclose(out[3], math.sin(t * angle / 2), rtol=2e-14, atol=0.0)


def test_quat_interpolate_nonidentity_start():
    p = np.array([1.0, 2.0, 3.0, 4.0]) / np.sqrt(30.0)
    tangent = np.array([2.0, -1.0, 4.0, -3.0]) / np.sqrt(30.0)
    theta = 0.7
    end = math.cos(theta) * p + math.sin(theta) * tangent
    for t in (0.0, 0.2, 0.8, 1.0):
        expected = math.cos(t * theta) * p + math.sin(t * theta) * tangent
        close(rotlib.quat_interpolate(p, -end, t), expected)


@pytest.mark.parametrize("axis, angle", CASES)
@pytest.mark.parametrize("world", [False, True])
def test_quaternion_derivative_matrices_and_wrappers(axis, angle, world):
    q = axis_quaternion(axis, angle)
    conjugate = q.copy()
    conjugate[1:] *= -1
    expected_from = np.column_stack(
        [
            0.5 * (hamilton(v, q) if world else hamilton(q, v))
            for v in np.eye(4)[1:]
        ]
    )
    expected_to = np.column_stack(
        [
            2
            * (hamilton(v, conjugate) if world else hamilton(conjugate, v))[1:]
            for v in np.eye(4)
        ]
    )
    close(rotlib.quat_deriv_from_angvel_mat(q, world=world), expected_from)
    close(rotlib.quat_deriv_to_angvel_mat(q, world=world), expected_to)
    for omega in (np.zeros(3), np.array([0.3, -0.4, 0.7])):
        qdot = rotlib.quat_deriv_from_angvel(q, omega, world=world)
        close(qdot, expected_from @ omega)
        close(np.dot(q, qdot), 0.0)
        close(rotlib.quat_deriv_to_angvel(q, qdot, world=world), omega)
        close(rotlib.quat_deriv_from_angvel(-q, omega, world=world), -qdot)


@pytest.mark.parametrize("axis, angle", CASES)
def test_quat_rotmat_and_shiftmat(axis, angle):
    q = axis_quaternion(axis, angle)
    expected = quaternion_matrix(q)
    close(rotlib.quat_rotmat(q), expected)
    close(rotlib.quat_rotmat(-q), expected)
    close(rotlib.quat_shiftmat(q, forward=True), expected.T)
    close(rotlib.quat_shiftmat(q), expected)


@pytest.mark.parametrize("v", LAYOUTS)
def test_quat_rotate_vectors_and_shift_vectors(v):
    q = axis_quaternion(MIXED_AXIS, 0.7)
    r = quaternion_matrix(q)
    out = rotlib.quat_rotate_vectors(v, q)
    assert out.shape == v.shape
    close(out, v @ r.T)
    for forward in (False, True):
        s = r.T if forward else r
        out = rotlib.quat_shift_vectors(v, q, forward=forward)
        assert out.shape == v.shape
        close(out, v @ s.T)


@pytest.mark.parametrize(
    "rank, function",
    [
        (2, rotlib.quat_shift_tensor2),
        (3, rotlib.quat_shift_tensor3),
    ],
)
@pytest.mark.parametrize("forward", [False, True])
def test_quat_shift_tensors(rank, function, forward):
    q = axis_quaternion(MIXED_AXIS, 0.7)
    r = quaternion_matrix(q)
    a = np.arange(3**rank, dtype=np.float64).reshape((3,) * rank) - 4.0
    close(
        function(a, q, forward=forward),
        tensor_reference(a, r.T if forward else r),
    )


@pytest.mark.parametrize("sign", [-1.0, 1.0])
@pytest.mark.parametrize("angle", [0.0, 1e-200, 0.7, np.pi])
def test_quat_to_aa(angle, sign):
    q = sign * axis_quaternion(MIXED_AXIS, angle)
    axis, actual = rotlib.quat_to_aa(q)
    check_axis_angle(axis, actual, quaternion_matrix(q))
    assert_allclose(actual, angle, rtol=2e-14, atol=0.0)
    if angle == 0.0:
        close(axis, [1.0, 0.0, 0.0], atol=0.0)


@pytest.mark.parametrize("axis, angle", CASES)
def test_quat_to_dcm(axis, angle):
    q = axis_quaternion(axis, angle)
    close(rotlib.quat_to_dcm(q), quaternion_matrix(q).T)


@pytest.mark.parametrize("seq, world", MODES)
def test_quat_to_euler(seq, world):
    q = axis_quaternion(MIXED_AXIS, 0.7)
    close(
        euler_matrix(rotlib.quat_to_euler(q, seq, world), seq, world),
        quaternion_matrix(q),
    )


def test_quat_to_euler_tolerance_and_default():
    e = np.array([0.7, np.pi / 2 - 5e-5, -0.4])
    q = euler_quaternion(e, "XYZ", True)
    close(
        euler_matrix(rotlib.quat_to_euler(q), "XYZ", True),
        euler_matrix(e, "XYZ", True),
    )
    assert rotlib.quat_to_euler(q, tol=1e-4)[0] == 0.0


@pytest.mark.parametrize(
    "representation", ["quat", "euler", "axis_angle", "dcm"]
)
def test_quat_from_any(representation):
    payload, matrix = orientation_payload(representation)
    q = rotlib.quat_from_any(payload)
    close(math.hypot(*q), 1.0)
    close(quaternion_matrix(q), matrix)


def test_quat_from_any_preserves_existing_sign():
    q = -axis_quaternion(MIXED_AXIS, 0.7)
    close(rotlib.quat_from_any({"repr": "quat", "quat": q}), q, atol=0.0)


def test_quat_from_any_unknown_representation():
    with pytest.raises(ValueError, match="Unrecognized orientation repr"):
        rotlib.quat_from_any({"repr": "unknown"})
