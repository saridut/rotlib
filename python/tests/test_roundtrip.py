"""Selected roundtrips adapted from xtest_main.py for the current API.

Compare orientations rather than unique parameters: q and -q are equivalent,
Euler triples need not be recovered, and an identity rotation has no unique
axis. Independent intermediate checks guard against matching inverse bugs.
Inputs are finite float64 arrays of the documented shapes.
"""

import math

import numpy as np
import pytest
from _rotation_cases import (
    MIXED_AXIS,
    VECTORS,
    axis_matrix,
    axis_quaternion,
    close,
    euler_matrix,
    euler_quaternion,
    quaternion_matrix,
    tensor_reference,
)
from numpy.testing import assert_allclose

import rotlib

# Exercise each possible middle angle and both fixed/moving-axis conventions
# without multiplying every conversion by every possible sequence.
EULER_MODES = [("XYZ", True), ("XZY", False), ("ZXY", True)]


@pytest.mark.parametrize("via", ["axis-angle", "dcm"])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
@pytest.mark.parametrize("angle", [0.0, 1e-200, 0.7, np.pi - 1e-10, np.pi])
def test_quaternion_roundtrip(via, sign, angle):
    q = sign * axis_quaternion(MIXED_AXIS, angle)
    expected = axis_matrix(MIXED_AXIS, angle)
    if via == "axis-angle":
        axis, recovered_angle = rotlib.quat_to_aa(q)
        close(axis_matrix(axis, recovered_angle), expected)
        back = rotlib.aa_to_quat(axis, recovered_angle)
    else:
        dcm = rotlib.quat_to_dcm(q)
        close(dcm, expected.T)
        back = rotlib.dcm_to_quat(dcm)
    if np.dot(q, back) < 0.0:
        back = -back
    close(back, q)
    close(quaternion_matrix(back), expected)
    if angle == 1e-200:
        # An absolute matrix tolerance alone would also accept lost rotation.
        assert_allclose(back[1:], q[1:], rtol=2e-14, atol=0.0)


@pytest.mark.parametrize("angle", [-0.7, 4.0, 2 * np.pi])
def test_axis_angle_dcm_roundtrip_canonicalizes_parameters(angle):
    dcm = rotlib.aa_to_dcm(MIXED_AXIS, angle)
    close(dcm, axis_matrix(MIXED_AXIS, angle).T)
    axis, recovered_angle = rotlib.dcm_to_aa(dcm)
    assert 0.0 <= recovered_angle <= np.pi
    close(math.hypot(*axis), 1.0)
    close(rotlib.aa_to_dcm(axis, recovered_angle), dcm)


@pytest.mark.parametrize("seq, world", EULER_MODES)
@pytest.mark.parametrize("via", ["quat", "dcm"])
@pytest.mark.parametrize(
    "middle",
    [None, np.pi / 2, -np.pi / 2, np.pi / 2 - 5e-15, np.pi / 2 + 1e-8],
)
def test_euler_representation_roundtrip(seq, world, via, middle):
    # Large outer angles intentionally start outside the principal ranges.
    original = np.array([3.8, -4.2, 7.1])
    if middle is not None:
        original["XYZ".index(seq[1])] = middle
    expected = euler_matrix(original, seq, world)
    if via == "quat":
        q = rotlib.euler_to_quat(original, seq, world)
        close(quaternion_matrix(q), expected)
        back = rotlib.quat_to_euler(q, seq, world)
    else:
        dcm = rotlib.euler_to_dcm(original, seq, world)
        close(dcm, expected.T)
        back = rotlib.dcm_to_euler(dcm, seq, world)
    # The default 1e-14 cutoff may discard a small outer-angle contribution.
    # This is an absolute error in dimensionless matrix entries, not angles.
    close(euler_matrix(back, seq, world), expected, atol=3e-14)


@pytest.mark.parametrize("seq, world", EULER_MODES)
@pytest.mark.parametrize("tol", [1e-6, 1e-4])
def test_euler_roundtrip_deliberate_near_lock_approximation(seq, world, tol):
    angles = np.array([0.7, -0.4, 0.3])
    angles["XYZ".index(seq[1])] = np.pi / 2 - tol / 2
    matrix = euler_matrix(angles, seq, world)
    back = rotlib.euler_factor_rotmat(matrix, seq, world, tol=tol)
    assert back["XYZ".index(seq[0])] == 0.0
    rebuilt = rotlib.euler_rotmat(back, seq, world)
    # For these inputs the discarded contribution is on the scale of tol.
    # Do not impose machine precision, or claim a universal tol error bound.
    error = np.max(np.abs(rebuilt - matrix))
    assert 1e-3 * tol < error < 2 * tol
    # Turning off the cutoff retains the small, resolvable contribution.
    uncollapsed = rotlib.euler_factor_rotmat(matrix, seq, world, tol=0.0)
    close(euler_matrix(uncollapsed, seq, world), matrix)


@pytest.mark.parametrize(
    "source, target",
    [
        (("XYZ", True), ("YZX", False)),
        (("XZY", False), ("ZXY", True)),
        (("ZYX", True), ("YXZ", False)),
    ],
)
@pytest.mark.parametrize("middle", [0.4, np.pi / 2, -np.pi / 2])
def test_euler_sequence_and_frame_roundtrip(source, target, middle):
    angles = np.array([0.3, -0.5, 0.7])
    angles["XYZ".index(source[0][1])] = middle
    expected = euler_matrix(angles, *source)
    converted = rotlib.euler_to_euler(angles, *source, *target)
    close(euler_matrix(converted, *target), expected)
    back = rotlib.euler_to_euler(converted, *target, *source)
    close(euler_matrix(back, *source), expected)


@pytest.mark.parametrize("representation", ["quat", "aa", "dcm", "euler"])
@pytest.mark.parametrize("rank", [1, 2, 3])
def test_frame_shift_roundtrip(representation, rank):
    angles = np.array([0.2, -0.3, 0.4])
    seq, world = "YZX", False
    matrix = euler_matrix(angles, seq, world)
    q = euler_quaternion(angles, seq, world)
    axis, angle = rotlib.quat_to_aa(q)
    args = {
        "quat": (q,),
        "aa": (axis, angle),
        "dcm": (matrix.T,),
        "euler": (angles, seq, world),
    }[representation]
    suffix = "vectors" if rank == 1 else f"tensor{rank}"
    shift = getattr(rotlib, f"{representation}_shift_{suffix}")
    original = (
        VECTORS.copy()
        if rank == 1
        else np.arange(3**rank, dtype=np.float64).reshape((3,) * rank) - 4
    )
    expected = (
        original @ matrix
        if rank == 1
        else tensor_reference(original, matrix.T)
    )
    # Tensor contractions accumulate error across indices; scale the absolute
    # budget by input magnitude so entries near zero are treated sensibly too.
    budget = 32 * rank * np.finfo(np.float64).eps * np.max(np.abs(original))
    in_b = shift(original.copy(), *args, forward=True)
    close(in_b, expected, atol=budget)
    close(shift(in_b.copy(), *args, forward=False), original, atol=budget)


@pytest.mark.parametrize("world", [False, True])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_quaternion_derivative_roundtrip_across_resolution_frames(world, sign):
    q = sign * axis_quaternion(MIXED_AXIS, 0.9)
    matrix = quaternion_matrix(q)
    omega = np.array([0.3, -0.7, 1.1])
    qdot = rotlib.quat_deriv_from_angvel(q, omega, world=world)
    close(np.dot(q, qdot), 0.0)
    close(rotlib.quat_deriv_to_angvel(q, qdot, world=world), omega)
    other = rotlib.quat_deriv_to_angvel(q, qdot, world=not world)
    close(other, matrix.T @ omega if world else matrix @ omega)
    close(rotlib.quat_deriv_from_angvel(q, other, world=not world), qdot)


@pytest.mark.parametrize("world", [False, True])
def test_tangent_quaternion_derivative_roundtrip(world):
    q = axis_quaternion(MIXED_AXIS, 0.9)
    candidate = np.array([0.7, -1.2, 0.4, 2.3])
    # A unit quaternion's derivative must be tangent to the unit sphere.
    qdot = candidate - np.dot(candidate, q) * q
    close(np.dot(q, qdot), 0.0)
    omega = rotlib.quat_deriv_to_angvel(q, qdot, world=world)
    close(rotlib.quat_deriv_from_angvel(q, omega, world=world), qdot)
