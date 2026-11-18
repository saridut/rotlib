"""Independent references shared by the public-API tests.

All inputs are finite float64 arrays of documented shapes. These helpers
never call Rotlib. Euler arrays always store X, Y, Z angles in that order.
"""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

ATOL = 2e-14
SEQUENCES = ("XYZ", "XZY", "YXZ", "YZX", "ZXY", "ZYX")
MODES = [(seq, world) for seq in SEQUENCES for world in (True, False)]
VECTORS = np.array([[2.0, -3.0, 4.0], [0.0, 0.0, 0.0], [-1.0, 0.5, 2.0]])
LAYOUTS = [
    pytest.param(VECTORS[0], id="vector"),
    pytest.param(VECTORS[:1], id="single-row"),
    pytest.param(VECTORS, id="batch"),
]
MIXED_AXIS = np.array([-1.0, 2.0, -3.0]) / np.sqrt(14.0)
CASES = [
    pytest.param(np.array([1.0, 0.0, 0.0]), 0.0, id="identity"),
    pytest.param(np.array([0.0, 0.0, 1.0]), np.pi / 2, id="quarter-turn-z"),
    pytest.param(MIXED_AXIS, 0.7, id="mixed-axis"),
    pytest.param(MIXED_AXIS, np.pi, id="mixed-half-turn"),
]


def close(actual, expected, atol=ATOL):
    assert_allclose(actual, expected, rtol=0.0, atol=atol)


def hamilton(p, q):
    """Scalar/vector definition, independent of Rotlib's product matrix."""
    return np.concatenate(
        (
            [p[0] * q[0] - np.dot(p[1:], q[1:])],
            p[0] * q[1:] + q[0] * p[1:] + np.cross(p[1:], q[1:]),
        )
    )


def axis_quaternion(axis, angle):
    return np.concatenate(([math.cos(angle / 2)], axis * math.sin(angle / 2)))


def quaternion_matrix(q):
    """Rotate the coordinate basis using quaternion sandwich products."""
    conjugate = q.copy()
    conjugate[1:] *= -1
    return np.column_stack(
        [
            hamilton(hamilton(q, np.concatenate(([0.0], v))), conjugate)[1:]
            for v in np.eye(3)
        ]
    )


def axis_matrix(axis, angle):
    return quaternion_matrix(axis_quaternion(axis, angle))


def elementary(k, angle):
    c, s = math.cos(angle), math.sin(angle)
    if k == 0:
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    if k == 1:
        return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def euler_matrix(angles, seq, world):
    result = np.eye(3)
    for letter in seq:
        k = "XYZ".index(letter)
        step = elementary(k, angles[k])
        result = step @ result if world else result @ step
    return result


def euler_quaternion(angles, seq, world):
    result = np.array([1.0, 0.0, 0.0, 0.0])
    for letter in seq:
        k = "XYZ".index(letter)
        step = axis_quaternion(np.eye(3)[k], angles[k])
        result = hamilton(step, result) if world else hamilton(result, step)
    return result / math.hypot(*result)


def tensor_reference(a, transform):
    """Transform a tensor as a sum of transformed basis tensor products."""
    result = np.zeros_like(a)
    for indices in np.ndindex(a.shape):
        term = transform[:, indices[0]]
        for k in indices[1:]:
            term = np.multiply.outer(term, transform[:, k])
        result += a[indices] * term
    return result


def check_axis_angle(axis, angle, expected):
    close(math.hypot(*axis), 1.0)
    assert 0.0 <= angle <= np.pi
    close(axis_matrix(axis, angle), expected)


def check_canonical(q):
    assert q[np.flatnonzero(q)[0]] > 0.0


def orientation_payload(representation):
    """Four descriptions of one ordinary, noncommuting orientation."""
    angles = np.array([0.2, -0.3, 0.4])
    seq, world = "YZX", False
    q = euler_quaternion(angles, seq, world)
    angle = 2 * math.atan2(math.hypot(*q[1:]), q[0])
    axis = q[1:] / math.hypot(*q[1:])
    matrix = euler_matrix(angles, seq, world)
    payloads = {
        "quat": {"quat": q},
        "euler": {"euler": angles, "seq": seq, "world": world},
        "axis_angle": {"axis": axis, "angle": angle},
        "dcm": {"dcm": matrix.T},
    }
    return {"repr": representation, **payloads[representation]}, matrix
