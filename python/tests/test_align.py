"""Active alignment in a shared fixed frame, using valid float64 inputs.

Expected rotations are specified independently of Rotlib's conversion
functions. No type, nonfinite-value, or input-mutation checks are made.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import rotlib


ATOL = 1e-14
IDENTITY = np.eye(3)
QUARTER_TURN_Z = np.array(
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
)
HALF_TURN_MIXED = np.array(
    [[-1.0, 2.0, 2.0], [2.0, -1.0, 2.0], [2.0, 2.0, -1.0]]
) / 3.0
# A right-handed orthonormal frame distinct from the coordinate basis.
OLD_FRAME = np.array(
    [[1.0, 2.0, 2.0], [2.0, 1.0, -2.0], [-2.0, 2.0, -1.0]]
) / 3.0
VECTORS = np.array(
    [[2.0, -3.0, 4.0], [0.0, 0.0, 0.0], [-1.0, 5.0, 0.5]]
)


def _assert_close(actual, expected):
    assert_allclose(actual, expected, rtol=0.0, atol=ATOL)


def _assert_proper_rotation(rotated_basis):
    # Rotated coordinate vectors are rows, so this array equals R.T.
    # Orthogonality alone would also accept a reflection such as -I.
    _assert_close(rotated_basis @ rotated_basis.T, IDENTITY)
    _assert_close(np.linalg.det(rotated_basis), 1.0)
    _assert_close(
        np.cross(rotated_basis[0], rotated_basis[1]), rotated_basis[2]
    )


class TestAlignOneReference:
    @pytest.mark.parametrize("ref_shape", [(3,), (1, 3), (3, 1)])
    @pytest.mark.parametrize(
        "v",
        [
            pytest.param(VECTORS[0], id="vector"),
            pytest.param(VECTORS[:1], id="single-row"),
            pytest.param(VECTORS, id="batch"),
        ],
    )
    def test_shortest_quarter_turn(self, ref_shape, v):
        # Both references lie in xy: the shortest rotation is +pi/2 about z.
        old = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
        new = np.array([-1.0, 1.0, 0.0]) / np.sqrt(2.0)
        old = old.reshape(ref_shape)
        new = new.reshape(ref_shape)

        out = rotlib.align(v, old, new)

        assert out.shape == v.shape
        _assert_close(out, v @ QUARTER_TURN_Z.T)
        _assert_close(rotlib.align(old.ravel(), old, new), new.ravel())
        # The perpendicular axis must remain fixed: no additional twist.
        _assert_close(rotlib.align(IDENTITY[2], old, new), IDENTITY[2])

    def test_parallel_references_leave_vectors_unchanged(self):
        old = np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)

        out = rotlib.align(VECTORS, old, old.copy())

        assert out.shape == VECTORS.shape
        _assert_close(out, VECTORS)

    @pytest.mark.parametrize(
        "old",
        [
            pytest.param(IDENTITY[0], id="x-axis"),
            pytest.param(IDENTITY[1], id="y-axis"),
            pytest.param(IDENTITY[2], id="z-axis"),
            pytest.param(
                np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0), id="mixed-axis"
            ),
        ],
    )
    def test_antiparallel_references_give_proper_half_turn(self, old):
        new = -old
        rotated_basis = rotlib.align(IDENTITY, old, new)
        out = rotlib.align(VECTORS, old, new)

        _assert_close(rotlib.align(old, old, new), new)
        _assert_proper_rotation(rotated_basis)
        _assert_close(np.trace(rotated_basis), -1.0)
        # Do not prescribe which perpendicular half-turn axis is selected.
        _assert_close(out, VECTORS @ rotated_basis)
        _assert_close(out @ out.T, VECTORS @ VECTORS.T)
        assert out.shape == VECTORS.shape

    @pytest.mark.parametrize(
        "angle",
        [
            pytest.param(1e-8, id="nearly-parallel"),
            pytest.param(np.pi - 1e-8, id="nearly-antiparallel"),
        ],
    )
    def test_nearly_collinear_references(self, angle):
        # These directions are close to collinear but still resolvable.
        c, s = np.cos(angle), np.sin(angle)
        old = IDENTITY[0]
        new = np.array([c, s, 0.0])
        expected_rotation = np.array(
            [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
        )

        rotated_basis = rotlib.align(IDENTITY, old, new)

        _assert_close(rotated_basis, expected_rotation.T)
        _assert_close(rotlib.align(old, old, new), new)
        _assert_proper_rotation(rotated_basis)


# Distribute vector layouts across known rotations without testing every
# combination. The dense starting frame exposes multiplication-order errors.
FRAME_CASES = [
    pytest.param(IDENTITY, VECTORS, id="unchanged-frame"),
    pytest.param(QUARTER_TURN_Z, VECTORS[0], id="quarter-turn-vector"),
    pytest.param(
        np.diag([-1.0, -1.0, 1.0]), VECTORS[:1], id="half-turn-single-row"
    ),
    pytest.param(HALF_TURN_MIXED, VECTORS, id="mixed-axis-half-turn-batch"),
]


def _check_frame_alignment(reference_count, rotation, v):
    new_frame = OLD_FRAME @ rotation.T
    old = OLD_FRAME[:reference_count].copy()
    new = new_frame[:reference_count].copy()

    out = rotlib.align(v, old, new)

    assert out.shape == v.shape
    _assert_close(out, v @ rotation.T)
    # Includes the completed third direction in the two-reference case.
    _assert_close(rotlib.align(OLD_FRAME, old, new), new_frame)
    rotated_basis = rotlib.align(IDENTITY, old, new)
    _assert_proper_rotation(rotated_basis)


class TestAlignTwoReferences:
    @pytest.mark.parametrize("rotation, v", FRAME_CASES)
    def test_known_rotation(self, rotation, v):
        _check_frame_alignment(2, rotation, v)


class TestAlignThreeReferences:
    @pytest.mark.parametrize("rotation, v", FRAME_CASES)
    def test_known_rotation(self, rotation, v):
        _check_frame_alignment(3, rotation, v)


def test_align_rejects_mismatched_reference_shapes():
    # Exercise the existing explicit rejection, without extending validation
    # to other inputs outside the documented preconditions.
    with pytest.raises(ValueError, match="differs from"):
        rotlib.align(VECTORS, IDENTITY[:2], IDENTITY)
