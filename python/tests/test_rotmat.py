"""Representative mathematical cases for rotation matrices.
"""

import numpy as np
import pytest

import rotlib

@pytest.mark.parametrize(
    "mat, expected", 
    [
        #Zero matrix
        pytest.param(np.zeros((3,3)), False, id="zero"),
        #Identity matrix
        pytest.param(np.eye(3), True, id="identity"),
        #Negative identity matrix is invalid
        pytest.param(-np.eye(3), False, id="-ve identity"),
        #180-degree rotation: Negative diagonal entries are valid
        pytest.param(np.diag([1, -1, -1]), True, id="rot_180"),
        #Half-turn (180 deg rotation) about [1,1,1]/sqrt(3)
        pytest.param(
            np.array( [[-1.0, 2.0, 2.0],
                       [2.0, -1.0, 2.0],
                       [2.0, 2.0, -1.0]] )/3.0,
            True, id="dense_half_turn"),
        # Orthogonal with determinant -1
        pytest.param(np.diag([-1.0, 1.0, 1.0]), False, id="reflection"),
        # Diagonal with determinant 1
        pytest.param(np.diag([2.0, 0.5, 1.0]), False, id="diag_det_1"),
        # Non-orthogonal with determinant +1
        pytest.param(
            np.array([[1.0, 0.2, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]]),
            False, id="determinant_one_shear"),
    ]
)
def test_mat_is_rotmat(mat, expected):
    assert rotlib.mat_is_rotmat(mat) == expected


@pytest.mark.parametrize(
    "eps, expected",
    [
        pytest.param(0.5, True, id="within_orthogonality_tol"),
        pytest.param(2, False, id="beyond_orthogonality_tol"),
    ],
)
def test_mat_is_rotmat_orthogonality_tol(eps, expected):
    # Choosing perturbations on either side of the cutoff.
    # Isolates the orthogonality criterion.
    tol = 1e-12
    delta = eps*tol
    mat = np.array([[1.0, delta, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0]])
    assert rotlib.mat_is_rotmat(mat, tol) == expected


@pytest.mark.parametrize(
    "eps, expected",
    [
        pytest.param(0.2, True, id="within_determinant_tol"),
        pytest.param(0.4, False, id="beyond_determinant_tol"),
    ],
)
def test_mat_is_rotmat_determinant_tol(eps, expected):
    # Choosing perturbations on either side of the cutoff.
    # Isolates the determinant criterion.
    tol = 1e-12
    delta = eps*tol
    mat = (1+delta)*np.eye(3)
    assert rotlib.mat_is_rotmat(mat, tol) == expected


# mat_is_dcm is the public DCM-validation wrapper. Keep these focused on its
# accepted matrices and tolerance forwarding; mat_is_rotmat is tested above.
@pytest.mark.parametrize("mat, expected", [
    pytest.param(np.eye(3), True, id="dcm-identity"),
    pytest.param(np.array([[-1.0, 2.0, 2.0],
                           [2.0, -1.0, 2.0],
                           [2.0, 2.0, -1.0]]) / 3, True, id="dcm-half-turn"),
    pytest.param(np.diag([-1.0, 1.0, 1.0]), False, id="dcm-reflection"),
    pytest.param(np.array([[1.0, 0.2, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]]), False, id="dcm-shear"),
])
def test_mat_is_dcm(mat, expected):
    assert rotlib.mat_is_dcm(mat) == expected


def test_mat_is_dcm_forwards_tolerance():
    mat = (1.0 + 1e-10) * np.eye(3)
    assert not rotlib.mat_is_dcm(mat)
    assert rotlib.mat_is_dcm(mat, tol=1e-8)
