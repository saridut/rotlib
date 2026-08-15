"""
**Euler angle ranges**

=============   ===============   ================  ================
 Sequence            phi               theta           psi
=============   ===============   ================  ================
  XYZ,   ZYX     [-pi,   pi]        [-pi/2, pi/2]     [-pi, pi]
  XZY,   YZX     [-pi,   pi]        [-pi,  pi]        [-pi/2, pi/2]
  ZXY,   YXZ     [-pi/2, pi/2]      [-pi,  pi]        [-pi, pi]
=============   ===============   ================  ================

Euler angle sequence: 'XYZ' (world). First rotation by *phi* about X, second
rotation by *theta* about Y, and the third rotation by *psi* about Z axis of
the world (i.e. fixed) frame.

This is the same as the sequence used in the `Blender blenlib code
<https://github.com/blender/blender/blob/
7d641fe96810cdc2598b2f37ec4f6192e133e6d8/
source/blender/blenlib/BLI_math_euler_types.hh>`__.

In contrast, the 'XYZ' sequence is understood in the Aerospace community as:
First rotation about Z-axis, second rotation about Y-axis, and the third
rotation about X-axis of the body frame.

*Reference:* D. Eberly, `Euler angle formulas.
<http://www.geometrictools.com/Documentation/EulerAngles.pdf>`__

"""

import math

import numpy as np

# *****************************************************************************#
#                               UTILITY FUNCTIONS                             #
# *****************************************************************************#


def align(v, old, new):
    """
    Rotate vectors such that a set of mutually orthogonal unit vectors aligns
    with another set of mutually orthogoanl unit vectors.

    The vectors `v` are rotated using a rotation matrix that, when applied to
    the unit vectors `old`, aligns them in the direction of the unit vectors
    `new`.

    Parameters
    ----------
    v : (3,) | (n,3) ndarray
        Vectors to align.
    old : (3,) | (m,3) ndarray
        A set of at most three mutually orthogonal unit vectors. A single
        vector can be of shape (3,) or (1,3). Else `m` must be 2 or 3.
    new : (3,) | (m,3) ndarray
        Another set of mutually orthogonal unit vectors of the same
        shape as `old`. The angles between the vectors must be the same as
        those of `old`.

    Returns
    -------
    out : ndarray
        Aligned vectors. `out` is of the same shape as `v`.
    """
    if np.any(old.shape != new.shape):
        raise ValueError(
            f"`old` of shape {old.shape} differs from`new` of shape {new.shape}."
        )
    n = old.size // 3
    if n == 1:
        # Angle between old and new
        cos_angle = np.dot(old, new)
        if (cos_angle > 1.0) or (cos_angle < -1.0):
            if np.isclose(abs(cos_angle), 1.0, rtol=1e-9, atol=1e-15):
                cos_angle = math.copysign(1.0, cos_angle)
            else:
                raise ValueError(
                    f"Cosine of angle(= {cos_angle}) between"
                    f" `old`(= {old}) and `new`(= {new}) out-of-range."
                )
        if np.isclose(cos_angle, 1.0, rtol=1e-9, atol=1e-15):
            # Angle = zero: No rotation needed
            out = v
        elif np.isclose(cos_angle, -1.0, rtol=1e-9, atol=1e-15):
            # Angle = pi: Rotation axis not unique. We flip the
            # vectors `v` to align.
            out = -v
        else:
            angle = math.acos(cos_angle)
            # Axis of rotation
            axis = np.cross(old, new)
            axis_nrm = np.linalg.norm(axis)
            axis /= axis_nrm
            out = aa_rotate_vectors(v, axis, angle)
        return out
    elif n == 2:
        # The third vector orthogonal to `old`
        z_old = np.cross(old[0, :], old[1, :])
        z_old_nrm = np.linalg.norm(z_old)
        z_old /= z_old_nrm
        # The third vector orthogonal to `new`
        z_new = np.cross(new[0, :], new[1, :])
        z_new_nrm = np.linalg.norm(z_new)
        z_new /= z_new_nrm
        # The `old` and new `frames`
        axes_old = np.vstack((old, z_old))
        axes_new = np.vstack((new, z_new))
        dcm = dcm_from_axes(axes_old, axes_new)
        return dcm_rotate_vectors(v, dcm)
    elif n == 3:
        dcm = dcm_from_axes(old, new)
        return dcm_rotate_vectors(v, dcm)


def mat_is_rotmat(mat):
    """
    Checks if `mat` is a rotation matrix or not.

    Parameters
    ----------
    mat : (3,3) ndarray
        Array to check.

    Returns
    -------
    bool
        ``True`` if `mat` is a rotation matrix, ``False`` otherwise.
    """
    det_is_one = math.isclose(
        np.linalg.det(mat), 1.0, abs_tol=1e-12, rel_tol=1e-12
    )
    is_orthogonal = np.allclose(np.dot(mat, mat.T), np.identity(3))
    return is_orthogonal and det_is_one


# QUATERNION-----------------------------------------------------------


def quat_rand(rng):
    """
    Returns a random unit quaternion.

    Parameters
    ----------
    rng : :py:class:`numpy.random.Generator`
        A random number generator.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    axis, angle = aa_rand(rng)
    q = aa_to_quat(axis, angle)
    return q


def quat_identity():
    """
    Returns the identity unit quaternion.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    return np.array([1.0, 0.0, 0.0, 0.0])


def quat_conjugated(q):
    """
    Conjugates a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Conjugated quaternion.
    """
    q[1:4] = -q[1:4]
    return q


def quat_inverted(q):
    """
    Inverts a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Inverted quaternion.
    """
    quat_conjugated(q)
    nrm = np.linalg.norm(q)
    q /= nrm
    return q


def quat_normalized(q):
    """
    Normalizes a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Normalized quaternion.
    """
    nrm = np.linalg.norm(q)
    q /= nrm
    return q


def quat_is_normalized(q):
    """
    Checks whether a quaternion is normalized, i.e. whether it is a unit
    quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    bool
        ``True`` if `q` is normalized, ``False`` otherwise.
    """
    norm = np.linalg.norm(q)
    return math.isclose(norm, 1.0, rel_tol=1e-14)


def quat_prod(p, q, normalize=True):
    """
    Returns the product of two quaternions.

    Parameters
    ----------
    p : (4,) ndarray
        Quaternion.
    q : (4,) ndarray
        Quaternion.
    normalize : bool
        Whether to normalize the product.

    Returns
    -------
    (4,) ndarray
        Product of two quaternions. If `normalize` is ``True``, this is also
        a unit quaternion.

    """
    p0, p1, p2, p3 = tuple(p)
    prod_mat = np.array([[p0, -p1, -p2, -p3],
                         [p1,  p0, -p3,  p2],
                         [p2,  p3,  p0, -p1],
                         [p3, -p2,  p1,  p0]])  # fmt: skip
    pq = np.dot(prod_mat, q)
    if normalize:
        quat_normalized(pq)
    return pq


def quat_angle_between(p, q):
    """
    Returns the angle between two unit quaternions p and q.

    Parameters
    ----------
    p : (4,) ndarray
        Unit quaternion.
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    float
        Angle in radian.
    """
    cos_angle = np.dot(p, q)
    if (cos_angle > 1.0) or (cos_angle < -1.0):
        if np.isclose(abs(cos_angle), 1.0, rtol=1e-9, atol=1e-15):
            cos_angle = math.copysign(1.0, cos_angle)
        else:
            raise ValueError(
                f"Cosine of angle(= {cos_angle}) between"
                f" `p`(= {p}) and `q`(= {q}) out-of-range."
            )
    angle = np.acos(cos_angle)
    return angle


def quat_interpolate(q1, q2, t):
    """
    Interpolate between two unit quaternions.

    Parameters
    ----------
    q1 : (4,) ndarray
        Unit quaternion.
    q2 : (4,) ndarray
        Unit quaternion.
    t : float
        A fraction between 0 and 1 (both inclusive) specifying the
        interpolation point.

    Returns
    -------
    (4,) ndarray
        Interpolated unit quaternion.
    """
    theta = quat_angle_between(q1, q2)
    q = (
        q1 * math.sin((1.0 - t) * theta) + q2 * math.sin(t * theta)
    ) / math.sin(theta)
    return quat_normalized(q)


def quat_deriv_to_angvel_mat(q):
    """
    Returns the matrix mapping the derivative of a unit quaternion to angular
    velocity.

    Parameters
    ----------
    qdot : (4,) ndarray
        Derivative of a unit quaternion.

    Returns
    -------
    (3,4) ndarray
        Angular velocity matrix.
    """
    q0, q1, q2, q3 = tuple(q)
    return 2*np.array([[-q1,  q0,  q3, -q2],
                       [-q2, -q3,  q0,  q1],
                       [-q3,  q2, -q1,  q0]])  # fmt: skip


def quat_deriv_to_angvel(q, qdot):
    """
    Calculates the angular velocity from a unit quaternion and its time
    derivative.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    qdot : (4,) ndarray
        Derivative of `q`. In general, this is not a unit vector.

    Returns
    -------
    (3,) ndarray
        Angular velocity.
    """
    mat = quat_deriv_to_angvel_mat(q)
    return np.dot(mat, qdot)


def quat_deriv_from_angvel_mat(q):
    """
    Returns the matrix mapping angular velocity to time derivative of a unit
    quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (4,3) ndarray
        Quaternion derivative matrix.
    """
    q0, q1, q2, q3 = tuple(q)
    return 0.5*np.array([[-q1, -q2, -q3],
                         [ q0, -q3,  q2],
                         [ q3,  q0, -q1],
                         [-q2,  q1,  q0]])  # fmt: skip


def quat_deriv_from_angvel(q, ang_vel):
    """
    Calculates the time derivative of a unit quaternion from angular velocity.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    ang_vel : (3,) ndarray
        Angular velocity.

    Returns
    -------
    (3,) ndarray
        Derivative of `q`. In general, this is not a unit vector.
    """
    mat = quat_deriv_from_angvel_mat(q)
    qdot = np.dot(mat, ang_vel)
    return qdot


def quat_rotmat(q):
    """Returns the rotation matrix corresponding to a unit quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.
    """
    rotmat = np.empty((3, 3))

    q0sq = q[0] * q[0]
    q1sq = q[1] * q[1]
    q2sq = q[2] * q[2]
    q3sq = q[3] * q[3]
    q0q1 = q[0] * q[1]
    q0q2 = q[0] * q[2]
    q0q3 = q[0] * q[3]
    q1q2 = q[1] * q[2]
    q1q3 = q[1] * q[3]
    q2q3 = q[2] * q[3]

    rotmat[0, 0] = 2 * (q0sq + q1sq) - 1.0
    rotmat[0, 1] = 2 * (q1q2 - q0q3)
    rotmat[0, 2] = 2 * (q1q3 + q0q2)
    rotmat[1, 0] = 2 * (q1q2 + q0q3)
    rotmat[1, 1] = 2 * (q0sq + q2sq) - 1.0
    rotmat[1, 2] = 2 * (q2q3 - q0q1)
    rotmat[2, 0] = 2 * (q1q3 - q0q2)
    rotmat[2, 1] = 2 * (q2q3 + q0q1)
    rotmat[2, 2] = 2 * (q0sq + q3sq) - 1.0
    return rotmat


def quat_shiftmat(q, forward=False):
    """
    Returns the shifter matrix corresponding to a unit quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
       Whether to shift forward, i.e., along the orientation or shift reverse.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """
    if forward:
        conj_q = q.copy()
        shiftmat = quat_rotmat(quat_conjugated(conj_q))
    else:
        shiftmat = quat_rotmat(q)
    return shiftmat


def quat_rotate_vectors(v, q):
    """
    Rotates vectors by a unit quaternion.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.
    """
    rotmat = quat_rotmat(q)
    return np.dot(v, rotmat.T)


def quat_shift_vectors(v, q, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a unit quaternion, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = quat_shiftmat(q, forward=forward)
    return np.dot(v, shiftmat.T)


def quat_shift_tensor2(a, quat, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a unit quaternion, shifts second-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor.
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second-order tensor.
    """
    shiftmat = quat_shiftmat(quat, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def shift_tensor3_quat(a, quat, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a unit quaternion, shifts third-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor.
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third-order tensor.
    """
    shiftmat = quat_shiftmat(quat, forward=forward)
    return np.einsum("ip,jq,kr,pqr", shiftmat, shiftmat, shiftmat, a)


def quat_to_aa(q):
    """
    Converts a unit quaternion to an *axis-angle* representation.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the axis of rotation.
    angle : float
        Angle in radian.
    """
    angle = 2 * math.acos(q[0])
    sin = math.sqrt(1.0 - q[0] ** 2)
    if angle > 0.0:
        if angle < math.pi:
            axis = q[1:4] / sin
        else:
            rotmat = quat_rotmat(q)
            axis, angle = aa_from_rotmat(rotmat)
    else:
        axis = np.array([1.0, 0.0, 0.0])
    return aa_fix(axis, angle, normalize=True)


def quat_to_dcm(q):
    """
    Converts a unit quaternion to a direction cosine matrix.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    return quat_shiftmat(q, forward=True)


def quat_to_euler(q, seq="XYZ", world=True):
    """
    Converts a unit quaternion to an Euler angle sequence.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = quat_rotmat(q)
    return euler_factor_rotmat(rotmat, seq=seq, world=world)


def quat_from_any(orientation):
    """Converts an orientation to unit quaternion from either of the
    following: (1) Quaternion, (2) Euler angles, (3) Axis-angle, or (4)
    Direction cosine matrix.

    Parameters
    ----------
    orientation : dict
        The keys and values are:

        - `'repr'`  = ``'quat'`` | ``'euler'`` | ``'axis_angle'`` | ``'dcm'``
        - `'quat'`  = (4,) *ndarray*
        - `'euler'` = (3,) *ndarray*
        - `'seq'` = ``'XYZ'`` | ``'XZY'`` | ``'YXZ'`` | ``'YZX'``
          | ``'ZXY'`` | ``'ZYX'``
        - `'world'` = ``True`` | ``False``
        - `'axis'` = (3,) *ndarray*
        - `'angle'` = *float*
        - `'dcm'` = (3,3) *ndarray*

        For any value of `'repr'`, only the relevant keys are accessed, the
        rest are ignored. E.g., if `'repr'` = ``'quat'``, only the `'quat'` key
        is necessary, but for `'repr'` = ``'euler'``, the required keys are
        `'euler'`, `'seq'`, and `'world'`.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    ori_repr = orientation["repr"]
    if ori_repr == "quat":
        quat = np.array(orientation["quat"])
    elif ori_repr == "euler":
        euler = np.array(orientation["euler"])
        seq = orientation["seq"]
        world = orientation["world"]
        quat = euler_to_quat(euler, seq=seq, world=world)
    elif ori_repr == "axis_angle":
        axis = np.array(orientation["axis"])
        angle = orientation["angle"]
        quat = aa_to_quat(axis, angle)
    elif ori_repr == "dcm":
        quat = dcm_to_quat(orientation["dcm"])
    else:
        raise ValueError(f"Unrecognized orientation repr {ori_repr}")
    return quat


# AXIS-ANGLE------------------------------------------------------------
def aa_fix(axis, angle, normalize=True):
    """
    Returns a copy of `axis` and `angle` by modifying their values
    to ensure a right handed rotation with `angle` in [0, *pi*).

    Parameters
    ----------
    axis : (3,) ndarray
        Axis of rotation. If this is not a unit vector, set `normalize` to
        ``True``.
    angle : float
        Angle in radian.
    normalize : bool
        Whether to normalize the axis to a unit vector.

    Returns
    -------
    axis : (3,) ndarray
        Modified axis of rotation, possibly normalized.
    angle : float
        Modified angle in radian.

    """
    if normalize:
        norm = np.linalg.norm(axis)
        if not math.isclose(norm, 1.0, abs_tol=1e-14, rel_tol=1e-14):
            axis /= norm
    angle = math.fmod(angle, 2 * math.pi)
    if angle < 0.0:
        angle = -angle
        axis = -axis
    if angle > math.pi:
        angle = 2 * math.pi - angle
        axis = -axis
    return (axis, angle)


def aa_rand(rng):
    """
    Generates a random orientation in *axis-angle* representation.

    The axis is a random vector drawn from a uniform distribution on the
    surface of a unit sphere. The current implementation in based on the
    algorithm from Allen & Tildesley p. 349.

    Parameters
    ----------
    rng : :py:class:`numpy.random.Generator`
        A random number generator.

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian.
    """
    axis = np.zeros((3,))
    # Generate angle: A uniform random number from [0.0, 2*pi)
    angle = 2.0 * math.pi * rng.random()
    while True:
        # Generate two uniform random numbers from [-1, 1)
        zeta1 = 2.0 * rng.random() - 1.0
        zeta2 = 2.0 * rng.random() - 1.0
        zetasq = zeta1**2 + zeta2**2
        if zetasq <= 1.0:
            break
    rt = np.sqrt(1.0 - zetasq)
    axis[0] = 2.0 * zeta1 * rt
    axis[1] = 2.0 * zeta2 * rt
    axis[2] = 1.0 - 2.0 * zetasq
    return aa_fix(axis, angle)


def aa_rotmat(axis, angle):
    """Returns the rotation matrix corresponding to an *axis-angle*
    representation.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. ``0 <= angle < pi``.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.

    """
    R = np.zeros((3, 3))
    sin = np.sin(angle)
    cos = np.cos(angle)
    icos = 1.0 - cos
    R[0, 0] = axis[0] * axis[0] * icos + cos
    R[0, 1] = axis[0] * axis[1] * icos - axis[2] * sin
    R[0, 2] = axis[0] * axis[2] * icos + axis[1] * sin
    R[1, 0] = axis[0] * axis[1] * icos + axis[2] * sin
    R[1, 1] = axis[1] * axis[1] * icos + cos
    R[1, 2] = axis[1] * axis[2] * icos - axis[0] * sin
    R[2, 0] = axis[2] * axis[0] * icos - axis[1] * sin
    R[2, 1] = axis[1] * axis[2] * icos + axis[0] * sin
    R[2, 2] = axis[2] * axis[2] * icos + cos
    return R


def aa_from_rotmat(rotmat):
    """
    Extracts axis and angle from a rotation matrix.

    Parameters
    ----------
    rotmat : (3,3) ndarray
        Rotation matrix (must be orthonormal).

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.

    """
    trace = np.trace(rotmat)
    angle = math.acos((trace - 1) / 2)
    if angle > 0:
        if angle < math.pi:
            u0 = rotmat[2, 1] - rotmat[1, 2]
            u1 = rotmat[0, 2] - rotmat[2, 0]
            u2 = rotmat[1, 0] - rotmat[0, 1]
        else:
            # Find the largest entry in the diagonal of rotmat
            k = np.argmax(np.diag(rotmat))
            if k == 0:
                u0 = (
                    math.sqrt(rotmat[0, 0] - rotmat[1, 1] - rotmat[2, 2] + 1)
                    / 2
                )
                s = 1.0 / (2 * u0)
                u1 = s * rotmat[0, 1]
                u2 = s * rotmat[0, 2]
            elif k == 1:
                u1 = (
                    math.sqrt(rotmat[1, 1] - rotmat[0, 0] - rotmat[2, 2] + 1)
                    / 2
                )
                s = 1.0 / (2 * u1)
                u0 = s * rotmat[0, 1]
                u2 = s * rotmat[1, 2]
            elif k == 2:
                u2 = (
                    math.sqrt(rotmat[2, 2] - rotmat[0, 0] - rotmat[1, 1] + 1)
                    / 2
                )
                s = 1.0 / (2 * u2)
                u0 = s * rotmat[0, 2]
                u1 = s * rotmat[1, 2]
    else:
        u0 = 1.0
        u1 = 0.0
        u2 = 0.0
    return aa_fix(np.array([u0, u1, u2]), angle, normalize=True)


def aa_rotate_vectors(v, axis, angle):
    """
    Rotates vectors about `axis` by `angle`.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. ``0 <= angle < pi``.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.

    """
    rotmat = aa_rotmat(axis, angle)
    return np.dot(v, rotmat.T)


def aa_shiftmat(axis, angle, forward=False):
    """Returns the shifter matrix corresponding to an *axis-angle*
    representation.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. ``0 <= angle < pi``.
    forward: bool
        Whether to shift forward, i.e., along the orientation or shift reverse.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.

    """
    shiftmat = aa_rotmat(-axis, angle)
    if not forward:
        shiftmat = shiftmat.T
    return shiftmat


def aa_shift_vectors(v, axis, angle, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by an axis-angle representation, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. ``0 <= angle < pi``.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = aa_shiftmat(axis, angle, forward=forward)
    return np.dot(v, shiftmat.T)


def aa_shift_tensor2(a, axis, angle, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by an axis-angle representation, shifts second order
    tensors from A to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. ``0 <= angle < pi``.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second order tensor.
    """
    shiftmat = aa_shiftmat(axis, angle, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def aa_shift_tensor3(a, axis, angle, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by an axis-angle representation, shifts third order
    tensors from A to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. ``0 <= angle < pi``.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third order tensor.
    """
    shiftmat = aa_shiftmat(axis, angle, forward=forward)
    return np.einsum("ip,jq,kr,pqr", shiftmat, shiftmat, shiftmat, a)


def aa_to_quat(axis, angle):
    """
    Converts an *axis-angle* representation to a unit quaternion.

    Parameters
    ----------
    axis : (3,) ndarray
        Axis of rotation. This must be a unit vector.
    angle : float
        Angle in radian. ``0 <= angle < pi``.

    Returns
    -------
    q : (4,) ndarray
        Unit quaternion.

    """
    w = math.cos(angle / 2)
    v = math.sin(angle / 2) * axis
    q = np.array([w, v[0], v[1], v[2]])
    return quat_normalized(q)


def aa_to_dcm(axis, angle):
    """Converts an *axis-angle* representation to a direction cosine matrix.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle in radian. ``0 <= angle < pi``.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    dcm = aa_shiftmat(axis, angle, forward=True)
    return dcm


def aa_to_euler(axis, angle, seq="XYZ", world=True):
    """Coverts an *axis-angle* representation to *Euler angles*.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle in radian. ``0 <= angle < pi``.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = aa_rotmat(axis, angle)
    euler = euler_factor_rotmat(rotmat, seq=seq, world=world)
    return euler


def aa_from_any(orientation):
    """Converts an orientation to *axis-angle* from any of the following: (1)
    Quaternion, (2) Euler angles, (3) Axis-angle, or (4) Direction cosine
    matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.any_to_quat`.

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian.

    """
    ori_repr = orientation["repr"]
    if ori_repr == "quat":
        quat = np.array(orientation["quat"])
        axis, angle = quat_to_aa(quat)
    elif ori_repr == "euler":
        euler = np.array(orientation["euler"])
        seq = orientation["seq"]
        world = orientation["world"]
        axis, angle = euler_to_aa(euler, seq=seq, world=world)
    elif ori_repr == "axis_angle":
        axis = np.array(orientation["axis"])
        angle = orientation["angle"]
    elif ori_repr == "dcm":
        axis, angle = dcm_to_aa(orientation["dcm"])
    else:
        raise ValueError(f"Unrecognized orientation repr {ori_repr}")
    return axis, angle


# DIRECTION COSINE MATRIX-----------------------------------------------
def mat_is_dcm(mat):
    """
    Checks if `mat` is a direction cosine matrix or not.

    Parameters
    ----------
    mat : (3,3) ndarray
        Array to check.

    Returns
    -------
    bool
        ``True`` if `mat` is a direction cosine matrix, ``False`` otherwise.
    """
    return mat_is_rotmat(mat)


def dcm_from_axes(A, B):
    """
    Returns the direction cosine matrix of axes(i.e. frame) B with respect to
    axes(i.e. frame) A.

    Parameters
    ----------
    A : (3,3) ndarray
        The rows of A represent the orthonormal basis vectors of frame A.

    B : (3,3) ndarray
        The rows of B represent the orthonormal basis vectors of frame B.

    Returns
    -------
    (3,3) ndarray
        The dcm of frame B w.r.t. frame A.

    """
    return np.dot(B, A.T)


def dcm_rotmat(dcm):
    """Returns the rotation matrix corresponding to a direction cosine matrix.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.
    """
    return dcm.T


def dcm_shiftmat(dcm, forward=False):
    """Returns the shifter matrix corresponding to a direction cosine matrix.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix.
    forward : bool
       Whether to shift forward, i.e., along the orientation or shift reverse.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """
    shiftmat = dcm
    if not forward:
        shiftmat = shiftmat.T
    return shiftmat


def dcm_rotate_vectors(v, dcm):
    """
    Rotates vectors by a direction cosine matrix.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    dcm : (3,3) ndarray
        Direction cosine matrix.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.
    """
    rotmat = dcm_rotmat(dcm)
    return np.dot(v, rotmat.T)


def dcm_shift_vectors(v, dcm, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a direction cosine matrix, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    dcm : (3,3) ndarray
        Direction cosine matrix.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = dcm_shiftmat(dcm, forward=forward)
    return np.dot(v, shiftmat.T)


def dcm_shift_tensor2(a, dcm, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a direction cosine matrix, shifts second-order
    tesnors from A to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor.
    dcm : (3,3) ndarray
        Direction cosine matrix.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second-order tensor.
    """
    shiftmat = dcm_shiftmat(dcm, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def dcm_shift_tensor3(a, dcm, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a direction cosine matrix, shifts third-order
    tesnors from A to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor.
    dcm : (3,3) ndarray
        Direction cosine matrix.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third-order tensor.
    """
    shiftmat = dcm_shiftmat(dcm, forward=forward)
    return np.einsum("ip,jq,kr,pqr", shiftmat, shiftmat, shiftmat, a)


def dcm_to_quat(dcm):
    """
    Converts a direction cosine matrix to a unit quaternion.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix

    Returns
    -------
    q : (4,) ndarray
        Unit quaternion

    """
    mat = dcm_rotmat(dcm)
    axis, angle = aa_from_rotmat(mat)
    return aa_to_quat(axis, angle)


def dcm_to_aa(dcm):
    """Converts a direction cosine matrix to an *axis-angle* representation.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian.
    """
    mat = dcm_rotmat(dcm)
    axis, angle = aa_from_rotmat(mat)
    return (axis, angle)


def dcm_to_euler(dcm, seq="XYZ", world=True):
    """
    Converts a direction cosine matrix to a Euler angles.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        The three Euler angles ``phi`` (rotation about X), ``theta`` (rotation
        about Y), and ``psi`` (rotation about Z).
    """
    mat = dcm_rotmat(dcm)
    euler = euler_factor_rotmat(mat, seq=seq, world=world)
    return euler


def dcm_from_any(orientation):
    """Converts an orientation to direction cosine matrix from any of the
    following: (1) Quaternion, (2) Euler angles, (3) Axis-angle, or (4)
    Direction cosine matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.any_to_quat`.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    ori_repr = orientation["repr"]
    if ori_repr == "quat":
        quat = np.array(orientation["quat"])
        dcm = quat_to_dcm(quat)
    elif ori_repr == "euler":
        euler = np.array(orientation["euler"])
        seq = orientation["seq"]
        world = orientation["world"]
        dcm = euler_to_dcm(euler, seq=seq, world=world)
    elif ori_repr == "axis_angle":
        axis = np.array(orientation["axis"])
        angle = orientation["angle"]
        dcm = aa_to_dcm(axis, angle)
    elif ori_repr == "dcm":
        dcm = dcm_to_quat(orientation["dcm"])
    else:
        raise ValueError(f"Unrecognized orientation repr {ori_repr}")
    return dcm


# EULER ANGLES-----------------------------------------------------------


def euler_rotmat(euler, seq="XYZ", world=True):
    """
    Returns the rotation matrix for a set of Euler angles.

    Parameters
    ----------
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    """
    return _rotmat_euler(euler, seq=seq, world=world)


def euler_factor_rotmat(rotmat, seq="XYZ", world=True):
    """Factorize a rotation matrix to obtain the three Euler angles.

    Parameters
    ----------
    rotmat : (3,3) ndarray
        Rotation matrix
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    """
    return _factor_rotmat(rotmat, seq=seq, world=world)


def euler_shiftmat(euler, seq="XYZ", world=True, forward=False):
    """Returns the shifter matrix for a set of Euler angles.

    Parameters
    ----------
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        Whether to shift forward, i.e., along the orientation or shift reverse.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """

    rotmat = euler_rotmat(euler, seq=seq, world=world)
    if forward:
        shiftmat = rotmat.T
    else:
        shiftmat = rotmat
    return shiftmat


def euler_rotate_vectors(v, euler, seq="XYZ", world=True):
    """
    Rotates vectors with a set of Euler angles.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    euler : (3,) ndarray
        Euler angles.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
        Whether the euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.
    """
    rotmat = euler_rotmat(euler, seq=seq, world=world)
    return np.dot(v, rotmat.T)


def euler_shift_vectors(v, euler, seq="XYZ", world=True, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by three Euler angles, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = euler_shiftmat(euler, seq=seq, world=world, forward=forward)
    return np.dot(v, shiftmat.T)


def euler_shift_tensor2(a, euler, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by three Euler angles, shifts second-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second order tensor.
    """
    shiftmat = euler_shiftmat(euler, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def euler_shift_tensor3(a, euler, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by three Euler angles, shifts third-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third-order tensor.
    """
    shiftmat = euler_shiftmat(euler, forward=forward)
    return np.einsum("ip,jq,kr,pqr", shiftmat, shiftmat, shiftmat, a)


def euler_to_quat(euler, seq="XYZ", world=True):
    """Convert Euler angles to a unit quaternion.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    axis, angle = euler_to_aa(euler, seq=seq, world=world)
    return aa_to_quat(axis, angle)


def euler_to_aa(euler, seq="XYZ", world=True):
    """Convert Euler angles to an *axis-angle* representation.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the direction of the axis.
    angle : float
        Angle in radian.
    """
    rotmat = euler_rotmat(euler, seq=seq, world=world)
    axis, angle = aa_from_rotmat(rotmat)
    return (axis, angle)


def euler_to_dcm(euler, seq="XYZ", world=True):
    """Convert Euler angles to a direction cosine matrix.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix
    """
    dcm = euler_shiftmat(euler, seq=seq, world=world, forward=True)
    return dcm


def euler_to_euler(euler, seq, world, to_seq, to_world):
    """Convert one set of Euler angles to another.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.
    to_seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Convert `euler` with sequence `seq` to the sequence `to_seq`.
    to_world : bool
        Whether the converted euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = euler_rotmat(euler, seq=seq, world=world)
    return euler_factor_rotmat(rotmat, seq=to_seq, world=to_world)


def euler_from_any(orientation, to_seq, to_world):
    """Converts to Euler angles from any of the following: (1)
    Quaternion, (2) Euler angles, (3) Axis-angle, or (4) Direction cosine
    matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.any_to_quat`.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    ori_repr = orientation["repr"]
    if ori_repr == "quat":
        quat = np.array(orientation["quat"])
        euler = quat_to_euler(quat, seq=to_seq, world=to_world)
    elif ori_repr == "euler":
        euler = np.array(orientation["euler"])
        seq = orientation["seq"]
        world = orientation["world"]
        euler = euler_to_euler(euler, seq, world, to_seq, to_world)
    elif ori_repr == "axis_angle":
        axis = np.array(orientation["axis"])
        angle = orientation["angle"]
        euler = aa_to_euler(axis, angle, seq=to_seq, world=to_world)
    elif ori_repr == "dcm":
        euler = dcm_to_euler(orientation["dcm"], seq=to_seq, world=to_world)
    else:
        raise ValueError(f"Unrecognized orientation repr {ori_repr}")
    return euler


def _rotmat_euler(euler, seq="XYZ", world=True):
    """
    Returns the rotation matrix for an Euler angle sequence.

    Parameters
    ----------
    euler : (3,) ndarray
        Euler angles
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    """
    rotmat_funcs = {
        "XYZ": _rotmat_XYZ,
        "XZY": _rotmat_XZY,
        "YXZ": _rotmat_YXZ,
        "YZX": _rotmat_YZX,
        "ZXY": _rotmat_ZXY,
        "ZYX": _rotmat_ZYX,
    }
    if not world:
        euler = -euler
    phi, theta, psi = tuple(euler)
    rotmat = rotmat_funcs[seq](phi, theta, psi)
    if not world:
        rotmat = rotmat.T
    return rotmat


def _factor_rotmat(rotmat, seq="XYZ", world=True):
    factor_rotmat_funcs = {
        "XYZ": _factor_rotmat_XYZ,
        "XZY": _factor_rotmat_XZY,
        "YXZ": _factor_rotmat_YXZ,
        "YZX": _factor_rotmat_YZX,
        "ZXY": _factor_rotmat_ZXY,
        "ZYX": _factor_rotmat_ZYX,
    }
    if not world:
        rotmat = rotmat.T
    factors = factor_rotmat_funcs[seq](rotmat)
    if not world:
        factors = -factors
    return factors


def _rotmat_XYZ(phi, theta, psi):
    rotmat = np.zeros((3, 3))
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0, 0] = cos_theta * cos_psi
    rotmat[0, 1] = sin_phi * sin_theta * cos_psi - cos_phi * sin_psi
    rotmat[0, 2] = cos_phi * sin_theta * cos_psi + sin_phi * sin_psi
    rotmat[1, 0] = cos_theta * sin_psi
    rotmat[1, 1] = sin_psi * sin_theta * sin_phi + cos_phi * cos_psi
    rotmat[1, 2] = cos_phi * sin_theta * sin_psi - sin_phi * cos_psi
    rotmat[2, 0] = -sin_theta
    rotmat[2, 1] = sin_phi * cos_theta
    rotmat[2, 2] = cos_phi * cos_theta
    return rotmat


def _factor_rotmat_XYZ(rotmat):
    if rotmat[2, 0] < 1.0:
        if rotmat[2, 0] > -1.0:
            theta = math.asin(-rotmat[2, 0])
            psi = math.atan2(rotmat[1, 0], rotmat[0, 0])
            phi = math.atan2(rotmat[2, 1], rotmat[2, 2])
        else:
            # Not unique: phi - psi = atan2(-rotmat[1,2], rotmat[1,1])
            theta = math.pi / 2
            psi = -math.atan2(-rotmat[1, 2], rotmat[1, 1])
            phi = 0.0
    else:
        # Not unique: phi + psi = atan2(-rotmat[1,2], rotmat[1,1])
        phi = 0.0
        theta = -math.pi / 2
        psi = math.atan2(-rotmat[1, 2], rotmat[1, 1])
    return np.array([phi, theta, psi])


def _rotmat_XZY(phi, theta, psi):
    rotmat = np.zeros((3, 3))
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0, 0] = cos_theta * cos_psi
    rotmat[0, 1] = sin_phi * sin_theta - cos_phi * cos_theta * sin_psi
    rotmat[0, 2] = cos_phi * sin_theta + sin_phi * cos_theta * sin_psi
    rotmat[1, 0] = sin_psi
    rotmat[1, 1] = cos_phi * cos_psi
    rotmat[1, 2] = -sin_phi * cos_psi
    rotmat[2, 0] = -sin_theta * cos_psi
    rotmat[2, 1] = sin_phi * cos_theta + cos_phi * sin_theta * sin_psi
    rotmat[2, 2] = cos_phi * cos_theta - sin_phi * sin_theta * sin_psi
    return rotmat


def _factor_rotmat_XZY(rotmat):
    if rotmat[1, 0] < 1.0:
        if rotmat[1, 0] > -1.0:
            phi = math.atan2(-rotmat[1, 2], rotmat[1, 1])
            theta = math.atan2(-rotmat[2, 0], rotmat[0, 0])
            psi = math.asin(rotmat[1, 0])
        else:
            # Not unique: phi - theta = atan2(rotmat[2,1], rotmat[2,2])
            phi = 0.0
            theta = -math.atan2(rotmat[2, 1], rotmat[2, 2])
            psi = -math.pi / 2
    else:
        # Not unique: phi + theta = atan2(rotmat[2,1], rotmat[2,2])
        phi = 0.0
        theta = math.atan2(rotmat[2, 1], rotmat[2, 1])
        psi = math.pi / 2
    return np.array([phi, theta, psi])


def _rotmat_YXZ(phi, theta, psi):
    rotmat = np.zeros((3, 3))
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0, 0] = cos_theta * cos_psi - sin_phi * sin_theta * sin_psi
    rotmat[0, 1] = -cos_phi * sin_psi
    rotmat[0, 2] = sin_theta * cos_psi + sin_phi * cos_theta * sin_psi
    rotmat[1, 0] = sin_phi * sin_theta * cos_psi + cos_theta * sin_psi
    rotmat[1, 1] = cos_phi * cos_psi
    rotmat[1, 2] = sin_theta * sin_psi - sin_phi * cos_theta * cos_psi
    rotmat[2, 0] = -cos_phi * sin_theta
    rotmat[2, 1] = sin_phi
    rotmat[2, 2] = cos_phi * cos_theta
    return rotmat


def _factor_rotmat_YXZ(rotmat):
    if rotmat[2, 1] < 1.0:
        if rotmat[2, 1] > -1.0:
            phi = math.asin(rotmat[2, 1])
            theta = math.atan2(-rotmat[2, 0], rotmat[2, 2])
            psi = math.atan2(-rotmat[0, 1], rotmat[1, 1])
        else:
            # Not unique: theta - psi = atan2(rotmat[0,2], rotmat[0,0])
            phi = -math.pi / 2
            theta = 0.0
            psi = -math.atan2(rotmat[0, 2], rotmat[0, 0])
    else:
        # Not unique: theta + psi = atan2(rotmat[0,2], rotmat[0,0])
        phi = math.pi / 2
        theta = 0.0
        psi = math.atan2(rotmat[0, 2], rotmat[0, 0])
    return np.array([phi, theta, psi])


def _rotmat_YZX(phi, theta, psi):
    rotmat = np.zeros((3, 3))
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0, 0] = cos_theta * cos_psi
    rotmat[0, 1] = -sin_psi
    rotmat[0, 2] = sin_theta * cos_psi
    rotmat[1, 0] = sin_phi * sin_theta + cos_phi * cos_theta * sin_psi
    rotmat[1, 1] = cos_phi * cos_psi
    rotmat[1, 2] = cos_phi * sin_theta * sin_psi - sin_phi * cos_theta
    rotmat[2, 0] = sin_phi * cos_theta * sin_psi - cos_phi * sin_theta
    rotmat[2, 1] = sin_phi * cos_psi
    rotmat[2, 2] = sin_phi * sin_theta * sin_psi + cos_phi * cos_theta
    return rotmat


def _factor_rotmat_YZX(rotmat):
    if rotmat[0, 1] < 1.0:
        if rotmat[0, 1] > -1.0:
            phi = math.atan2(rotmat[2, 1], rotmat[1, 1])
            theta = math.atan2(rotmat[0, 2], rotmat[0, 0])
            psi = math.asin(-rotmat[0, 1])
        else:
            # Not unique: theta - phi = atan2(-rotmat[2,0], rotmat[2,2])
            phi = -math.atan2(-rotmat[2, 0], rotmat[2, 2])
            theta = 0.0
            psi = math.pi / 2
    else:
        # Not unique: theta + phi = atan2(-rotmat[2,0], rotmat[2,2])
        phi = math.atan2(-rotmat[2, 0], rotmat[2, 2])
        theta = 0.0
        psi = -math.pi / 2
    return np.array([phi, theta, psi])


def _rotmat_ZXY(phi, theta, psi):
    rotmat = np.zeros((3, 3))
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0, 0] = cos_theta * cos_psi + sin_phi * sin_theta * sin_psi
    rotmat[0, 1] = sin_phi * sin_theta * cos_psi - cos_theta * sin_psi
    rotmat[0, 2] = cos_phi * sin_theta
    rotmat[1, 0] = cos_phi * sin_psi
    rotmat[1, 1] = cos_phi * cos_psi
    rotmat[1, 2] = -sin_phi
    rotmat[2, 0] = sin_phi * cos_theta * sin_psi - sin_theta * cos_psi
    rotmat[2, 1] = sin_phi * cos_theta * cos_psi + sin_theta * sin_psi
    rotmat[2, 2] = cos_phi * cos_theta
    return rotmat


def _factor_rotmat_ZXY(rotmat):
    if rotmat[1, 2] < 1.0:
        if rotmat[1, 2] > -1.0:
            phi = math.asin(-rotmat[1, 2])
            theta = math.atan2(rotmat[0, 2], rotmat[2, 2])
            psi = math.atan2(rotmat[1, 0], rotmat[1, 1])
        else:
            # Not unique: psi - theta = atan2(-rotmat[0,1], rotmat[0,0])
            phi = math.pi / 2
            theta = -math.atan2(-rotmat[0, 1], rotmat[0, 0])
            psi = 0.0
    else:
        # Not unique: psi + theta = atan2(-rotmat[0,1], rotmat[0,0])
        phi = -math.pi / 2
        theta = math.atan2(-rotmat[0, 1], rotmat[0, 0])
        psi = 0.0
    return np.array([phi, theta, psi])


def _rotmat_ZYX(phi, theta, psi):
    rotmat = np.zeros((3, 3))
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0, 0] = cos_theta * cos_psi
    rotmat[0, 1] = -cos_theta * sin_psi
    rotmat[0, 2] = sin_theta
    rotmat[1, 0] = sin_phi * sin_theta * cos_psi + cos_phi * sin_psi
    rotmat[1, 1] = cos_phi * cos_psi - sin_phi * sin_theta * sin_psi
    rotmat[1, 2] = -sin_phi * cos_theta
    rotmat[2, 0] = sin_phi * sin_psi - cos_phi * sin_theta * cos_psi
    rotmat[2, 1] = sin_phi * cos_psi + cos_phi * sin_theta * sin_psi
    rotmat[2, 2] = cos_phi * cos_theta
    return rotmat


def _factor_rotmat_ZYX(rotmat):
    if rotmat[0, 2] < 1.0:
        if rotmat[0, 2] > -1.0:
            phi = math.atan2(-rotmat[1, 2], rotmat[2, 2])
            theta = math.asin(rotmat[0, 2])
            psi = math.atan2(-rotmat[0, 1], rotmat[0, 0])
        else:
            # Not unique: psi - phi = atan2(rotmat[1,0], rotmat[1,1])
            phi = -math.atan2(rotmat[1, 0], rotmat[1, 1])
            theta = -math.pi / 2
            psi = 0.0
    else:
        # Not unique: psi + phi = atan2(rotmat[1,0], rotmat[1,1])
        phi = math.atan2(rotmat[1, 0], rotmat[1, 1])
        theta = math.pi / 2
        psi = 0.0
    return np.array([phi, theta, psi])
