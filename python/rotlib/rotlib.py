import math
import numpy as np

def _normalize(v):
    """
    Normalizes a vector in-place such that its magnitude is unity.

    Parameters
    ----------
    v : (n,) ndarray
        On entry, a finite non-zero vector. On exit, the normalized vector.

    Returns
    -------
    None
    """
    scale = np.max(np.abs(v))
    if scale == 0.0:
        raise ValueError("Cannot normalize a zero vector.")
    v /= scale
    v /= np.linalg.norm(v)

# *****************************************************************************#
#                               UTILITY FUNCTIONS                             #
# *****************************************************************************#


def align(v, old, new):
    """
    Rotate vectors such that a set of at most three mutually orthogonal unit
    vectors aligns with another set of mutually orthogonal unit vectors. The
    two sets of mutually orthogonal unit vectors must contain the same number
    of unit vectors.

    `v`, `old`, and `new` store their constituent vectors row-wise, with
    all vectors resolved in the same fixed underlying basis.

    The vectors in `v` are rotated using a rotation matrix R such that
    ``new[i,:] = old[i,:] @ R.T``, for each unit vector `i` (stored in the `i` th
    row), with the single vector (`n=1`) case treated equivalently.

    If `old` is a single vector, a unique rotation cannot be defined.  In this
    case, the shortest rotation mapping `old` to `new` is chosen. For
    anti-parallel directions, this is a half-turn about a chosen axis
    perpendicular to `old`.

    Parameters
    ----------
    v : (3,) | (m,3) ndarray
        Vectors to align.
    old : (3,) | (1,3) | (3,1) | (n,3) ndarray
        A set of at most three mutually orthogonal unit vectors. A single unit
        vector can be of shape (3,), (1,3), or (3,1). Else `n` must be 2 or 3.
        If `n = 3`, the three unit vectors must form a right-handed basis.
    new : (3,) | (1,3) | (3,1) | (n,3) ndarray
        Another set of mutually orthogonal unit vectors of the same shape as
        `old`. If `n = 3`, the three unit vectors must form a right-handed
        basis.

    Returns
    -------
    out : ndarray
        Aligned vectors. `out` is of the same shape as `v`. The components in
        `out` are resolved in the same underlying basis as the vectors in `v`.
    """
    if np.any(old.shape != new.shape):
        raise ValueError(
            f"`old` of shape {old.shape} differs from`new` of shape {new.shape}."
        )
    n = old.size // 3
    if n == 1:
        old_ = old.ravel()
        new_ = new.ravel()
        # Angle between old and new
        cos_angle = np.dot(old_, new_)
        axis = np.cross(old_, new_)
        sin_angle = math.hypot(*axis) #always non-negative
        #angle will be in [0, pi] as sin_angle >= 0.
        angle = math.atan2(sin_angle, cos_angle)
        #try normalizing the axis
        try:
            _normalize(axis)
            estat = 0
        except ValueError:
            estat = 1

        if estat == 0:
            #Non-collinear vectors: standard case
            out = aa_rotate_vectors(v, axis, angle)
        else:
            #Collinear vectors: parallel or anti-parallel
            if cos_angle > 0.0:
                #angle = 0; no rotation necessary
                out = v.copy()
            else:
                #angle = pi: Rotation axis not unique. Pick a vector by
                #crossing old_ with a cartesian direction least parallel to it.
                #1. Find cartesian direction least parallel to old_.
                min_idx = np.argmin(abs(old_))
                #2. Create that basis vector
                vtmp = np.zeros((3,), dtype=np.float64)
                vtmp[min_idx] = 1.0
                #3. Cross old_ with vtmp and normalize to get the rotation
                # axis
                axis = np.cross(old_, vtmp) #axis will always have +ve norm.
                axis, angle = aa_fix(axis, angle, True)
                #Rotate
                out = aa_rotate_vectors(v, axis, angle)
        return out
    elif n == 2:
        # The third vector orthogonal to `old`
        z_old = np.cross(old[0,:], old[1,:])
        _normalize(z_old)
        # The third vector orthogonal to `new`
        z_new = np.cross(new[0,:], new[1,:])
        _normalize(z_new)
        # The `old` and new `frames`
        axes_old = np.vstack((old, z_old))
        axes_new = np.vstack((new, z_new))
        mat = np.dot(axes_new.T, axes_old)
        return np.dot(v, mat.T)
    elif n == 3:
        mat = np.dot(new.T, old)
        return np.dot(v, mat.T)


def mat_is_rotmat(mat, tol=1e-14):
    """
    Checks if `mat` is a rotation matrix or not.

    Parameters
    ----------
    mat : (3,3) ndarray
        Array to check.
    tol : float
        Absolute tolerance for comparing determinant and orthogonality.

    Returns
    -------
    bool
        ``True`` if `mat` is a rotation matrix, ``False`` otherwise.
    """
    det_mat = np.linalg.det(mat)
    det_is_one = abs(det_mat-1.0) <= tol
    is_orthogonal = np.abs(mat @ mat.T - np.identity(3)).max() <= tol
    return is_orthogonal and det_is_one


# QUATERNION-----------------------------------------------------------

def quat_canonicalized(q):
    """
    Modifies a non-canonical unit quaternion in-place to bring it to the
    canonical form.

    No modifications are performed if the input is already in canonical form.

    Parameters
    ----------
    q : (4,) ndarray
        On entry, a unit quaternion. If `q` is not in canonical form, on return, 
        `q` is modified to bring it to the canonical form, else `q` is left as
        is.

    Returns
    -------
    None

    """
    inz = np.nonzero(q)[0][0]
    if q[inz] < 0.0:
        q[:] = -q 


def quat_rand(rng):
    """
    Returns a unit quaternion from a uniform distribution on the space of
    rotations in 3D.

    Parameters
    ----------
    rng : :py:class:`numpy.random.Generator`
        A random number generator.

    Returns
    -------
    (4,) ndarray
        A unit quaternion in canonical form.
    """
    u = rng.random(size=(3,), dtype=np.float64)
    sqrt_u0 = math.sqrt(u[0])
    sqrt_1mu0 = math.sqrt(1.0 - u[0])
    cos_u1 = math.cos(2*math.pi*u[1])
    sin_u1 = math.sin(2*math.pi*u[1])
    cos_u2 = math.cos(2*math.pi*u[2])
    sin_u2 = math.sin(2*math.pi*u[2])
    q = np.array([sqrt_1mu0*cos_u1, sqrt_1mu0*sin_u1, sqrt_u0*sin_u2,
                  sqrt_u0*cos_u2], dtype=np.float64)
    quat_canonicalized(q)
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
        On entry, a quaternion. On return, the conjugated quaternion.

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
        On entry, a non-zero quaternion. On return, the inverted quaternion.

    Returns
    -------
    (4,) ndarray
        Inverted quaternion.

    """
    scale = np.max(np.abs(q))
    if scale == 0.0:
        raise ValueError("Cannot invert a zero quaternion.")
    q_scaled = q/scale
    d = np.dot(q_scaled, q_scaled)
    quat_conjugated(q_scaled)
    q[:] = (q_scaled/d)/scale
    return q


def quat_normalized(q, canonicalize=False):
    """
    Normalizes a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        On entry, a non-zero quaternion. On return, the normalized quaternion
        (in canonical form if `canonicalize` is ``True``).
    canonicalize : bool
        Whether to canonicalize the output.

    Returns
    -------
    (4,) ndarray
        Normalized quaternion (in canonical form if `canonicalize` is ``True``).
    """
    _normalize(q)
    if canonicalize:
        quat_canonicalized(q)
    return q


def quat_is_normalized(q, tol=1e-14):
    """
    Checks whether a quaternion is normalized, i.e. whether it is a unit
    quaternion, within a tolerance `tol`.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion
    tol : float
        Absolute allowable deviation of the norm from one.

    Returns
    -------
    bool
        ``True`` if `q` is normalized within tolerance, ``False`` otherwise.
    """
    norm = math.hypot(*q)
    return abs(norm - 1.0) <= tol


def quat_prod(p, q):
    """
    Returns the product of two quaternions following Hamilton's convention. In
    case of unit quaternions describing rotations in a common frame, `q`'s
    rotation is applied first, followed by that of `p`.

    Parameters
    ----------
    p : (4,) ndarray
        Quaternion.
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Product of two quaternions.

    """
    p0, p1, p2, p3 = tuple(p)
    prod_mat = np.array([[p0, -p1, -p2, -p3],
                         [p1,  p0, -p3,  p2],
                         [p2,  p3,  p0, -p1],
                         [p3, -p2,  p1,  p0]])  # fmt: skip
    pq = np.dot(prod_mat, q)
    return pq


def quat_angle_between(p, q):
    """
    Returns the smallest rotation angle between the orientations represented by
    two unit quaternions.

    Parameters
    ----------
    p : (4,) ndarray
        Unit quaternion.
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    float
        Relative rotation angle in radians, in the interval [0, pi]. 
        Opposite-sign quaternions representing the same orientation have zero
        separation.
    """
    q_ = q.copy()
    cos_angle = np.dot(p, q_)
    if cos_angle < 0.0:
        #Negate q_ to ensure both p and q_ lie on the same hemisphere
        q_[:] = -q_[:]
    d_mag = math.hypot(*(p-q_))
    s_mag = math.hypot(*(p+q_))
    angle = 4*math.atan2(d_mag, s_mag)
    return angle


def quat_interpolate(q1, q2, t):
    """
    Performs spherical linear interpolation between two unit quaternions
    following the shortest path.

    The second endpoint may be negated to select the shortest path.

    At a relative orientation (physical) of exactly *pi*, two equally short
    paths exist; the supplied endpoint signs determine which is selected.

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
    q2_ = q2.copy()
    cos_angle = np.dot(q1, q2_)
    if cos_angle < 0.0:
        #Negate q2_ to ensure both q1 and q2_ lie on the same hemisphere
        q2_[:] = -q2_[:]
    d_mag = math.hypot(*(q1 - q2_))
    s_mag = math.hypot(*(q1 + q2_))
    theta = 2*math.atan2(d_mag, s_mag)
    #SLERP
    c = np.sinc(theta/np.pi)
    a = (1 - t) * np.sinc( (1-t)*theta/np.pi ) / c
    b = t * np.sinc(t*theta/np.pi) / c
    q = a*q1 + b*q2_
    return quat_normalized(q, canonicalize=False)


def quat_deriv_to_angvel_mat(q, world=False):
    """
    Returns the matrix mapping the derivative of a unit quaternion to angular
    velocity.

    Let A (world frame) and B (body frame) be reference frames such that an
    active rotation in A by a unit quaternion ``q`` aligns its basis vectors
    with those of frame B.  Then the matrix M returned is such that ``omega_i =
    M_ij * qdot_j``, where ``omega_i`` is the angular velocity of B relative to
    A resolved in A (or B) and ``qdot_j`` is the time derivative of ``q``.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    world : bool
        If ``True``, the angular velocity is resolved in frame A (world frame),
        else, it is resolved in frame B (body frame). Default is ``False``.

    Returns
    -------
    (3,4) ndarray
        Angular velocity matrix.
    """
    q0, q1, q2, q3 = tuple(q)
    if world:
        M = 2*np.array([[-q1,  q0, -q3,  q2],
                        [-q2,  q3,  q0, -q1],
                        [-q3, -q2,  q1,  q0]])  # fmt: skip
    else:
        M = 2*np.array([[-q1,  q0,  q3, -q2],
                        [-q2, -q3,  q0,  q1],
                        [-q3,  q2, -q1,  q0]])  # fmt: skip
    return M


def quat_deriv_to_angvel(q, qdot, world=False):
    """
    Calculates the angular velocity from a unit quaternion and its time
    derivative.

    Let A (world frame) and B (body frame) be reference frames such that an
    active rotation in A by a unit quaternion ``q`` aligns its basis vectors
    with those of frame B. Let ``qdot`` be the time derivative of ``q``. The
    quantity returned is the angular velocity of B relative to A resolved in A
    (or B).

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    qdot : (4,) ndarray
        Time derivative of `q`. In general, this is not a unit vector. Note that
        `qdot` is orthogonal to `q`, i.e. ``q . qdot = 0``.
    world : bool
        If ``True``, the angular velocity is resolved in frame A (world frame),
        else, it is resolved in frame B (body frame). Default is ``False``.

    Returns
    -------
    (3,) ndarray
        Angular velocity resolved in A (world frame) if `world=True`, or in B
        (body frame) otherwise.
    """
    mat = quat_deriv_to_angvel_mat(q, world)
    return np.dot(mat, qdot)


def quat_deriv_from_angvel_mat(q, world=False):
    """
    Returns the matrix mapping angular velocity to time derivative of a unit
    quaternion.

    Let A (world frame) and B (body frame) be reference frames such that an
    active rotation in A by a unit quaternion ``q`` aligns its basis vectors
    with those of frame B. Then the matrix M returned is such that ``qdot_i =
    M_ij * omega_j``, where ``omega_j`` is the angular velocity of B relative
    to A resolved in A (or B) and ``qdot_i`` is the time derivative of ``q``.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    world : bool
        If ``True``, the angular velocity is resolved in frame A (world frame),
        else, it is resolved in frame B (body frame). Default is ``False``.

    Returns
    -------
    (4,3) ndarray
        Quaternion derivative matrix.
    """
    q0, q1, q2, q3 = tuple(q)
    if world:
        M = 0.5*np.array([[-q1, -q2, -q3],
                          [ q0,  q3, -q2],
                          [-q3,  q0,  q1],
                          [ q2, -q1,  q0]])  # fmt: skip
    else:
        M = 0.5*np.array([[-q1, -q2, -q3],
                          [ q0, -q3,  q2],
                          [ q3,  q0, -q1],
                          [-q2,  q1,  q0]])  # fmt: skip
    return M


def quat_deriv_from_angvel(q, ang_vel, world=False):
    """
    Calculates the time derivative of a unit quaternion from angular velocity.

    Let A (world frame) and B (body frame) be reference frames such that an
    active rotation in A by a unit quaternion ``q`` aligns its basis vectors
    with those of frame B. Let ``ang_vel`` be the angular velocity of B
    relative to A resolved in A (or B). The quantity returned is the time
    derivative of ``q``.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion. An active rotation of frame A by ``q`` aligns it with
        frame B.
    ang_vel : (3,) ndarray
        Angular velocity of frame B relative to frame A resolved in A if
        ``world=True``, or in B otherwise.
    world : bool
        If ``True``, the angular velocity is resolved in frame A (world frame),
        else, it is resolved in frame B (body frame). Default is ``False``.

    Returns
    -------
    (4,) ndarray
        Derivative of `q`. It is orthogonal to `q` and not a unit vector in
        general.
    """
    mat = quat_deriv_from_angvel_mat(q, world)
    qdot = np.dot(mat, ang_vel)
    return qdot


def quat_rotmat(q):
    """
    Returns the matrix for an active rotation corresponding to a unit
    quaternion.

    For a column vector v of length 3, the rotated vector ``v' = R v``, where
    ``R`` is the returned matrix.  The components of ``v`` and ``v'`` are
    expressed in the same fixed orthonormal basis. For a row vector, the
    equivalent rotated vector ``v' = v R^T``.

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

    The orientation of frame B relative to frame A is specified by a unit
    quaternion ``q``, such that an active rotation in A of its basis vectors
    aligns them with those of B. The shifter matrix ``S`` transforms the
    components of a fixed vector resolved in A to those resolved in B or
    vice-versa according to ``v_B_i = S_ij * v_A_j`` or ``v_A_i = S_ij *
    v_B_j``.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

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

    The components of `v` and the rotated vectors are all expressed in the same
    fixed basis.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows.
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
    Transform vector components between frames A and B.
    
    The orientation of B relative to A is specified by a unit quaternion `q`.
    This orientation is used for both transformation directions.  The vectors
    remain unchanged; only their components are expressed in a different frame.
    
    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows.
        Expressed in A if `forward` is ``True``, or in B otherwise.
    q : (4,) ndarray
        Unit quaternion expressed in A, regardless of `forward`.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.
    
    Returns
    -------
    (3,) or (n,3) ndarray
        Components of the same vectors expressed in B if `forward` is ``True``,
        or in A otherwise. Has the same shape as `v`.
    """
    shiftmat = quat_shiftmat(q, forward=forward)
    return np.dot(v, shiftmat.T)


def quat_shift_tensor2(a, q, forward=False):
    """
    Transform second-order tensor components between frames A and B.
    
    The orientation of B relative to A is specified by a unit quaternion `q`.
    This orientation is used for both transformation directions.  The tensor
    remains unchanged; only its components are expressed in a different frame.
    
    Parameters
    ----------
    a : (3,3) ndarray
        Components of a second-order tensor.  Expressed in A if `forward` is
        ``True``, or in B otherwise.
    q : (4,) ndarray
        Unit quaternion expressed in A, regardless of `forward`.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.
    
    Returns
    -------
    (3,3) ndarray
        Components of the same second-order tensor expressed in B if `forward`
        is ``True``, or in A otherwise.
    """
    shiftmat = quat_shiftmat(q, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def quat_shift_tensor3(a, q, forward=False):
    """
    Transform third-order tensor components between frames A and B.
    
    The orientation of B relative to A is specified by a unit quaternion `q`.
    This orientation is used for both transformation directions. The tensor
    remains unchanged; only its components are expressed in a different frame.
    
    Parameters
    ----------
    a : (3,3,3) ndarray
        Components of a third-order tensor.  Expressed in A if `forward` is
        ``True``, or in B otherwise.
    q : (4,) ndarray
        Unit quaternion expressed in A, regardless of `forward`.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.
    
    Returns
    -------
    (3,3,3) ndarray
        Components of the same third-order tensor expressed in B if `forward`
        is ``True``, or in A otherwise.
    """
    shiftmat = quat_shiftmat(q, forward=forward)
    return np.einsum("ip,jq,kr,pqr", shiftmat, shiftmat, shiftmat, a)


def quat_to_aa(q):
    """
    Converts a unit quaternion to an *axis-angle* representation.
    Note that an identity quaternion will return axis [1,0,0].

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the axis of rotation.
    angle : float
        Angle in radian (lies in [0, *pi*]).
    """
    q_ = q.copy()
    if q_[0] < 0.0:
        q_[:] = -q_
    s = np.max(np.abs(q_[1:4]))
    if s == 0:
        #Identity rotation (angle = 0)
        axis = np.array([1.0, 0.0, 0.0])
        angle = 0.0
    else:
        #Usual case
        p = q_[1:4]/s
        r = math.hypot(*p)
        axis = p/r
        angle = 2*math.atan2(s*r, q_[0])
    return (axis, angle)


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


def quat_to_euler(q, seq="XYZ", world=True, tol=1e-14):
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
    tol : float
        Tolerance (must be >= 0) for setting one outer angle to zero.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = quat_rotmat(q)
    return euler_factor_rotmat(rotmat, seq=seq, world=world, tol=tol)


def quat_from_any(orientation):
    """Converts an orientation to unit quaternion from either of the
    following: (1) Quaternion, (2) Euler angles, (3) Axis-angle, or (4)
    Direction cosine matrix.

    Parameters
    ----------
    orientation : dict
        The keys and values are:

        * `'repr'`  = ``'quat'`` | ``'euler'`` | ``'axis_angle'`` | ``'dcm'``
        * `'quat'`  = (4,) *ndarray* representing a unit quaternion.
        * `'euler'` = (3,) *ndarray*
        * `'seq'` = ``'XYZ'`` | ``'XZY'`` | ``'YXZ'`` | ``'YZX'``
          | ``'ZXY'`` | ``'ZYX'``
        * `'world'` = ``True`` | ``False``
        * `'axis'` = (3,) *ndarray* representing a unit vector.
        * `'angle'` = *float*
        * `'dcm'` = (3,3) *ndarray* representing a direction cosine matrix.

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
    For a given orientation in *axis-angle* representation, returns a modified
    `axis` and `angle` to ensure that `angle` lies in [0, *pi*] while
    preserving the rotation. 

    Parameters
    ----------
    axis : (3,) ndarray
        A non-zero vector representing the axis of rotation.
    angle : float
        Angle in radian.
    normalize : bool
        Whether to normalize the axis to a unit vector.

    Returns
    -------
    ax : (3,) ndarray
        Modified axis of rotation. If `normalize`=``True``, `ax` is normalized,
        else, it has the same magnitude as `axis` and may need to be normalized
        before using other axis-angle functions.
    ang : float
        Modified angle in radian (lies in [0, *pi*]).

    """
    ax = axis.copy()
    if normalize:
        _normalize(ax)
    #Get principal signed angle in [-pi, pi]
    ang = math.atan2(math.sin(angle), math.cos(angle))
    #Make the angle non-negative
    if ang < 0.0:
        ang = -ang
        ax[:] = -ax
    return (ax, ang)


def aa_rand(rng):
    """
    Generates a random orientation in *axis-angle* representation drawn from a
    uniform distribution on the space of rotations in 3D.

    Parameters
    ----------
    rng : :py:class:`numpy.random.Generator`
        A random number generator.

    Returns
    -------
    axis : (3,) ndarray
        A unit vector along the axis of rotation.
    angle : float
        Angle in radian (lies in [0, *pi*]).
    """
    q = quat_rand(rng)
    axis, angle = quat_to_aa(q) 
    return axis, angle


def aa_rotmat(axis, angle):
    """
    Returns the matrix for an active rotation about `axis` through `angle`,
    with positive angles following the right-hand rule.

    For a column vector v of length 3, the rotated vector v' = R v, where R is
    the returned matrix.  The components of v, v', and `axis` are expressed in
    the same fixed orthonormal basis. For a row vector, the equivalent rotated
    vector v' = v R^T.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. Positive angles follow the right-hand
        rule.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.

    """
    R = np.zeros((3, 3))
    sin = math.sin(angle)
    cos = math.cos(angle)
    sin_2 = math.sin(angle/2)
    icos = 2*sin_2*sin_2 # equals 1.0 - cos
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

    If `rotmat` is an identity matrix, the axis returned is [1,0,0]. If the
    angle is exactly *pi*, either sign of the axis is valid.

    Parameters
    ----------
    rotmat : (3,3) ndarray
        Rotation matrix (must be orthonormal with determinant equal to one).

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian within [0, *pi*].

    """
    axis = np.zeros((3,), dtype=np.float64)
    #Cos of the angle from the trace
    trace = np.trace(rotmat)
    cos_angle = (trace-1)/2
    #Clamp cos_angle to [-1, 1]
    cos_angle = max(-1, min(1, cos_angle))
    #Sin of the angle from the anti-symmetric vector w
    w = 0.5*np.array([ rotmat[2,1]-rotmat[1,2],
                       rotmat[0,2]-rotmat[2,0],
                       rotmat[1,0]-rotmat[0,1] ])
    sin_angle = math.hypot(*w)
    #Calculate the angle
    angle = math.atan2(sin_angle, cos_angle)

    if cos_angle >= 0.0:
        #Angle lies in [0, pi/2]
        if sin_angle > 0.0:
            axis[:] = w/sin_angle
        else:
            #angle = 0, axis non-unique
            angle = 0.0
            axis[:] = np.array([1.0, 0.0, 0.0])
    else:
        # Find the largest entry on the diagonal of rotmat
        k = np.argmax(np.diag(rotmat))
        if k == 0:
            axis[0] = math.sqrt((rotmat[0,0]-cos_angle)/(1.0-cos_angle))
            s = 1.0/(2*(1.0-cos_angle)*axis[0])
            axis[1] = s*(rotmat[1,0]+rotmat[0,1])
            axis[2] = s*(rotmat[2,0]+rotmat[0,2])
        elif k == 1:
            axis[1] = math.sqrt((rotmat[1,1]-cos_angle)/(1.0-cos_angle))
            s = 1.0/(2*(1.0-cos_angle)*axis[1])
            axis[0] = s*(rotmat[0,1]+rotmat[1,0])
            axis[2] = s*(rotmat[2,1]+rotmat[1,2])
        elif k == 2:
            axis[2] = math.sqrt((rotmat[2,2]-cos_angle)/(1.0-cos_angle))
            s = 1.0/(2*(1.0-cos_angle)*axis[2])
            axis[0] = s*(rotmat[0,2]+rotmat[2,0])
            axis[1] = s*(rotmat[1,2]+rotmat[2,1])

    #Resolve axis sign
    if np.dot(axis, w) < 0.0:
        axis[:] = -axis
    return aa_fix(axis, angle, normalize=True)


def aa_rotate_vectors(v, axis, angle):
    """
    Rotates vectors about `axis` by `angle`.

    The components of `axis`, `v`, and the rotated vectors are all expressed in
    the same fixed basis.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. Positive angles follow the right-hand rule.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.

    """
    rotmat = aa_rotmat(axis, angle)
    return np.dot(v, rotmat.T)


def aa_shiftmat(axis, angle, forward=False):
    """
    Returns the shifter matrix corresponding to an *axis-angle* representation.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian. Positive angles follow the right-hand rule.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

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
    Transform vector components between frames A and B.
    
    The orientation of B relative to A is specified by `axis` and `angle`,
    with `axis` expressed in A and positive angles following the right-hand
    rule. This orientation is used for both transformation directions.
    The vectors remain unchanged; only their components are expressed in a
    different frame.
    
    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows.
        Expressed in A if `forward` is ``True``, or in B otherwise.
    axis : (3,) ndarray
        Unit axis vector expressed in A, regardless of `forward`.
    angle : float
        Angle specifying the orientation of B relative to A, in radians.
        Positive angles follow the right-hand rule.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.
    
    Returns
    -------
    (3,) or (n,3) ndarray
        Components of the same vectors expressed in B if `forward` is ``True``,
        or in A otherwise. Has the same shape as `v`.
    """
    shiftmat = aa_shiftmat(axis, angle, forward=forward)
    return np.dot(v, shiftmat.T)


def aa_shift_tensor2(a, axis, angle, forward=False):
    """
    Transforms second-order tensor components between frames A and B.
    
    The orientation of B relative to A is specified by `axis` and `angle`,
    with `axis` expressed in A and positive angles following the right-hand
    rule. This orientation is used for both transformation directions.
    The tensor remains unchanged; only its components are expressed in a
    different frame.
    
    Parameters
    ----------
    a : (3,3) ndarray
        Components of a second-order tensor, expressed in A if `forward` is
        ``True``, or in B otherwise.
    axis : (3,) ndarray
        Unit axis vector expressed in A, regardless of `forward`.
    angle : float
        Angle specifying the orientation of B relative to A, in radians.
        Positive angles follow the right-hand rule.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.
    
    Returns
    -------
    (3,3) ndarray
        Components of the second-order tensor expressed in B if `forward` is
        ``True``, or in A otherwise.

    """
    shiftmat = aa_shiftmat(axis, angle, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def aa_shift_tensor3(a, axis, angle, forward=False):
    """
    Transforms third-order tensor components between frames A and B.
    
    The orientation of B relative to A is specified by `axis` and `angle`,
    with `axis` expressed in A and positive angles following the right-hand
    rule. This orientation is used for both transformation directions.
    The tensor remains unchanged; only its components are expressed in a
    different frame.
    
    Parameters
    ----------
    a : (3,3,3) ndarray
        Components of a third-order tensor, expressed in A if `forward` is
        ``True``, or in B otherwise.
    axis : (3,) ndarray
        Unit axis vector expressed in A, regardless of `forward`.
    angle : float
        Angle specifying the orientation of B relative to A, in radians.
        Positive angles follow the right-hand rule.
    forward : bool, optional
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.
    
    Returns
    -------
    (3,3,3) ndarray
        Components of the third-order tensor expressed in B if `forward` is
        ``True``, or in A otherwise.

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
        Angle in radian.

    Returns
    -------
    q : (4,) ndarray
        Unit quaternion in the canonical form.

    """
    w = math.cos(angle / 2)
    v = math.sin(angle / 2) * axis
    q = np.array([w, v[0], v[1], v[2]])
    quat_normalized(q, canonicalize=True)
    return q


def aa_to_dcm(axis, angle):
    """Converts an *axis-angle* representation to a direction cosine matrix.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle in radian.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    dcm = aa_shiftmat(axis, angle, forward=True)
    return dcm


def aa_to_euler(axis, angle, seq="XYZ", world=True, tol=1e-14):
    """Converts an *axis-angle* representation to *Euler angles*.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    tol : float
        Tolerance (must be >= 0) for setting one outer angle to zero.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = aa_rotmat(axis, angle)
    euler = euler_factor_rotmat(rotmat, seq=seq, world=world, tol=tol)
    return euler


def aa_from_any(orientation):
    """Converts an orientation to *axis-angle* from any of the following: (1)
    Quaternion, (2) Euler angles, (3) Axis-angle, or (4) Direction cosine
    matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.quat_from_any`.

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian (lies in [0, *pi*]).

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
        axis, angle = aa_fix(axis, angle, normalize=True)
    elif ori_repr == "dcm":
        axis, angle = dcm_to_aa(orientation["dcm"])
    else:
        raise ValueError(f"Unrecognized orientation repr {ori_repr}")
    return axis, angle


# DIRECTION COSINE MATRIX-----------------------------------------------
def mat_is_dcm(mat, tol=1e-14):
    """
    Checks if `mat` is a direction cosine matrix or not.

    Parameters
    ----------
    mat : (3,3) ndarray
        Array to check.
    tol : float
        Absolute tolerance for comparing determinant and orthogonality.

    Returns
    -------
    bool
        ``True`` if `mat` is a direction cosine matrix, ``False`` otherwise.
    """
    return mat_is_rotmat(mat, tol=tol)


def dcm_from_frames(A, B):
    """
    Returns the direction cosine matrix of frame B with respect to frame A,
    where the bases of the frames are resolved in an underlying fixed frame F.
    If F coincides with A, then `A` reduces to an identity matrix.

    Parameters
    ----------
    A : (3,3) ndarray
        The rows of A represent the orthonormal basis vectors of a right-handed 
        frame A resolved in F.

    B : (3,3) ndarray
        The rows of B represent the orthonormal basis vectors of a right-handed
        frame B resolved in F.

    Returns
    -------
    (3,3) ndarray
        The dcm of frame B with respect to frame A.

    """
    return np.dot(B, A.T)


def dcm_rotmat(dcm):
    """
    Returns the rotation matrix corresponding to a direction cosine matrix.

    Let D be the direction cosine matrix of frame B relative to frame A.  The
    rotation matrix R corresponding to D rotates vectors resolved in A such
    that the bases of A align with those of B, when all the bases are resolved
    in A.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.
    """
    return dcm.copy().T


def dcm_shiftmat(dcm, forward=False):
    """
    Returns the shifter matrix corresponding to a direction cosine matrix.

    Let D be the direction cosine matrix of frame B relative to frame A. The
    shifter matrix S corresponding to D transforms components of vectors (or
    tensors) resolved in A to those in B (or vice-versa).

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix.
    forward : bool, optional
       If ``True``, transform components from A to B. If ``False``, transform
       components from B to A. Default is ``False``.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """
    shiftmat = dcm.copy()
    if not forward:
        shiftmat = shiftmat.T
    return shiftmat


def dcm_rotate_vectors(v, dcm):
    """
    Rotates vectors by a direction cosine matrix.

    Let D be the direction cosine matrix of frame B relative to frame A.
    Vectors resolved in A are rotated such that the bases of A align with those
    of B, when all the bases are resolved in A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows.
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
    Transform vector components between frames A and B.

    The orientation of B relative to A is specified by `dcm`. This orientation
    is used for both transformation directions. The vectors remain unchanged;
    only their components are expressed in a different frame.
    
    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows. Expressed
        in A if `forward` is ``True``, or in B otherwise.
    dcm : (3,3) ndarray
        Direction cosine matrix of frame B with respect to frame A.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

    Returns
    -------
    (3,) or (n,3) ndarray
        Components of the same vectors expressed in B if `forward` is ``True``,
        or in A otherwise. Has the same shape as `v`.
    """
    shiftmat = dcm_shiftmat(dcm, forward=forward)
    return np.dot(v, shiftmat.T)


def dcm_shift_tensor2(a, dcm, forward=False):
    """
    Transform second-order tensor components between frames A and B.

    The orientation of B relative to A is specified by `dcm`. This orientation
    is used for both transformation directions. The tensor remains unchanged;
    only its components are expressed in a different frame.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor, resolved in A if `forward` is ``True``, or in B
        otherwise.
    dcm : (3,3) ndarray
        Direction cosine matrix of frame B with respect to frame A.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

    Returns
    -------
    (3,3) ndarray
        Components of the same tensor expressed in B if `forward` is ``True``,
        or in A otherwise.
    """
    shiftmat = dcm_shiftmat(dcm, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def dcm_shift_tensor3(a, dcm, forward=False):
    """
    Transform third-order tensor components between frames A and B.

    The orientation of B relative to A is specified by `dcm`.  This orientation
    is used for both transformation directions. The tensor remains unchanged;
    only its components are expressed in a different frame.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor, resolved in A if `forward` is ``True``, or in B
        otherwise.
    dcm : (3,3) ndarray
        Direction cosine matrix of frame B with respect to frame A.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

    Returns
    -------
    (3,3,3) ndarray
        Components of the same tensor expressed in B if `forward` is ``True``,
        or in A otherwise.
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
        Angle in radian (lies in [0, *pi*]).
    """
    mat = dcm_rotmat(dcm)
    axis, angle = aa_from_rotmat(mat)
    return (axis, angle)


def dcm_to_euler(dcm, seq="XYZ", world=True, tol=1e-14):
    """
    Converts a direction cosine matrix to *Euler angles*.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    tol : float
        Tolerance (must be >= 0) for setting one outer angle to zero.

    Returns
    -------
    (3,) ndarray
        The three Euler angles ``phi`` (rotation about X), ``theta`` (rotation
        about Y), and ``psi`` (rotation about Z).
    """
    mat = dcm_rotmat(dcm)
    euler = euler_factor_rotmat(mat, seq=seq, world=world, tol=tol)
    return euler


def dcm_from_any(orientation):
    """Converts an orientation to direction cosine matrix from any of the
    following: (1) Quaternion, (2) Euler angles, (3) Axis-angle, or (4)
    Direction cosine matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.quat_from_any`.

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
        dcm = orientation["dcm"].copy()
    else:
        raise ValueError(f"Unrecognized orientation repr {ori_repr}")
    return dcm


# EULER ANGLES-----------------------------------------------------------


def euler_rotmat(euler, seq="XYZ", world=True):
    """
    Returns the rotation matrix for a set of Euler angles.

    Parameters
    ----------
    euler : (3,) ndarray
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix
    """
    return _rotmat_euler(euler, seq=seq, world=world)


def euler_factor_rotmat(rotmat, seq="XYZ", world=True, tol=1e-14):
    """Factorize a rotation matrix to obtain the three Euler angles.

    Parameters
    ----------
    rotmat : (3,3) ndarray
        Rotation matrix (must be orthonormal with determinant equal to one).
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    tol : float
        Tolerance (must be >= 0) for setting one outer angle to zero.

    Returns
    -------
    (3,) ndarray
        Euler angles in radian.
    """
    return _factor_rotmat(rotmat, seq=seq, world=world, tol=tol)


def euler_shiftmat(euler, seq="XYZ", world=True, forward=False):
    """Returns the shifter matrix corresponding to a set of Euler angles.

    The orientation of frame B relative to frame A is specified by a set of
    Euler angles. The shifter matrix S transforms components of vectors (or
    tensors) resolved in A to those in B (or vice-versa).

    Parameters
    ----------
    euler : (3,) ndarray
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
       If ``True``, transform components from A to B. If ``False``, transform
       components from B to A. Default is ``False``.

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

    The components of `v` and the rotated vectors are all expressed in the same
    fixed basis.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows.
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
    Transform vector components between frames A and B.

    The orientation of B relative to A is specified by Euler angles. This
    orientation is used for both transformation directions. The vectors remain
    unchanged; only their components are expressed in a different frame.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        Components of a single vector or n vectors stored as rows. Expressed
        in A if `forward` is ``True``, or in B otherwise.
    euler : (3,) ndarray
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

    Returns
    -------
    (3,) or (n,3) ndarray
        Components of the same vectors expressed in B if `forward` is ``True``,
        or in A otherwise. Has the same shape as `v`.
    """
    shiftmat = euler_shiftmat(euler, seq=seq, world=world, forward=forward)
    return np.dot(v, shiftmat.T)


def euler_shift_tensor2(a, euler, seq="XYZ", world=True, forward=False):
    """
    Transform second-order tensor components between frames A and B.

    The orientation of B relative to A is specified by Euler angles. This
    orientation is used for both transformation directions. The tensor remains
    unchanged; only its components are expressed in a different frame.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor. Expressed in A if `forward` is ``True``, or in B
        otherwise.
    euler : (3,) ndarray
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

    Returns
    -------
    (3,3) ndarray
        Components of the same tensor expressed in B if `forward` is ``True``,
        or in A otherwise.
    """
    shiftmat = euler_shiftmat(euler, seq=seq, world=world, forward=forward)
    return np.einsum("ip,jq,pq", shiftmat, shiftmat, a)


def euler_shift_tensor3(a, euler, seq="XYZ", world=True, forward=False):
    """
    Transform third-order tensor components between frames A and B.

    The orientation of B relative to A is specified by Euler angles. This
    orientation is used for both transformation directions. The tensor remains
    unchanged; only its components are expressed in a different frame.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor. Expressed in A if `forward` is ``True``, or in B
        otherwise.
    euler : (3,) ndarray
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, transform components from A to B. If ``False``, transform
        components from B to A. Default is ``False``.

    Returns
    -------
    (3,3,3) ndarray
        Components of the same tensor expressed in B if `forward` is ``True``,
        or in A otherwise.
    """
    shiftmat = euler_shiftmat(euler, seq=seq, world=world, forward=forward)
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
        Unit quaternion in canonical form.
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
        Angle in radian (lies in [0, *pi*]).
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


def euler_to_euler(euler, seq, world, to_seq, to_world, tol=1e-14):
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
        Whether the converted euler angles are with respect to the *world*
        frame or not.
    tol : float
        Tolerance (must be >= 0) for setting one outer angle to zero.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = euler_rotmat(euler, seq=seq, world=world)
    return euler_factor_rotmat(rotmat, seq=to_seq, world=to_world, tol=tol)


def euler_from_any(orientation, to_seq, to_world, tol=1e-14):
    """Converts to Euler angles from any of the following: (1)
    Quaternion, (2) Euler angles, (3) Axis-angle, or (4) Direction cosine
    matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.quat_from_any`.
    to_seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of the returned Euler angles.
    to_world : bool
        Whether the returned Euler angles are with respect to the *world* frame
        or not.
    tol : float
        Tolerance (must be >= 0) for setting one outer angle to zero.

    Returns
    -------
    euler : (3,) ndarray
        Euler angles.
    """
    ori_repr = orientation["repr"]
    if ori_repr == "quat":
        quat = np.array(orientation["quat"])
        euler = quat_to_euler(quat, seq=to_seq, world=to_world, tol=tol)
    elif ori_repr == "euler":
        euler = np.array(orientation["euler"])
        seq = orientation["seq"]
        world = orientation["world"]
        euler = euler_to_euler(euler, seq, world, to_seq, to_world, tol=tol)
    elif ori_repr == "axis_angle":
        axis = np.array(orientation["axis"])
        angle = orientation["angle"]
        euler = aa_to_euler(axis, angle, seq=to_seq, world=to_world, tol=tol)
    elif ori_repr == "dcm":
        euler = dcm_to_euler(orientation["dcm"], seq=to_seq, world=to_world,
                             tol=tol)
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


def _factor_rotmat(rotmat, seq="XYZ", world=True, tol=1e-14):
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
    factors = factor_rotmat_funcs[seq](rotmat, tol=tol)
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


def _factor_rotmat_XYZ(rotmat, tol=1e-14):
    # Recover the middle angle theta
    s = -rotmat[2,0]
    c = math.hypot(rotmat[0,0], rotmat[1,0]) # = abs(cos theta))
    theta = math.atan2(s,c) #Returns in [-pi/2, pi/2]
    #Calculate phi (accounting for gimbal lock)
    if c > tol:
        phi = math.atan2(rotmat[2,1], rotmat[2,2])
    else:
        phi = 0.0
    #Calculate psi using phi
    sphi = math.sin(phi)
    cphi = math.cos(phi)
    psi = math.atan2(sphi*rotmat[0,2]-cphi*rotmat[0,1], 
                     cphi*rotmat[1,1]-sphi*rotmat[1,2])
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


def _factor_rotmat_XZY(rotmat, tol=1e-14):
    #Recover angle psi
    s = rotmat[1,0]
    c = math.hypot(rotmat[0,0], rotmat[2,0]) # = abs(cos psi))
    psi = math.atan2(s,c) #Returns in [-pi/2, pi/2]
    #Calculate phi (accounting for gimbal lock)
    if c > tol:
        phi = math.atan2(-rotmat[1,2], rotmat[1,1])
    else:
        phi = 0.0
    #Calculate theta using phi
    sphi = math.sin(phi)
    cphi = math.cos(phi)
    theta = math.atan2(sphi*rotmat[0,1]+cphi*rotmat[0,2], 
                       sphi*rotmat[2,1]+cphi*rotmat[2,2])
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


def _factor_rotmat_YXZ(rotmat, tol=1e-14):
    #Recover middle angle phi
    s = rotmat[2,1]
    c = math.hypot(rotmat[0,1], rotmat[1,1]) # = abs(cos phi))
    phi = math.atan2(s,c) #Returns in [-pi/2, pi/2]
    #Calculate theta (accounting for gimbal lock)
    if c > tol:
        theta = math.atan2(-rotmat[2,0], rotmat[2,2])
    else:
        theta = 0.0
    #Calculate psi using theta
    stheta = math.sin(theta)
    ctheta = math.cos(theta)
    psi = math.atan2(ctheta*rotmat[1,0]+stheta*rotmat[1,2], 
                     ctheta*rotmat[0,0]+stheta*rotmat[0,2])
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


def _factor_rotmat_YZX(rotmat, tol=1e-14):
    #Recover middle angle psi
    s = -rotmat[0,1]
    c = math.hypot(rotmat[1,1], rotmat[2,1]) # = abs(cos psi))
    psi = math.atan2(s,c) #Returns in [-pi/2, pi/2]
    #Calculate theta (accounting for gimbal lock)
    if c > tol:
        theta = math.atan2(rotmat[0,2], rotmat[0,0])
    else:
        theta = 0.0
    #Calculate phi using theta
    stheta = math.sin(theta)
    ctheta = math.cos(theta)
    phi = math.atan2(stheta*rotmat[1,0]-ctheta*rotmat[1,2], 
                     ctheta*rotmat[2,2]-stheta*rotmat[2,0])
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


def _factor_rotmat_ZXY(rotmat, tol=1e-14):
    #Recover middle angle phi
    s = -rotmat[1,2]
    c = math.hypot(rotmat[0,2], rotmat[2,2]) # = abs(cos phi))
    phi = math.atan2(s,c) #Returns in [-pi/2, pi/2]
    #Calculate psi (accounting for gimbal lock)
    if c > tol:
        psi = math.atan2(rotmat[1,0], rotmat[1,1])
    else:
        psi = 0.0
    #Calculate theta using psi
    spsi = math.sin(psi)
    cpsi = math.cos(psi)
    theta = math.atan2(spsi*rotmat[2,1]-cpsi*rotmat[2,0], 
                       cpsi*rotmat[0,0]-spsi*rotmat[0,1])
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


def _factor_rotmat_ZYX(rotmat, tol=1e-14):
    # Recover the middle angle theta
    s = rotmat[0,2]
    c = math.hypot(rotmat[1,2], rotmat[2,2]) # = abs(cos theta))
    theta = math.atan2(s,c) #Returns in [-pi/2, pi/2]
    #Calculate psi (accounting for gimbal lock)
    if c > tol:
        psi = math.atan2(-rotmat[0,1], rotmat[0,0])
    else:
        psi = 0.0
    #Calculate phi using psi
    spsi = math.sin(psi)
    cpsi = math.cos(psi)
    phi = math.atan2(spsi*rotmat[2,0]+cpsi*rotmat[2,1], 
                     spsi*rotmat[1,0]+cpsi*rotmat[1,1])
    return np.array([phi, theta, psi])
