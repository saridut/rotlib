/** @file */
#ifndef  ROTLIB_H
#define ROTLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

/** @brief Euler angle sequences. `ES_ABC` signifies first rotation about the
 * `A` axis, second rotation about the `B` axis, and the third rotation about
 * the `C` axis.
 */
enum EulangSeq {ES_XYZ, ES_XZY, ES_YXZ, ES_YZX, ES_ZXY, ES_ZYX};

/******************************************************************************/
                            /* UTILITY FUNCTIONS */
/******************************************************************************/

/**
 * @defgroup util Utility
 * @{
 */

/**
 * @brief Rotates vectors such that a set of mutually orthogonal unit vectors
 *  aligns with another set of mutually orthogoanl unit vectors. The
 *  two sets of mutually orthogonal unit vectors must contain the same number
 *  of unit vectors.
 *
 *  `v`, `old`, and `new` store their constituent vectors row-wise, with
 *  all vectors resolved in the same fixed underlying basis.
 *
 *  The vectors in `v` are rotated using a rotation matrix R such that
 *  `new[i,:] = old[i,:] @ R.T`, for each unit vector `i` (stored in the `i`th
 *  row), with the single vector (`n=1`) case treated equivalently.
 *
 *  If `old` is a single vector, a unique rotation cannot be defined.  In this
 *  case, the shortest rotation mapping `old` to `new` is chosen. For
 *  anti-parallel directions, this is a half-turn about a chosen axis
 *  perpendicular to `old`.
 *
 * `valigned` must not overlap with `v`.
 *
 * @param[in] m Number of rows of `v`.
 * @param[in] v A `m x 3` array whose rows are the vectors to align.
 * @param[in] n Number of rows of `old`, same as the number of rows of `new`.
 * @param[in] old A `n x 3` array, where `n` is at most 3. The rows of `old`
 *  form a set of mutually orthogonal unit vectors. If `n = 3`, the three unit
 *  vectors must form a right-handed basis.
 * @param[in] new A `n x 3` array, where `n` is at most 3. The rows of `new`
 * form a set of mutually orthogonal unit vectors. If `n = 3`, the three unit
 * vectors must form a right-handed basis.
 * @param[out] valigned A `m x 3` array whose rows are the aligned vectors. The
 * components in `out` are resolved in the same underlying basis as the vectors
 *  in `v`.
 */
void align(const int m, const double* v, const int n, const double* old, 
        const double* new, double* valigned);

/**
 * @brief Checks whether a 2D array is a rotation matrix.
 * @param[in] mat A 3x3 array.
 * @param[in] tol Absolute tolerance for comparing determinant and orthogonality.
 * @returns `true` if `mat` is a rotation matrix, `false` otherwise.
 */
bool mat_is_rotmat(const double* mat, const double tol);

/**
 * @}
 */

/******************************************************************************/
                               /* QUATERNIONS */
/******************************************************************************/

/**
 * @defgroup quat Quaternion
 * @{
 */

/**
 * @brief Modifies a non-canonical unit quaternion in-place to bring it to the
 * canonical form.  No modifications are performed if the input is already in
 * canonical form.
 *
 * @param [in out]  q  On entry, a unit quaternion. If `q` is not in canonical
 * form, on return, `q` is modified to bring it to the canonical form, else `q`
 * is left as is.
 */
void quat_canonicalized(double q[4]);

/**
 * @brief Creates a unit quarternion from a uniform distribution on the space of
    rotations in 3D.
 * @param [out]  q  A unit quaternion in canonical form.
 */
void quat_rand(double q[4]);

/**
 * @brief Creates an identity quarternion.
 * @param [out] q Identity quaternion.
 */
void quat_identity(double q[4]);

/**
 * @brief Performs in-place conjugation of a quaternion.
 * @param[in,out] q `in:` A quaternion.\n`out:` The quaternion `q` conjugated.
 */
void quat_conjugated(double q[4]);

/**
 * @brief Performs in-place inversion of a quaternion. In case of a zero
 * quaternion, `q` is not changed and `EXIT_FAILURE` is returned.
 * @param[in,out] q `in:` A non-zero quaternion.\n`out:` The quaternion 
 * `q` inverted.
 * @returns `EXIT_FAILURE` if it is a zero-quaternion, else `EXIT_SUCCESS`.
 */
int quat_inverted(double q[4]);

/**
 * @brief Performs in-place normalization of a quaternion. In case of a zero
 * quaternion, `q` is not changed and `EXIT_FAILURE` is returned.
 * @param[in,out] q `in:` A non-zero quaternion.\n`out:` The normalized
 * quaternion `q` (in canonical form if `canonicalize` is ``true``).
 * @param[in] canonicalize  Whether to canonicalize the output.
 * @returns `EXIT_FAILURE` if it is a zero-quaternion, else `EXIT_SUCCESS`.
 */
int quat_normalized(double q[4], const bool canonicalize);

/**
 * @brief Returns `true` if a quaternion is normalized, i.e. if it is a unit
 * quaternion, within a tolerance `tol`; `false` otherwise.
 * @param[in] q A quaternion.
 * @param[in] tol Absolute allowable deviation of the norm from one.
 * @returns `true` if `q` is normalized within tolerance, `false` otherwise.
 */
bool quat_is_normalized(const double q[4], const double tol);

/** @brief Calculates the product of two quaternions following Hamilton's
 * convention. In case of unit quaternions describing rotations in a common
 * frame, `q`'s rotation is applied first, followed by that of `p`.
 *
 * `pq` must not overlap with `p` or `q`.
 *
 * @param[in] p A quaternion.
 * @param[in] q Another quaternion.
 * @param[out] pq Product of `p` and `q`.
 */
void quat_prod(const double p[4], const double q[4], double pq[4]);

/** @brief Returns the smallest rotation angle between the orientations
 * represented by two _unit_ quaternions.
 * @param[in] p A unit quaternion.
 * @param[in] q Another unit quaternion.
 * @returns Relative rotation angle in radians, in the interval [0, *pi*]. 
 * Opposite-sign quaternions representing the same orientation have zero
 * separation.
 */
double quat_angle_between(const double p[4], const double q[4]);

/** @brief Performs spherical linear interpolation between two unit quaternions
 * following the shortest path. The second endpoint may be negated to select
 * the shortest path. At a relative orientation (physical) of exactly *pi*, two
 * equally short paths exist; the supplied endpoint signs determine which is
 *  selected.
 *
 * `q` may coincide exactly with `q1`, but must not partially overlap it.
 *
 * @param[in] q1 A unit quaternion.
 * @param[in] q2 Another unit quaternion.
 * @param[in] t Interpolation factor, `0 <= t <= 1` specifying the
 * interpolation point..
 * @param[out] q Interpolated unit quaternion.
 */
void quat_interpolate(const double q1[4], const double q2[4], const double t,
        double q[4]);

/** @brief Returns the matrix mapping the derivative of a unit quaternion to
 * angular velocity.
 *
 *  Let A (world frame) and B (body frame) be reference frames such that an
 *  active rotation in A by a unit quaternion ``q`` aligns its basis vectors
 *  with those of frame B.  Then the matrix M returned is such that ``omega_i =
 *  M_ij * qdot_j``, where ``omega_i`` is the angular velocity of B relative to
 *  A resolved in A (or B) and ``qdot_j`` is the time derivative of ``q``.
 *  
 *
 * `mat` must not overlap with `q`.
 *
 * @param[in] q A unit quaternion.
 * @param[in] world If ``true``, the angular velocity is resolved in frame A
 * (world frame), else, it is resolved in frame B (body frame).
 * @param[out] mat The `3 x 4` angular velocity matrix.
*/
void quat_deriv_to_angvel_mat(const double q[4], const bool world, double* mat);

/** @brief Calculates the angular velocity from a unit quaternion and its time
 *  derivative.
 *
 *  Let A (world frame) and B (body frame) be reference frames such that an
 *  active rotation in A by a unit quaternion ``q`` aligns its basis vectors
 *  with those of frame B. Let ``qdot`` be the time derivative of ``q``. The
 *  quantity returned is the angular velocity of B relative to A resolved in A
 *  (or B).
 *
 * `angvel` must not overlap with `qdot`.
 *
 * @param[in] q A unit quaternion.
 * @param[in] qdot Time derivative of `q`.  In general, this is not a unit
 * vector. Note that `qdot` is orthogonal to `q`, i.e. ``q . qdot = 0``.
 * @param[in] world If ``true``, the angular velocity is resolved in frame A
 * (world frame), else, it is resolved in frame B (body frame).
 * @param[out] angvel Angular velocity.
 */
void quat_deriv_to_angvel(const double q[4], const double qdot[4], 
        const bool world, double angvel[3]);

/** @brief Returns the matrix mapping the angular velocity to the time
 * derivative of a unit quaternion.
 *
 * Let A (world frame) and B (body frame) be reference frames such that an
 * active rotation in A by a unit quaternion ``q`` aligns its basis vectors
 * with those of frame B. Then the matrix M returned is such that ``qdot_i =
 * M_ij * omega_j``, where ``omega_j`` is the angular velocity of B relative
 * to A resolved in A (or B) and ``qdot_i`` is the time derivative of ``q``.
 *
 * `mat` must not overlap with `q`.
 *
 * @param[in] q A unit quaternion.
 * @param[in] world If ``true``, the angular velocity is resolved in frame A
 * (world frame), else, it is resolved in frame B (body frame).
 * @param[out] mat The `4 x 3` quaternion derivative matrix.
*/
void quat_deriv_from_angvel_mat(const double q[4], const bool world, double* mat);

/** @brief Calculates the time derivative of a unit quaternion from angular
 * velocity.
 *
 * Let A (world frame) and B (body frame) be reference frames such that an
 * active rotation in A by a unit quaternion ``q`` aligns its basis vectors
 * with those of frame B. Let ``ang_vel`` be the angular velocity of B
 * relative to A resolved in A (or B). The quantity returned is the time
 * derivative of ``q``.
 *
 * `qdot` must not overlap with `angvel`.
 *
 * @param[in] q A unit quaternion. An active rotation of frame A by ``q``
 * aligns it with frame B.
 * @param[in] angvel Angular velocity of frame B relative to frame A resolved
 * in A if ``world=true``, or in B otherwise.
 * @param[in] world If ``true``, the angular velocity is resolved in frame A
 * (world frame), else, it is resolved in frame B (body frame).
 * @param[out] qdot Time derivative of `q`.
 */
void quat_deriv_from_angvel(const double q[4], const double angvel[3],
        const bool world, double qdot[4]);

/** @brief Calculates the rotation matrix corresponding to a unit quaternion.
 *
 * For a column vector v of length 3, the rotated vector ``v' = R v``, where
 * ``R`` is the returned matrix.  The components of ``v`` and ``v'`` are
 * expressed in the same fixed orthonormal basis. For a row vector, the
 * equivalent rotated vector ``v' = v R^T``.
 * 
 * @param[in] q A unit quaternion.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void quat_rotmat(const double q[4], double* rotmat);

/** @brief Calculates the shifter matrix corresponding to a unit quaternion.
 * @param[in] q A unit quaternion.
 *
 * The orientation of frame B relative to frame A is specified by a unit
 * quaternion ``q``, such that an active rotation in A of its basis vectors
 * aligns them with those of B. The shifter matrix ``S`` transforms the
 * components of a fixed vector resolved in A to those resolved in B or
 * vice-versa according to ``v_B_i = S_ij * v_A_j`` or ``v_A_i = S_ij *
 * v_B_j``.
 *
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void quat_shiftmat(const double q[4], const bool forward, double* shiftmat);

/** @brief Rotates vectors by a unit quaternion. The components of `v` and the
 * rotated vectors are all expressed in the same fixed basis.
 *
 * `vrot` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the components of vectors to rotate.
 * @param[in] q A unit quaternion.
 * @param[out] vrot `n x 3` array whose rows are the components of the rotated vectors.
 */
void quat_rotate_vectors(const int n, const double* v, const double q[4],
        double* vrot);

/** @brief Transform vector components between frames A and B.
 *
 *  The orientation of B relative to A is specified by a unit quaternion `q`.
 *  This orientation is used for both transformation directions.  The vectors
 *  remain unchanged; only their components are expressed in a different frame.
 * 
 *
 * `vshift` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to shift. Expressed in
 * A if `forward` is ``true``, or in B otherwise.
 * @param[in] q A unit quaternion expressed in A, regardless of `forward`.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] vshift `n x 3` array whose rows are the components of the same vectors
 *  expressed in B if `forward` is ``true``, or in A otherwise.
 */
void quat_shift_vectors(const int n, const double* v, const double q[4],
        const bool forward, double* vshift);

/** @brief Transform second-order tensor components between frames A and B.
 *  
 *  The orientation of B relative to A is specified by a unit quaternion `q`.
 *  This orientation is used for both transformation directions.  The tensor
 *  remains unchanged; only its components are expressed in a different frame.
 *  
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A `3 x 3` second-order tensor. Expressed in A if `forward` is
 * ``true``, or in B otherwise.
 * @param[in] q A unit quaternion expressed in A, regardless of `forward`.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] Ashift `3 x 3` matrix containing the components of the same 
 * second-order tensor expressed in B if `forward` is ``true``, or in A otherwise.
 */
void quat_shift_tensor2(const double* A, const double q[4],
        const bool forward, double* Ashift);

/** @brief Transform third-order tensor components between frames A and B.
 *  
 *  The orientation of B relative to A is specified by a unit quaternion `q`.
 *  This orientation is used for both transformation directions.  The tensor
 *  remains unchanged; only its components are expressed in a different frame.
 *  
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A `3 x 3 x 3` third-order tensor. Expressed in A if `forward` is
 * ``true``, or in B otherwise.
 * @param[in] q A unit quaternion expressed in A, regardless of `forward`.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] Ashift `3 x 3 x 3` matrix containing the components of the same 
 * second-order tensor expressed in B if `forward` is ``true``, or in A otherwise.
 */
void quat_shift_tensor3(const double* A, const double q[4],
        const bool forward, double* Ashift);

/** @brief Converts a unit quaternion to an _axis-angle_ representation. Note
 * that an identity quaternion will return axis [1,0,0].
 *
 * `axis` and `angle` must not overlap with each other.
 *
 * @param[in] q A unit quaternion.
 * @param[out] axis Unit vector along the axis of rotation.
 * @param[out] angle Angle in radian (lies in [0, *pi*]).
 */
void quat_to_aa(const double q[4], double axis[3], double* angle);

/** @brief Converts a unit quaternion to a direction cosine matrix.
 * @param[in] q A unit quaternion.
 * @param[out] dcm `3 x 3` direction cosine matrix. 
 */
void quat_to_dcm(const double q[4], double* dcm);

/** @brief Converts a unit quaternion to an Euler angle sequence.
 *
 * @param[in] q A unit quaternion.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Euler angle sequence.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] tol Tolerance (must be >= 0) for setting one outer angle to zero.
 * @param[out] euler Euler angles.
 */
void quat_to_euler(const double q[4], enum EulangSeq seq, const bool world,
        const double tol, double euler[3]);
/**
 *@}
 */

/******************************************************************************/
                                /* AXIS ANGLE */
/******************************************************************************/

/**
 * @defgroup aa Axis Angle
 * @{
 */

/** @brief For a given orientation in *axis-angle* representation, returns a
 * modified `axis` and `angle` to ensure that `angle` lies in [0, *pi*] while
 * preserving the rotation.
 *
 * `axis` and `angle` must not overlap with each other.
 *
 * @param[in,out] axis `in:` A non-zero 3-vector representing the axis of
 * rotation.\n`out:` If `normalize_axis=true`, the normalized `axis`. If
 * `normalize_axis=false`, it has the same magnitude as `axis` and may need to
 * be normalized before using other axis-angle functions.
 * @param[in,out] angle `in:` Angle in radian.\n`out:` Modified angle in radian
 * (lies in [0, *pi*]).
 * @param[in] normalize_axis Whether to normalize `axis` to a unit vector.
 * @returns `EXIT_FAILURE` if `axis` is a zero-vector, `EXIT_SUCCESS`
 * otherwise. No changes to `axis` or `angle` is made is case of
 * `EXIT_FAILURE`.
 */
int aa_fix(double axis[3], double* angle, const bool normalize_axis);

/** @brief Generates a random orientation in *axis-angle* representation drawn
 * from a uniform distribution on the space of rotations in 3D.
 *
 * `axis` and `angle` must not overlap with each other.
 *
 * @param[out] axis A 3-vector representing a unit vector along the axis of
 *  rotation.
 * @param[out] angle Angle in radian (lies in [0, *pi*]).
 */
void aa_rand(double axis[3], double* angle);

/** @brief Calculates the matrix for an active rotation about `axis` through
 * `angle`, with positive angles following the right-hand rule.
 *
 *  For a column vector v of length 3, the rotated vector v' = R v, where R is
 *  the returned matrix.  The components of v, v', and `axis` are expressed in
 *  the same fixed orthonormal basis. For a row vector, the equivalent rotated
 *  vector v' = v R^T.
 * 
 *
 * `rotmat` must not overlap with `axis`.
 *
 * @param[in] axis A 3-vector representing a unit vector along the axis of
 *  rotation.
 * @param[in] angle Angle in radian. Positive angles follow the right-hand
 * rule.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void aa_rotmat(const double axis[3], const double angle,
        double* rotmat);

/** @brief Extracts axis and angle from a rotation matrix.
 *
 *  If `rotmat` is an identity matrix, the axis returned is [1,0,0]. If the
 *  angle is exactly *pi*, either sign of the axis is valid.
 *
 * `axis` and `angle` must not overlap with each other or with `rotmat`.
 *
 * @param[in] rotmat `3 x 3` rotation matrix. Must be orthonormal with
 * determinant equal to one.
 * @param[out] axis [3,] Unit vector along the axis.
 * @param[out] angle Angle of roation in radian within [0, *pi*].
 */
void aa_from_rotmat(const double* rotmat, double axis[3], double* angle);

/** @brief Rotates vectors about `axis` by `angle`. The components of `axis`,
 * `v`, and the rotated vectors are all expressed in the same fixed basis.
 *
 * `vrot` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to rotate.
 * @param[in] axis [3,] Unit vector along the axis.
 * @param[in] angle Angle in radian. Positive angles follow the right-hand
 * rule.
 * @param[out] vrot `n x 3` array whose rows are the rotated vectors.
 */
void aa_rotate_vectors(const int n, const double* v, const double axis[3],
        const double angle, double* vrot);

/** @brief Calculates the shifter matrix corresponding to an _axis-angle_
 * representation.
 * @param[in] axis [3,] Unit vector along the axis.
 * @param[in] angle Angle in radian. Positive angles follow the right-hand rule.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void aa_shiftmat(const double axis[3], const double angle,
        const bool forward, double* shiftmat);

/** @brief Transform vector components between frames A and B.
 *
 *  The orientation of B relative to A is specified by `axis` and `angle`,
 *  with `axis` expressed in A and positive angles following the right-hand
 *  rule. This orientation is used for both transformation directions.
 *  The vectors remain unchanged; only their components are expressed in a
 *  different frame.
 *
 * `vshift` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the components of vectors to
 *  shift. Expressed in A if `forward` is ``true``, or in B otherwise.
 * @param[in] axis [3,] Unit vector expressed in A, regardless of `forward`.
 * @param[in] angle Angle in radian specifying the orientation of B relative to
 * A, in radians. Positive angles follow the right-hand rule.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] vshift `n x 3` array whose rows are the components of the same
 * vectors expressed in B if `forward` is ``true``, or in A otherwise.
 */
void aa_shift_vectors(const int n, const double* v, 
        const double axis[3], const double angle, const bool forward,
        double* vshift);

/** @brief Transforms second-order tensor components between frames A and B.
 *
 * The orientation of B relative to A is specified by `axis` and `angle`,
 * with `axis` expressed in A and positive angles following the right-hand
 * rule. This orientation is used for both transformation directions.  The
 * tensor remains unchanged; only its components are expressed in a different
 * frame.
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A Components of a `3 x 3` second-order tensor, expressed in A if
 * `forward` is ``true``, or in B otherwise.
 * @param[in] axis [3,] Unit vector along the axis expressed in A, regardless
 * of `forward`.
 * @param[in] angle Angle specifying the orientation of B relative to A, in
 * radians. Positive angles follow the right-hand rule.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] Ashift `3 x 3` array containing the components of the
 * second-order tensor expressed in B if `forward` is ``true``, or in A
 * otherwise.
 */
void aa_shift_tensor2(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift);

/** @brief Transforms third-order tensor components between frames A and B.
 *
 * The orientation of B relative to A is specified by `axis` and `angle`,
 * with `axis` expressed in A and positive angles following the right-hand
 * rule. This orientation is used for both transformation directions.  The
 * tensor remains unchanged; only its components are expressed in a different
 * frame.
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A Components of a `3 x 3 x 3` third-order tensor, expressed in A if
 * `forward` is ``true``, or in B otherwise.
 * @param[in] axis [3,] Unit vector along the axis expressed in A, regardless
 * of `forward`.
 * @param[in] angle Angle specifying the orientation of B relative to A, in
 * radians. Positive angles follow the right-hand rule.
 * @param[in] forward If ``true``, transform components from A to B. If
 * ``false``, transform components from B to A.
 * @param[out] Ashift `3 x 3 x 3` array containing the components of the
 * third-order tensor expressed in B if `forward` is ``true``, or in A
 * otherwise.
 */
void aa_shift_tensor3(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift);

/** @brief Converts an _axis-angle_ representation to a unit quaternion.
 *
 * `q` must not overlap with `axis`.
 *
 * @param[in] axis [3,] Unit vector along the axis.
 * @param[in] angle Angle in radian.
 * @param[out] q A unit quaternion in the canonical form.
 */
void aa_to_quat(const double axis[3], const double angle, double q[4]);

/** @brief Converts an _axis-angle_ representation to a direction cosine
 * matrix.
 * @param[in] axis [3,] Unit vector along the axis.
 * @param[in] angle Angle in radian.
 * @param[out] dcm `3 x 3` direction cosine matrix. 
 */
void aa_to_dcm(const double axis[3], const double angle, double* dcm);

/** @brief Converts an _axis-angle_ representation to an Euler angle sequence.
 * @param[in] axis [3,] Unit vector along the axis.
 * @param[in] angle Angle in radian.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Euler angle sequence.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] tol Tolerance (must be >= 0) for setting one outer angle to zero.
 * @param[out] euler Euler angles.
 */
void aa_to_euler(const double axis[3], const double angle, enum EulangSeq seq,
        const bool world, const double tol, double euler[3]);
/**
 * @}
 */

/******************************************************************************/
                         /* DIRECTION COSINE MATRIX */
/******************************************************************************/

/**
 * @defgroup dcm Direction cosine matrix
 * @{
 */

/** @brief Checks if a `3 x 3` array is a direction cosine matrix.
 * @param[in] mat `3 x 3` array.
 * @param[in] tol Absolute tolerance for comparing determinant and orthogonality.
 * @returns  `true` if `mat` is a direction cosine matrix, `false` otherwise.
 */ 
bool mat_is_dcm(const double* mat, const double tol);

/** @brief Returns the direction cosine matrix of frame B with respect to frame A,
 *  where the bases of the frames are resolved in an underlying fixed frame F.
 *  If F coincides with A, then `A` reduces to an identity matrix.
 *
 * `dcm` must not overlap with `A` or `B`.
 *
 * @param[in] A `3 x 3` array. The rows of A represent the orthonormal basis
 * vectors of a right-handed frame A resolved in F.
 * @param[in] B `3 x 3` array. The rows of B represent the orthonormal basis
 * vectors of a right-handed frame B resolved in F.
 * @param[out] dcm `3 x 3` direction cosine matrix of frame B with respect to
 * frame A.
 */
void dcm_from_axes(const double* A, const double* B, double* dcm);

/** @brief Returns the rotation matrix corresponding to a direction cosine
 * matrix.
 *
 *  Let D be the direction cosine matrix of frame B relative to frame A.  The
 *  rotation matrix R corresponding to D rotates vectors resolved in A such
 *  that the bases of A align with those of B, when all the bases are resolved
 *  in A.
 *
 * `rotmat` must not overlap with `dcm`.
 *
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void dcm_rotmat(const double* dcm, double* rotmat);

/** @brief Returns the shifter matrix corresponding to a direction cosine
 * matrix.
 *
 * Let D be the direction cosine matrix of frame B relative to frame A. The
 * shifter matrix S corresponding to D transforms components of vectors (or
 * tensors) resolved in A to those in B (or vice-versa).
 *
 * `shiftmat` must not overlap with `dcm`, regardless of `forward`.
 *
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] forward If `true`, transform components from A to B. If
 * `false`, transform components from B to A.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void dcm_shiftmat(const double* dcm, const bool forward, double* shiftmat);

/** @brief Rotates vectors by a direction cosine matrix.
 *
 * Let D be the direction cosine matrix of frame B relative to frame A.
 * Vectors resolved in A are rotated such that the bases of A align with
 * those of B, when all the bases are resolved in A.
 *
 * `vrot` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the components of vectors to
 * rotate.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[out] vrot `n x 3` array whose rows are the components of rotated
 * vectors.
 */
void dcm_rotate_vectors(const int n, const double* v, const double* dcm,
        double* vrot);

/** @brief Transform vector components between frames A and B.
 *
 * The orientation of B relative to A is specified by `dcm`. This orientation
 * is used for both transformation directions. The vectors remain unchanged;
 * only their components are expressed in a different frame.
 *
 * `vshift` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the components of vectors. Expressed
 *      in A if `forward` is `true`, or in B otherwise.
 * @param[in] dcm `3 x 3` direction cosine matrix of frame B with respect to
 * frame A.
 * @param[in] forward If `true`, transform components from A to B. If
 * `false`, transform components from B to A.
 * @param[out] vshift `n x 3` array whose rows are the components of the same
 * vectors expressed in B if `forward` is `true`, or in A otherwise.
 */
void dcm_shift_vectors(const int n, const double* v, const double* dcm,
        const bool forward, double* vshift);

/** @brief Transform second-order tensor components between frames A and B.
 *
 * The orientation of B relative to A is specified by `dcm`. This orientation
 * is used for both transformation directions. The tensor remains unchanged;
 * only its components are expressed in a different frame.
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A `3 x 3` second-order tensor. Resolved in A if `forward` is
 * `true`, or in B otherwise.
 * @param[in] dcm `3 x 3` direction cosine matrix  of frame B with respect to
 * frame A.
 * @param[in] forward If `true`, transform components from A to B. If
 * `false`, transform components from B to A.
 * @param[out] Ashift `3 x 3` second-order tensor whose components are
 * expressed in B if `forward` is `true`, or in A otherwise.
 */
void dcm_shift_tensor2(const double* A, const double* dcm, const bool forward,
        double* Ashift);

/** @brief Transform third-order tensor components between frames A and B.
 *
 * The orientation of B relative to A is specified by `dcm`. This orientation
 * is used for both transformation directions. The tensor remains unchanged;
 * only its components are expressed in a different frame.
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A `3 x 3 x 3` third-order tensor. Resolved in A if `forward` is
 * `true`, or in B otherwise.
 * @param[in] dcm `3 x 3` direction cosine matrix  of frame B with respect to
 * frame A.
 * @param[in] forward If `true`, transform components from A to B. If
 * `false`, transform components from B to A.
 * @param[out] Ashift `3 x 3 x 3` third-order tensor whose components are
 * expressed in B if `forward` is `true`, or in A otherwise.
 */
void dcm_shift_tensor3(const double* A, const double* dcm, const bool forward,
        double* Ashift);

/** @brief Converts a direction cosine matrix to a unit quaternion.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[out] q A unit quaternion.
 */
void dcm_to_quat(const double* dcm, double q[4]);

/** @brief Converts a direction cosine matrix to an _axis-angle_
 * representation.
 *
 * `axis` and `angle` must not overlap with each other.
 *
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[out] axis [3,] Unit vector along the axis.
 * @param[out] angle Angle in radian (lies in [0, *pi*]).
 */
void dcm_to_aa(const double* dcm, double axis[3], double* angle);

/** @brief Converts a direction cosine matrix to _Euler angles_.
 *
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Euler angle sequence.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] tol Tolerance (must be >= 0) for setting one outer angle to zero.
 * @param[out] euler Euler angles.
 */
void dcm_to_euler(const double* dcm, enum EulangSeq seq, const bool world,
        const double tol, double euler[3]);
/**
 * @}
 */

/******************************************************************************/
                               /* EULER ANGLES */
/******************************************************************************/

/**
 * @defgroup euler Euler angles
 * @{
 */

/** @brief Calculates the rotation matrix for a set of Euler angles.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void euler_rotmat(const double euler[3], enum EulangSeq seq, const bool world,
        double* rotmat);

/** @brief Calculates the shifter matrix for a set of Euler angles.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void euler_shiftmat(const double euler[3], enum EulangSeq seq, const bool world,
        const bool forward, double* shiftmat);

/** @brief Rotates vectors with a set of Euler angles.
 *
 * `vrot` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to rotate.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] vrot `n x 3` array whose rows are the rotated vectors.
 */
void euler_rotate_vectors(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, double* vrot);

/** @brief Transform vector components between frames A and B.
 *
 * The orientation of B relative to A is specified by Euler angles. This
 * orientation is used for both transformation directions. The vectors remain
 * unchanged; only their components are expressed in a different frame.
 *
 * `vshift` must not overlap with `v`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the components of vectors.
 * Expressed in A if `forward` is `true`, or in B otherwise.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward If `true`, transform components from A to B. If `false`,
 * transform components from B to A.
 * @param[out] vshift `n x 3` array. Components of the same vectors expressed
 * in B if `forward` is `true`, or in A otherwise.
 */
void euler_shift_vectors(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* vshift);

/** @brief Transform second-order tensor components between frames A and B.
 *
 *  The orientation of B relative to A is specified by Euler angles. This
 *  orientation is used for both transformation directions. The tensor remains
 *  unchanged; only its components are expressed in a different frame.
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A `3 x 3` second-order tensor. Expressed in A if `forward` is
 * `true`, or in B otherwise.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward If `true`, transform components from A to B. If `false`,
 * transform components from B to A.
 * @param[out] Ashift `3 x 3` array containing the components of the same
 * tensor expressed in B if `forward` is `true`, or in A otherwise.
 */
void euler_shift_tensor2(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward, double* Ashift);

/** @brief Transform third-order tensor components between frames A and B.
 *
 *  The orientation of B relative to A is specified by Euler angles. This
 *  orientation is used for both transformation directions. The tensor remains
 *  unchanged; only its components are expressed in a different frame.
 *
 * `Ashift` must not overlap with `A`.
 *
 * @param[in] A `3 x 3 x 3` third-order tensor. Expressed in A if `forward` is
 * `true`, or in B otherwise.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward If `true`, transform components from A to B. If `false`,
 * transform components from B to A.
 * @param[out] Ashift `3 x 3 x 3` array containing the components of the same
 * tensor expressed in B if `forward` is `true`, or in A otherwise.
 */
void euler_shift_tensor3(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward, double* Ashift);

/** @brief Convert Euler angles to a unit quaternion.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] q A unit quaternion in canonical form.
 */
void euler_to_quat(const double euler[3], enum EulangSeq seq, const bool world,
        double q[4]);

/** @brief Convert Euler angles to an _axis-angle_ representation.
 *
 * `axis` and `angle` must not overlap with each other.
 *
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] axis [3,] Unit vector along the axis of rotation.
 * @param[out] angle Angle of rotation in radian (lies in [0, *pi*]).
 */
void euler_to_aa(const double euler[3], enum EulangSeq seq, const bool world,
        double axis[3], double* angle);

/** @brief Convert Euler angles to a direction cosine matrix.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] dcm `3 x 3` direction cosine matrix.
 */
void euler_to_dcm(const double euler[3], enum EulangSeq seq, const bool world,
        double* dcm);

/** @brief
 */
/** @brief Converts Convert one set of Euler angles to another.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] to_seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Convert `euler` with sequence `seq` to `to_euler` with sequence `to_seq`.
 * @param[in] to_world Whether the converted euler angles are with respect to
 * the _world_ frame or not.
 * @param[in] tol Tolerance (must be >= 0) for setting one outer angle to zero.
 * @param[out] to_euler Converted Euler angles.
 */
void euler_to_euler(const double euler[3], enum EulangSeq seq, const bool world,
        enum EulangSeq to_seq, const bool to_world, const double tol, 
        double to_euler[3]);
/**
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif /* ROTLIB_H */
