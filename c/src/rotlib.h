/** @file */
#ifndef  ROTLIB_H
#define ROTLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

/** @note
* Euler angle ranges:                                                         
*                                                                             
* XYZ, ZYX: phi in [-pi, pi],     theta in [-pi/2, pi/2], psi in [-pi, pi]    
* XZY, YZX: phi in [-pi, pi],     theta in [-pi, pi],     psi in [-pi/2, pi/2]
* ZXY, YXZ: phi in [-pi/2, pi/2], theta in [-pi, pi],     psi in [-pi, pi]    
*                                                                             
* Euler angle sequence: XYZ (world). First rotation about X, second rotation  
* about Y, and the third rotation about Z axis of the world(i.e. fixed) frame.
* This is the same as the sequence used in Blender.                           
* In contrast, the XYZ sequence is understood in the Aerospace community as:  
* First rotation about Z-axis, second rotation about Y-axis, and the third    
* rotation about X-axis of the body frame.                                    
*                                                                             
* @see http://www.geometrictools.com/Documentation/EulerAngles.pdf            
*/

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
 *  aligns with another set of mutually orthogoanl unit vectors.
 *
 *  The vectors `v` are rotated using a rotation matrix that, when applied
 *  to the unit vectors `old`, aligns them in the direction of the unit
 *  vectors `new`.
 *
 * @param[in] n Number of rows of `v`.
 * @param[in] v A `n x 3` array whose rows are the vectors to align.
 * @param[in] m Number of rows of `old`, same as the number of rows of `new`.
 * @param[in] old A `m x 3` array, where `m` is at most 3. The rows of `old` form a
 * set of mutually orthogonal unit vectors.
 * @param[in] new A `m x 3` array, where `m` is at most 3. The rows of `new` form a
 *  set of mutually orthogonal unit vectors. The angles between the vectors must
 *  be the same as those of `old`.
 * @param[out] valigned A `n x 3` array whose rows are the aligned vectors.
 */
void align(const int n, const double* v, const int m, const double* old, 
        const double* new, double* valigned);

/**
 * @brief Checks whether a 2D array is a rotation matrix.
 * @param[in] mat A 3x3 array.
 * @returns `true` if `mat` is a rotation matrix, `false` otherwise.
 */
bool mat_is_rotmat(const double* mat);

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
 * @brief Creates a unit quarternion with random orientation.
 * @param [out]  q  The unit quaternion.
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
 * @brief Performs in-place inversion of a quaternion.
 * @param[in,out] q `in:` A quaternion.\n`out:` The quaternion `q` inverted.
 */
void quat_inverted(double q[4]);

/**
 * @brief Performs in-place normalization of a quaternion.
 * @param[in,out] q `in:` A quaternion.\n`out:` The quaternion `q` normalized.
 */
void quat_normalized(double q[4]);

/**
 * @brief Returns `true` if a quaternion is normalized, i.e. if it is a unit
 * quaternion.
 * @param[in] q A quaternion.
 * @returns `true` if `q` is normalized, `false` otherwise.
 */
bool quat_is_normalized(const double q[4]);

/** @brief Calculates the product of two unit quaternions.
 * @param[in] p A quaternion.
 * @param[in] q Another quaternion.
 * @param[out] pq Product of `p` and `q`.
 * @param[in] normalize Whether to normalize the product.
 */
void quat_prod(const double p[4], const double q[4], double pq[4],
        bool normalize);

/** @brief Returns the angle between two _unit_ quaternions.
 * @param[in] p A quaternion.
 * @param[in] q Another quaternion.
 * @returns Angle in radian.
 */
double quat_angle_between(const double p[4], const double q[4]);

/** @brief Interpolates between two _unit_ quaternions.
 * @param[in] q1 A unit quaternion.
 * @param[in] q2 Another unit quaternion.
 * @param[in] t Interpolation factor, `0 <= t <= 1`.
 * @param[out] q Interpolated unit quaternion.
 */
void quat_interpolate(const double q1[4], const double q2[4], const double t,
        double q[4]);

/** @brief Returns the matrix mapping the derivative of a unit quaternion to
 * angular velocity.
 * @param[in] q Derivative of a unit quaternion.
 * @param[out] mat The `3 x 4` angular velocity matrix.
*/
void quat_deriv_to_angvel_mat(const double q[4], double* mat);

/** @brief Calculates the angular velocity from a unit quaternion and its time 
 * derivative.
 * @param[in] q A unit quaternion.
 * @param[in] qdot Time derivative of `q`.
 * @param[out] angvel Angular velocity.
 */
void quat_deriv_to_angvel(const double q[4], const double qdot[4], 
        double angvel[3]);

/** @brief Returns the matrix mapping the angular velocity to the time
 * derivative of a unit quaternion.
 * @param[in] q A unit quaternion.
 * @param[out] mat The `4 x 3` quaternion derivative matrix.
*/
void quat_deriv_from_angvel_mat(const double q[4], double* mat);

/** @brief Calculates the time derivative of a unit quaternion from angular
 * velocity.
 * @param[in] q A unit quaternion.
 * @param[in] angvel Angular velocity.
 * @param[out] qdot Time derivative of `q`.
 */
void quat_deriv_from_angvel(const double q[4], const double angvel[3],
        double qdot[4]);

/** @brief Calculates the rotation matrix corresponding to a unit quaternion.
 * @param[in] q A unit quaternion.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void quat_rotmat(const double q[4], double* rotmat);

/** @brief Calculates the shifter matrix corresponding to a unit quaternion.
 * @param[in] q A unit quaternion.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void quat_shiftmat(const double q[4], const bool forward,
        double* shiftmat);

/** @brief Rotates vectors by a unit quaternion.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to rotate.
 * @param[in] q A unit quaternion.
 * @param[out] vrot `n x 3` array whose rows are the rotated vectors.
 */
void quat_rotate_vectors(const int n, const double* v, const double q[4],
        double* vrot);

/** @brief Given two frames A and B such that the orientation of frame B with
 *  respect to frame A is given by a unit quaternion, shifts vectors from A
 *  to B or B to A.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to shift.
 * @param[in] q A unit quaternion.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] vshift `n x 3` array whose rows are the rotated vectors.
 */
void quat_shift_vectors(const int n, const double* v, const double q[4],
        const bool forward, double* vshift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by a unit quaternion, shifts second-order
 * tensors from A to B or B to A.
 * @param[in] A `3 x 3` second-order tensor.
 * @param[in] q A unit quaternion.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3` shifted second-order tensor.
 */
void quat_shift_tensor2(const double* A, const double q[4],
        const bool forward, double* Ashift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by a unit quaternion, shifts third-order
 * tensors from A to B or B to A.
 * @param[in] A `3 x 3 x 3` third-order tensor.
 * @param[in] q A unit quaternion.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3 x 3` shifted third-order tensor.
 */
void quat_shift_tensor3(const double* A, const double q[4],
        const bool forward, double* Ashift);

/** @brief Converts a unit quaternion to an _axis-angle_ representation.
 * @param[in] q A unit quaternion.
 * @param[out] axis Unit vector along the axis of rotation.
 * @param[out] angle Angle in radian.
 */
void quat_to_aa(const double q[4], double axis[3], double* angle);

/** @brief Converts a unit quaternion to a direction cosine matrix.
 * @param[in] q A unit quaternion.
 * @param[out] dcm `3 x 3` direction cosine matrix. 
 */
void quat_to_dcm(const double q[4], double* dcm);

/** @brief Converts a unit quaternion to an Euler angle sequence.
 * @param[in] q A unit quaternion.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Euler angle sequence.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] euler Euler angles.
 */
void quat_to_euler(const double q[4], enum EulangSeq seq, const bool world,
        double euler[3]);
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

/** @brief Modifies `axis` and `angle` to ensure a right handed rotation with
 * `angle` \f$\in [0, \pi)\f$.
 * @param[in,out] axis Axis of rotation. If this is not a unit vector, set `normalize` to 
 * `true`.
 * @param[in,out] angle Angle in radian.
 * @param[in] normalize Whether to normalize `axis` to a unit vector.
 */
void aa_fix(double axis[3], double* angle, const bool normalize);

/** @brief Generates a random orientation in *axis-angle* representation.
 *
 * The axis is a random vector drawn from a uniform distribution on the surface
 * of a unit sphere. The current implementation in based on the algorithm from
 * Allen & Tildesley p. 349.
 * @param[out] axis Axis of rotation.
 * @param[out] angle Angle in radian.
 */
void aa_rand(double axis[3], double* angle);

/** @brief Calculates the rotation matrix corresponding to an *axis-angle*
 * representation.
 * @param[in] axis Axis of rotation.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void aa_rotmat(const double axis[3], const double angle,
        double* rotmat);

/** @brief Extracts axis and angle from a rotation matrix.
 * @param[in] rotmat `3 x 3` rotation matrix.
 * @param[out] axis Unit vector along the axis.
 * @param[out] angle Angle in radian.
 */
void aa_from_rotmat(const double* rotmat, double axis[3], double* angle);

/** @brief Rotates vectors about `axis` by `angle`.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to rotate.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[out] vrot `n x 3` array whose rows are the rotated vectors.
 */
void aa_rotate_vectors(const int n, const double* v, const double axis[3],
        const double angle, double* vrot);

/** @brief Calculates the shifter matrix corresponding to an _axis-angle_
 * representation.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void aa_shiftmat(const double axis[3], const double angle,
        const bool forward, double* shiftmat);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by an _axis-angle_ representation, shifts
 * vectors from A to B or B to A.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to shift.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] vshift `n x 3` array whose rows are the rotated vectors.
 */
void aa_shift_vectors(const int n, const double* v, 
        const double axis[3], const double angle, const bool forward,
        double* vshift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by _axis-angle_ representation, shifts
 * second-order tensors from A to B or B to A.
 * @param[in] A `3 x 3` second-order tensor.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3` shifted second-order tensor.
 */
void aa_shift_tensor2(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by _axis-angle_ representation, shifts
 * third-order tensors from A to B or B to A.
 * @param[in] A `3 x 3 x 3` third-order tensor.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3 x 3` shifted second-order tensor.
 */
void aa_shift_tensor3(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift);

/** @brief Converts an _axis-angle_ representation to a unit quaternion.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[out] q A unit quaternion.
 */
void aa_to_quat(const double axis[3], const double angle, double q[4]);

/** @brief Converts an _axis-angle_ representation to a direction cosine
 * matrix.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[out] dcm `3 x 3` direction cosine matrix. 
 */
void aa_to_dcm(const double axis[3], const double angle, double* dcm);

/** @brief Converts an _axis-angle_ representation to an Euler angle sequence.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian. `0 <= angle < pi`.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Euler angle sequence.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] euler Euler angles.
 */
void aa_to_euler(const double axis[3], const double angle, 
        enum EulangSeq seq, const bool world, double euler[3]);
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

/** @brief Checks is a `3 x 3` array is a direction cosine matrix.
 * @param[in] mat `3 x 3` array.
 * @returns  `true` if a (3,3) matrix is a direction cosine matrix.
 */ 
bool mat_is_dcm(const double* mat);

/** @brief Returns the direction cosine matrix of axes(i.e. frame) B with
 * respect to axes(i.e. frame) A.
 * @param[in] A `3 x 3` array. The rows of A represent the orthonormal basis
 * vectors of frame A.
 * @param[in] B `3 x 3` array. The rows of B represent the orthonormal basis
 * vectors of frame B.
 * @param[out] dcm `3 x 3` direction cosine matrix.
 */
void dcm_from_axes(const double* A, const double* B, double* dcm);

/** @brief Returns the rotation matrix corresponding to a direction cosine
 * matrix.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[out] rotmat `3 x 3` rotation matrix.
 */
void dcm_rotmat(const double* dcm, double* rotmat);

/** @brief Returns the shifter matrix corresponding to a direction cosine
 * matrix.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] shiftmat `3 x 3` shifter matrix.
 */
void dcm_shiftmat(const double* dcm, const bool forward, double* shiftmat);

/** @brief Rotates vectors by a direction cosine matrix.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to rotate.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[out] vrot `n x 3` array whose rows are the rotated vectors.
 */
void dcm_rotate_vectors(const int n, const double* v, const double* dcm,
        double* vrot);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by a direction cosine matrix, shifts
 * vectors from A to B or B to A.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to shift.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] vshift `n x 3` array whose rows are the rotated vectors.
 */
void dcm_shift_vectors(const int n, const double* v, const double* dcm,
        const bool forward, double* vshift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by a direction cosine matrix, shifts
 * second-order tensors from A to B or B to A.
 * @param[in] A `3 x 3` second-order tensor.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3` shifted second-order tensor.
 */
void dcm_shift_tensor2(const double* A, const double* dcm, const bool forward,
        double* Ashift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by a direction cosine matrix, shifts
 * third-order tensors from A to B or B to A.
 * @param[in] A `3 x 3 x 3` third-order tensor.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3 x 3` shifted third-order tensor.
 */
void dcm_shift_tensor3(const double* A, const double* dcm, const bool forward,
        double* Ashift);

/** @brief Converts a direction cosine matrix to a unit quaternion.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] q A unit quaternion.
 */
void dcm_to_quat(const double* dcm, double q[4]);

/** @brief Converts a direction cosine matrix to an _axis-angle_
 * representation.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] axis Unit vector along the axis.
 * @param[in] angle Angle in radian.
 */
void dcm_to_aa(const double* dcm, double axis[3], double* angle);

/** @brief Converts a direction cosine matrix to Euler angles.
 * @param[in] dcm `3 x 3` direction cosine matrix.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Euler angle sequence.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] euler Euler angles.
 */
void dcm_to_euler(const double* dcm, enum EulangSeq seq, const bool world,
        double euler[3]);
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

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by three Euler angles, shifts
 * vectors from A to B or B to A.
 * @param[in] n Number of rows of `v`.
 * @param[in] v `n x 3` array whose rows are the vectors to shift.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] vshift `n x 3` array whose rows are the rotated vectors.
 */
void euler_shift_vectors(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* vshift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by three Euler angles, shifts second-order
 * tensors from A to B or B to A.
 * @param[in] A `3 x 3` second-order tensor.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3` shifted second-order tensor.
 */
void euler_shift_tensor2(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward, double* Ashift);

/** @brief Given two frames A and B such that the orientation of frame B with
 * respect to frame A is given by three Euler angles, shifts third-order
 * tensors from A to B or B to A.
 * @param[in] A `3 x 3 x 3` third-order tensor.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[in] forward Whether to shift forward, i.e., along the orientation or
 * shift reverse.
 * @param[out] Ashift `3 x 3 x 3` shifted third-order tensor.
 */
void euler_shift_tensor3(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward, double* Ashift);

/** @brief Convert Euler angles to a unit quaternion.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] q A unit quaternion.
 */
void euler_to_quat(const double euler[3], enum EulangSeq seq, const bool world,
        double q[4]);

/** @brief Convert Euler angles to an _axis-angle_ representation.
 * @param[in] euler Euler angles.
 * @param[in] seq `ES_XYZ` | `ES_XZY` | `ES_YXZ` | `ES_YZX` | `ES_ZXY` | `ES_ZYX`.\n
 * Sequence for `euler`.
 * @param[in] world Whether the Euler angles are with respect to the _world_
 * frame or not.
 * @param[out] axis Unit vector along the axis of rotation.
 * @param[out] angle Angle of rotation in radian.
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
 * @param[out] to_euler Converted Euler angles.
 */
void euler_to_euler(const double euler[3], enum EulangSeq seq, const bool world,
        enum EulangSeq to_seq, const bool to_world, double to_euler[3]);
/**
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif /* ROTLIB_H */
