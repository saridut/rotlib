module rotlib_m

use constants_m, only: math_pi, math_pi_2
use utils_math_m

implicit none

public

private :: factor_rotmat, shift_tensor2, shift_tensor3,      &
    rotmat_xyz, rotmat_xzy, rotmat_yxz,                      &
    rotmat_yzx, rotmat_zxy, rotmat_zyx,                      &
    factor_rotmat_xyz, factor_rotmat_xzy, factor_rotmat_yxz, &
    factor_rotmat_yzx, factor_rotmat_zxy, factor_rotmat_zyx

contains

!*****************************************************************************!
!                               UTILITY FUNCTIONS                             !
!*****************************************************************************!

subroutine align(v, old, new, valigned)
    !! Rotate vectors such that a set of at most three mutually orthogonal unit
    !! vectors aligns with another set of mutually orthogonal unit vectors. The
    !! two sets of mutually orthogonal unit vectors must contain the same number
    !! of unit vectors.
    !!
    !! `v`, `old`, and `new` store their constituent vectors column-wise, with
    !! all vectors resolved in the same fixed underlying basis.
    !!
    !! The vectors in `v` are rotated using a rotation matrix `R` such that
    !! `new(:,i) = R @ old(:,i)`, for each unit vector `i` (stored in the
    !! `i` th column), with the single vector (`n=1`) case treated equivalently.
    !! Here `@` represents matrix multiplication.
    !!
    !! If `old` is a single vector, a unique rotation cannot be defined. In this
    !! case, the shortest rotation mapping `old` to `new` is chosen. For
    !! anti-parallel directions, this is a half-turn about a chosen axis
    !! perpendicular to `old`.

    real(rp), dimension(:,:), intent(in) :: v
        !! `(3,m)`. Vectors to align.
    real(rp), dimension(:,:), intent(in) :: old
        !! `(3,n)`. A set of at most three mutually orthogonal unit vectors.
        !! `n` must be one of 1, 2 or 3. If `n = 3`, the three unit vectors
        !! must form a right-handed basis.
    real(rp), dimension(:,:), intent(in) :: new
        !! `(3,n)`. Another set of mutually orthogonal unit vectors of the same
        !! shape as `old`. If `n = 3`, the three unit vectors must form a right-handed
        !! basis.
    real(rp), dimension(:,:), intent(out) :: valigned
        !! `(3,m)`. Aligned vectors, of the same shape as `v`. The components in
        !! `out` are resolved in the same underlying basis as the vectors in `v`.
    real(rp) :: axis(3), vtmp(3), z_old(3), z_new(3), axes_old(3,3), &
                axes_new(3,3), mat(3,3)
    real(rp) :: cos_angle, sin_angle, angle
    integer :: n, estat, min_idx

    if ( any(shape(old) /= shape(new)) ) then
        write(*,*) 'Shape mismatch between `old` and `new`.'
        write(*,*) 'Shape of `old` = ', shape(old)
        write(*,*) 'Shape of `new` = ', shape(new)
        error stop
    end if
    n = size(old, dim=2)
    if (n == 1) then
        ! Angle between old and new
        cos_angle = dot_product(old(:,1), new(:,1))
        call cross(old(:,1), new(:,1), axis)
        sin_angle = norm(axis) !Always non-negative
        ! angle in [0, pi] as sin_angle > 0.
        angle = atan2(sin_angle, cos_angle) 
        ! try normalizing the axis
        call normalized(axis, estat)

        if ( estat == 0 ) then
            !Non-collinear vectors: standard case
            call aa_rotate_vectors(v, axis, angle, valigned)
        else
            !Collinear vectors: parallel or anti-parallel
            if ( cos_angle > 0.0_rp ) then
                !angle = 0; no rotation necessary
                valigned = v
            else
                !angle = pi: Rotation axis not unique. Pick a vector by
                !crossing old(:,1) with a cartesian direction least parallel
                !to it.
                !1. Find cartesian direction least parallel to old(:,1)
                min_idx = minloc(abs(old(:,1)), 1)
                !2. Create that basis vector
                vtmp = 0.0_rp
                vtmp(min_idx) = 1.0_rp
                !3. Cross old(:,1) with vtmp and unitize to get the rotation
                ! axis
                call cross(old(:,1), vtmp, axis) !axis will always have +ve norm.
                call aa_fix(axis, angle, .true., estat)
                !Rotate
                call aa_rotate_vectors(v, axis, angle, valigned)
            end if
        end if
    else if (n == 2) then
        !The third vector orthogonal to `old`
        call cross( old(:,1), old(:,2), z_old )
        call normalized(z_old, estat)
        !The third vector orthogonal to `new`
        call cross( new(:,1), new(:,2), z_new )
        call normalized(z_new, estat)
        !The `old` and new `frames`
        axes_old(:,1:2) = old; axes_old(:,3) = z_old
        axes_new(:,1:2) = new; axes_new(:,3) = z_new
        mat = matmul(axes_new, transpose(axes_old))
        valigned = matmul(mat, v)
    else if (n == 3) then
        mat = matmul(new, transpose(old))
        valigned = matmul(mat, v)
    end if

    end subroutine

!************************************************************************

logical function mat_is_rotmat(mat, tol)
    !! Returns `.true.` if a `(3,3)` matrix is a rotation matrix, `.false.`
    !! otherwise.

    real(rp), dimension(3,3), intent(in) :: mat
        !! `(3,3)`. Matrix to check.
    real(rp), intent(in) :: tol
        !! Absolute tolerance for comparing determinant and orthogonality.
    logical :: det_is_one, is_orthogonal
    real(rp) :: eye(3,3), mat_prod(3,3), det_mat

    call identity(eye)
    mat_prod = matmul(mat, transpose(mat))

    det_mat = det(mat)
    det_is_one = abs(det_mat-1.0_rp) <= tol
    is_orthogonal = maxval(abs(mat_prod - eye)) <= tol
    mat_is_rotmat = (is_orthogonal .and. det_is_one)

    end function

!******************************************************************************!
!                                  QUATERNION                                  !
!******************************************************************************!

subroutine quat_canonicalized(q)
    !! Modifies a non-canonical unit quaternion in-place to bring it to the
    !! canonical form.
    !!
    !! No modifications are performed if the input is already in canonical form.
    real(rp), dimension(4), intent(in out) :: q
        !! On entry, a unit quaternion. If `q` is not in canonical form,
        !! on return, `q` is modified to bring it to the canonical form, 
        !! else `q` is left as is.
    integer :: i

    do i = 1, 4
        if (q(i) > 0.0_rp) then
            return
        else if (q(i) < 0.0_rp) then
            q(1) = -q(1); q(2) = -q(2); q(3) = -q(3); q(4) = -q(4)
            return
        end if
    end do

    end subroutine

!******************************************************************************!

subroutine quat_rand(q)
    !! Draws a unit quaternion from a uniform distribution on the space of
    !! rotations in 3D.
    real(rp), dimension(4), intent(out) :: q
        !! Unit quaternion in canonical form.
    real(rp) :: u(3), sqrt_u1, sqrt_1mu1, cos_u2, sin_u2, cos_u3, sin_u3

    call random_number(u)

    sqrt_u1 = sqrt(u(1))
    sqrt_1mu1 = sqrt(1.0 - u(1))

    cos_u2 = cos(2.0_rp*math_pi*u(2))
    sin_u2 = sin(2.0_rp*math_pi*u(2))
    cos_u3 = cos(2.0_rp*math_pi*u(3))
    sin_u3 = sin(2.0_rp*math_pi*u(3))

    q(1) = sqrt_1mu1*cos_u2
    q(2) = sqrt_1mu1*sin_u2
    q(3) = sqrt_u1*sin_u3
    q(4) = sqrt_u1*cos_u3

    call quat_canonicalized(q)

    end subroutine 

!************************************************************************

subroutine quat_identity(q)
    !! Returns the identity unit quaternion.
    real(rp), dimension(4), intent(out) :: q
        !! Unit quaternion.
    q(1) = 1.0_rp; q(2) = 0.0_rp; q(3) = 0.0_rp; q(4) = 0.0_rp

    end subroutine 

!************************************************************************

subroutine quat_conjugated(q)
    !! Conjugates a quaternion in-place.
    real(rp), dimension(4), intent(in out) :: q
        !! On entry, a quaternion. On return, the conjugated quaternion.

    q(2:4) = -q(2:4)

    end subroutine

!************************************************************************

subroutine quat_inverted(q, estat)
    !! Inverts a quaternion in-place.
    real(rp), dimension(4), intent(in out) :: q
        !! On entry, a quaternion. On return, the inverted quaternion.
    integer, intent(out), optional :: estat
        !! Error code. Equals `1` for a zero-quaternion, else `0`.
    real(rp) :: d, q_scaled(4), s

    if (present(estat)) estat = 0

    s = maxval(abs(q))
    if (s == 0.0_rp) then
        if (present(estat)) estat = 1
        return
    else
        q_scaled = q/s
        d = dot_product(q_scaled, q_scaled)
        call quat_conjugated(q_scaled)
        q = (q_scaled/d)/s
    end if

    end subroutine

!************************************************************************

subroutine quat_normalized(q, canonicalize, estat)
    !! Normalizes a quaternion in-place. In case of a zero quaternion,
    !! `q` is not changed and an error code of `1` is returned.
    real(rp), dimension(4), intent(in out) :: q
        !! On entry, a quaternion. On return, the normalized quaternion
        !! (in canonical form if `canonicalize` is `.true.`).
    logical, intent(in), optional :: canonicalize
        !! Whether to canonicalize the output.
    integer, intent(out), optional :: estat
        !! Error code. Returns 1 for a zero quaternion, else 0.
    integer :: estat_

    call normalized(q, estat_)
    if (estat_ == 1) then
        if (present(estat)) estat = estat_
        return
    end if

    if (present(canonicalize)) then
        if (canonicalize) then
            call quat_canonicalized(q)
        end if
    end if

    if (present(estat)) estat = 0

    end subroutine

!************************************************************************

logical function quat_is_normalized(q, tol)
    !! Returns `.true.` if a quaternion is normalized, i.e. if it is a unit
    !! quaternion, within a tolerance `tol`; `.false.` otherwise.
    real(rp), dimension(4), intent(in) :: q
        !! A quaternion
    real(rp), intent(in) :: tol
        !! Absolute allowable deviation of the norm from one.
    real(rp) :: nrm

    nrm = norm(q)
    quat_is_normalized = abs(nrm - 1.0_rp) <= tol

    end function

!************************************************************************

subroutine quat_prod (p, q, pq)
    !! Calculates the product of two quaternions following Hamilton's
    !! convention. In case of unit quaternions describing rotations in a
    !! common frame, `q`'s rotation is applied first, followed by that of `p`.
    real(rp), dimension(4), intent(in)  :: p
        !! Quaternion.
    real(rp), dimension(4), intent(in)  :: q
        !! Quaternion.
    real(rp), dimension(4), intent(out) :: pq
        !! Product of `p` and `q`.

    pq(1) = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)
    pq(2) = p(2)*q(1) + p(1)*q(2) - p(4)*q(3) + p(3)*q(4)
    pq(3) = p(3)*q(1) + p(4)*q(2) + p(1)*q(3) - p(2)*q(4)
    pq(4) = p(4)*q(1) - p(3)*q(2) + p(2)*q(3) + p(1)*q(4)

    end subroutine

!************************************************************************

function quat_angle_between(p, q) result(angle)
    !! Returns the smallest rotation angle between the orientations represented
    !! by two unit quaternions.
    real(rp), dimension(4), intent(in) :: p
        !! Unit quaternion.
    real(rp), dimension(4), intent(in) :: q
        !! Unit quaternion.
    real(rp) :: angle
        !! Relative rotation angle in radians, in the interval `[0, pi]`. 
        !! Opposite-sign quaternions representing the same orientation have
        !! zero separation
    real(rp) :: q_(4), cos_angle, d_mag, s_mag

    q_ = q
    cos_angle = dot_product(p, q_)
    if (cos_angle < 0.0_rp) then
        ! Negate q_ to ensure both p and q_ lie on the same hemisphere
        q_ = -q_
    end if

    d_mag = norm(p - q_)
    s_mag = norm(p + q_)
    angle = 4.0_rp*atan2(d_mag, s_mag)

    end function

!************************************************************************

subroutine quat_interpolate(q1, q2, t, q)
    !! Performs spherical linear interpolation between two unit quaternions
    !! following the shortest path.
    !!
    !! The second endpoint may be negated to select the shortest path.
    !!
    !! At a relative orientation (physical) of exactly `pi`, two equally short
    !! paths exist; the supplied endpoint signs determine which is selected.
    real(rp), dimension(4), intent(in) :: q1
        !! Unit quaternion.
    real(rp), dimension(4), intent(in) :: q2
        !! Unit quaternion.
    real(rp),               intent(in)  :: t
        !! A fraction between 0 and 1 (both inclusive) specifying the
        !! interpolation point.
    real(rp), dimension(4), intent(out) :: q
        !! Interpolated unit quaternion.
    real(rp) :: q2_(4), cos_angle, d_mag, s_mag, theta, a, b, c

    q2_ = q2
    cos_angle = dot_product(q1, q2_)
    if (cos_angle < 0.0_rp) then
        ! Negate q2_ to ensure both q1 and q2_ lie on the same hemisphere
        q2_ = -q2_
    end if
    d_mag = norm(q1 - q2_)
    s_mag = norm(q1 + q2_)
    theta = 2.0_rp*atan2(d_mag, s_mag)

    ! SLERP
    c = sinc(theta)
    a = (1.0_rp - t) * sinc( (1.0_rp-t)*theta ) / c
    b = t * sinc(t*theta) / c
    q = a*q1 + b*q2_

    call quat_normalized(q, canonicalize=.false.)

    end subroutine

!************************************************************************

subroutine quat_deriv_to_angvel_mat(q, world, mat)
    !! Returns the matrix mapping the derivative of a unit quaternion to
    !! angular velocity.
    !!
    !! Let A (world frame) and B (body frame) be reference frames such that an
    !! active rotation in A by a unit quaternion `q` aligns its basis vectors
    !! with those of frame B.  Then the matrix M returned is such that `omega_i =
    !! M_ij * qdot_j`, where `omega_i` is the angular velocity of B relative to
    !! A resolved in A (or B) and `qdot_j` is the time derivative of `q`.
    !!
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    logical, intent(in) :: world
        !! If `.true.`, the angular velocity is resolved in frame A (world
        !! frame), else, it is resolved in frame B (body frame). 
    real(rp), dimension(3,4), intent(out) :: mat
        !! Angular velocity matrix.

    if (world) then
        mat = reshape([-q(2), -q(3), -q(4), &
                        q(1),  q(4), -q(3), &
                       -q(4),  q(1),  q(2), &
                        q(3), -q(2),  q(1) ], [3,4])
    else
        mat = reshape([-q(2), -q(3), -q(4), &
                        q(1), -q(4),  q(3), &
                        q(4),  q(1), -q(2), &
                       -q(3),  q(2),  q(1) ], [3,4])
    end if

    mat = 2.0_rp*mat

    end subroutine

!************************************************************************

subroutine quat_deriv_to_angvel (q, qdot, world, angvel)
    !! Calculates the angular velocity from a unit quaternion and its time
    !! derivative.
    !!
    !! Let A (world frame) and B (body frame) be reference frames such that an
    !! active rotation in A by a unit quaternion `q` aligns its basis vectors
    !! with those of frame B. Let `qdot` be the time derivative of `q`. The
    !! quantity returned is the angular velocity of B relative to A resolved in A
    !! (or B).
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(4), intent(in)  :: qdot
        !! Time derivative of `q`. In general, this is not a unit vector.
        !! Note that `qdot` is orthogonal to `q`, i.e. `q . qdot = 0`.
    logical, intent(in) :: world
        !! If `.true.`, the angular velocity is resolved in frame A (world frame),
        !! else, it is resolved in frame B (body frame).
    real(rp), dimension(3), intent(out) :: angvel
        !! Angular velocity.
    real(rp) :: mat(3,4)

    call quat_deriv_to_angvel_mat(q, world, mat)
    angvel = matmul(mat, qdot)

    end subroutine

!************************************************************************

subroutine quat_deriv_from_angvel_mat (q, world, mat)
    !! Returns the matrix mapping angular velocity to time derivative of a 
    !! unit quaternion.
    !!
    !! Let A (world frame) and B (body frame) be reference frames such that an
    !! active rotation in A by a unit quaternion `q` aligns its basis vectors
    !! with those of frame B. Then the matrix M returned is such that `qdot_i =
    !! M_ij * omega_j`, where `omega_j` is the angular velocity of B relative
    !! to A resolved in A (or B) and `qdot_i` is the time derivative of `q`.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    logical, intent(in) :: world
        !! If `.true.`, the angular velocity is resolved in frame A (world
        !! frame), else, it is resolved in frame B (body frame).
    real(rp), dimension(4,3), intent(out) :: mat
        !! Quaternion derivative matrix.
    
    if (world) then
        mat = reshape([-q(2),  q(1), -q(4),  q(3), &
                       -q(3),  q(4),  q(1), -q(2), &
                       -q(4), -q(3),  q(2),  q(1)], [4,3])
    else
        mat = reshape([-q(2),  q(1),  q(4), -q(3), &
                       -q(3), -q(4),  q(1),  q(2), &
                       -q(4),  q(3), -q(2),  q(1)], [4,3])
    end if

    mat = 0.5_rp*mat

    end subroutine

!************************************************************************

subroutine quat_deriv_from_angvel(q, angvel, world, qdot)
    !! Calculates the time derivative of a unit quaternion from angular velocity.
    !!
    !! Let A (world frame) and B (body frame) be reference frames such that an
    !! active rotation in A by a unit quaternion `q` aligns its basis vectors
    !! with those of frame B. Let `angvel` be the angular velocity of B
    !! relative to A resolved in A (or B). The quantity returned is the time
    !! derivative of `q`.
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion. An active rotation of frame A by `q` aligns it with
        !! frame B.
    real(rp), dimension(3), intent(in)  :: angvel
        !! Angular velocity of frame B relative to frame A resolved in A if
        !! `world=.true.`, or in B otherwise.
    logical, intent(in) :: world
        !! If `.true.`, the angular velocity is resolved in frame A (world frame),
        !! else, it is resolved in frame B (body frame).
    real(rp), dimension(4), intent(out) :: qdot
        !! Time derivative of `q`. It is orthogonal to `q` and not a unit vector in
        !! general.
    real(rp) :: mat(4,3)

    call quat_deriv_from_angvel_mat(q, world, mat)
    qdot = matmul(mat, angvel)

    end subroutine

!************************************************************************

subroutine quat_rotmat (q, rotmat)
    !! Returns the rotation matrix corresponding to a unit quaternion.
    !!
    !! For a column vector `v` of length 3, the rotated vector `v' = R v`, where
    !! `R` is the returned matrix.  The components of `v` and `v'` are
    !! expressed in the same fixed orthonormal basis.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(3,3), intent(out) :: rotmat
        !! Rotation matrix.
    real(rp) :: q1sq, q2sq, q3sq, q4sq
    real(rp) :: q1q2, q1q3, q1q4, q2q3, q2q4, q3q4

    q1sq = q(1)*q(1)
    q2sq = q(2)*q(2)
    q3sq = q(3)*q(3)
    q4sq = q(4)*q(4)

    q1q2 = q(1)*q(2)
    q1q3 = q(1)*q(3)
    q1q4 = q(1)*q(4)
    q2q3 = q(2)*q(3)
    q2q4 = q(2)*q(4)
    q3q4 = q(3)*q(4)

    rotmat(1,1) = 2.0_rp*(q1sq + q2sq) - 1.0_rp
    rotmat(2,1) = 2.0_rp*(q2q3 + q1q4)
    rotmat(3,1) = 2.0_rp*(q2q4 - q1q3)

    rotmat(1,2) = 2.0_rp*(q2q3 - q1q4)
    rotmat(2,2) = 2.0_rp*(q1sq + q3sq) - 1.0_rp
    rotmat(3,2) = 2.0_rp*(q3q4 + q1q2)

    rotmat(1,3) = 2.0_rp*(q2q4 + q1q3)
    rotmat(2,3) = 2.0_rp*(q3q4 - q1q2)
    rotmat(3,3) = 2.0_rp*(q1sq + q4sq) - 1.0_rp

    end subroutine

!************************************************************************

subroutine quat_shiftmat(q, forward, shiftmat)
    !! Returns the shifter matrix corresponding to a unit quaternion.
    !!
    !! The orientation of frame B relative to frame A is specified by a unit
    !! quaternion `q`, such that an active rotation in A of its basis vectors
    !! aligns them with those of B. The shifter matrix `S` transforms the
    !! components of a fixed vector resolved in A to those resolved in B or
    !! vice-versa according to `v_B_i = S_ij * v_A_j` or `v_A_i = S_ij *
    !! v_B_j`.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: shiftmat
        !! Shifter matrix.
    real(rp), dimension(4) :: conj_q
    
    if (forward) then
        conj_q = q
        call quat_conjugated(conj_q)
        call quat_rotmat(conj_q, shiftmat)
    else
        call quat_rotmat(q, shiftmat)
    end if

    end subroutine

!************************************************************************

subroutine quat_rotate_vectors(v, q, vrot)
    !! Rotates vectors by a unit quaternion. The components of `v` and the 
    !! rotated vectors are all expressed in the same fixed basis.
    real(rp), dimension(:,:), intent(in)  :: v
        !! `(3,n)`. Vectors to rotate, stored as columns.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! `(3,n)` Rotated vectors.
    real(rp), dimension(3,3) :: rotmat

    call quat_rotmat(q, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine quat_shift_vectors(v, q, forward, vshift)
    !! Transform vector components between frames A and B.
    !!  
    !! The orientation of B relative to A is specified by a unit quaternion `q`.
    !! This orientation is used for both transformation directions.  The vectors
    !! remain unchanged; only their components are expressed in a different frame.

    real(rp), dimension(:,:), intent(in)  :: v
        !! `(3,n)`. Components of vectors stored as columns. Expressed in A if
        !! `forward` is `.true.`, or in B otherwise.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion expressed in A, regardless of `forward`.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`, 
        !! transform components from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! `(3,n)` Components of the same vectors expressed in B if `forward`
        !! is `.true.`, or in A otherwise. Has the same shape as `v`.
    real(rp), dimension(3,3) :: shiftmat

    call quat_shiftmat(q, forward, shiftmat)
    vshift = matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine quat_shift_tensor2(a, q, forward, ashift)
    !! Transform second-order tensor components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by a unit quaternion `q`.
    !! This orientation is used for both transformation directions.  The tensor
    !! remains unchanged; only its components are expressed in a different frame.

    real(rp), dimension(3,3), intent(in)  :: a
        !! Components of a second-order tensor.  Expressed in A if `forward` is
        !! `.true.`, or in B otherwise.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion expressed in A, regardless of `forward`.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`, 
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Components of the same second-order tensor expressed in B if
        !! `forward` is `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call quat_shiftmat(q, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine quat_shift_tensor3(a, q, forward, ashift)
    !! Transform third-order tensor components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by a unit quaternion `q`.
    !! This orientation is used for both transformation directions.  The tensor
    !! remains unchanged; only its components are expressed in a different frame.

    real(rp), dimension(3,3,3), intent(in)  :: a
        !! Components of a third-order tensor.  Expressed in A if `forward` is
        !! `.true.`, or in B otherwise.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion expressed in A, regardless of `forward`.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`, 
        !! transform components from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Components of the same third-order tensor expressed in B if
        !! `forward` is `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call quat_shiftmat(q, forward, shiftmat)
    call shift_tensor3(a, shiftmat, ashift)

    end subroutine 

!************************************************************************

subroutine quat_to_aa (q, axis, angle)
    !! Converts a unit quaternion to an *axis-angle* representation. An
    !! identity quaternion will return [1,0,0] as the axis.

    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(3), intent(out) :: axis
        !! Unit vector along the axis of rotation.
    real(rp),               intent(out) :: angle
        !! Angle in radian (lies in `[0, pi]`).
    real(rp) :: p(3), q_(4), r, s

    q_ = q
    if (q_(1) < 0.0_rp) q_ = -q_
    s = maxval( abs( q_(2:4) ) )

    if (s == 0.0_rp) then
        ! Identity rotation (angle = 0)
        axis(1) = 1.0_rp; axis(2) = 0.0_rp; axis(3) = 0.0_rp
        angle = 0.0_rp
    else
        ! Usual case
        p = q_(2:4)/s
        r = norm(p)
        axis = p/r
        angle = 2.0_rp*atan2(s*r, q_(1))
    end if

    end subroutine

!************************************************************************

subroutine quat_to_dcm (q, dcm)
    !! Converts a unit quaternion to a direction cosine matrix.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(3,3), intent(out) :: dcm
        !! Direction cosine matrix.

    call quat_shiftmat(q, .true., dcm)

    end subroutine

!************************************************************************

subroutine quat_to_euler (q, seq, world, tol, euler)
    !! Converts a unit quaternion to an Euler angle sequence.
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    character(3),           intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp),               intent(in) :: tol
        !! Tolerance (must be >= 0) for setting one outer angle to zero.
    real(rp), dimension(3), intent(out) :: euler
        !! Euler angles.
    real(rp), dimension(3,3) :: rotmat

    call quat_rotmat(q, rotmat)
    call factor_rotmat(rotmat, seq, world, tol, euler)

    end subroutine

!*****************************************************************************!
!                                  AXIS ANGLE                                 !
!*****************************************************************************!

subroutine aa_fix (axis, angle, normalize_axis, estat)
    !! For a given orientation in *axis-angle* representation, returns a
    !! modified `axis` and `angle` to ensure that `angle` lies in `[0, pi]`
    !! while preserving the rotation.
    real(rp), dimension(3), intent(in out) :: axis
        !! Axis of rotation. On entry, must be a non-zero vector. On return,
        !! if `normalize_axis=.true.`, the normalized `axis`. 
        !! If `normalize_axis=.false.`, it has the same magnitude as `axis`
        !! and may need to be normalized before using other axis-angle functions.
    real(rp), intent(in out) :: angle
        !! On entry, angle in radian. On return, modified angle in radian
        !! within `[0, pi]`.
    logical, intent(in) :: normalize_axis
        !! Whether to normalize the axis to a unit vector.
    integer, intent(out) :: estat
        !! Error code. Equals `1` if `axis` is a zero vector, else `0`.

    estat = 0

    ! Check if axis is a zero vector
    if ((axis(1)==0.0_rp) .and. (axis(2)==0.0_rp) .and. (axis(3)==0.0_rp)) then
        estat = 1; return 
    end if

    if (normalize_axis) call normalized(axis, estat)
    if (estat == 1) return
    ! Get principal signed angle in [-pi, pi]
    angle = atan2(sin(angle), cos(angle))
    ! Make the angle non-negative
    if (angle < 0.0_rp) then
        angle = -angle
        axis = -axis
    end if

    end subroutine

!************************************************************************

subroutine aa_rand (axis, angle)
    !! Generates a random orientation in *axis-angle* representation drawn from
    !! a uniform distribution on the space of rotations in 3D.
    real(rp), dimension(3), intent(out) :: axis
    real(rp), intent(out) :: angle
    real(rp) :: q(4)

    call quat_rand(q)
    call quat_to_aa(q, axis, angle)

    end subroutine

!************************************************************************

subroutine aa_rotmat (axis, angle, rotmat)
    !! Returns the matrix for an active rotation about `axis` through `angle`,
    !! with positive angles following the right-hand rule.
    !!
    !! For a column vector `v` of length 3, the rotated vector `v' = R v`, where
    !! `R` is the returned matrix. The components of `v`, `v'`, and `axis` are
    !! expressed in the same fixed orthonormal basis.

    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. Positive angles follow the right-hand
        !! rule.
    real(rp), dimension(3,3), intent(out) :: rotmat
        !! Rotation matrix.
    real(rp) :: sine, sine_2, cosine, icos

    sine = sin(angle)
    cosine = cos(angle)
    sine_2 = sin(angle/2.0_rp)
    icos = 2.0_rp*sine_2*sine_2 ! equals (1 - cosine)

    rotmat(1,1) = axis(1)*axis(1)*icos + cosine
    rotmat(2,1) = axis(1)*axis(2)*icos + axis(3)*sine
    rotmat(3,1) = axis(3)*axis(1)*icos - axis(2)*sine

    rotmat(1,2) = axis(1)*axis(2)*icos - axis(3)*sine
    rotmat(2,2) = axis(2)*axis(2)*icos + cosine
    rotmat(3,2) = axis(2)*axis(3)*icos + axis(1)*sine

    rotmat(1,3) = axis(1)*axis(3)*icos + axis(2)*sine
    rotmat(2,3) = axis(2)*axis(3)*icos - axis(1)*sine
    rotmat(3,3) = axis(3)*axis(3)*icos + cosine

    end subroutine

!************************************************************************

subroutine aa_from_rotmat (rotmat, axis, angle)
    !! Extracts axis and angle from a rotation matrix.
    !! If `rotmat` is an identity matrix, the axis returned is [1,0,0]. If the
    !! angle is exactly `pi`, either sign of the axis is valid.
    !! 
    real(rp), dimension(3,3), intent(in) :: rotmat
        !! Rotation matrix (must be orthonormal with determinant equal to one).
    real(rp), dimension(3), intent(out) :: axis
        !! Unit vector along the axis.
    real(rp), intent(out) :: angle
        !! Angle of rotation in radian within `[0, pi]`.
    real(rp) :: cos_angle, sin_angle, s, tr, w(3)
    integer  :: k, estat

    ! Cos of the angle from the trace
    tr = rotmat(1,1) + rotmat(2,2) + rotmat(3,3)
    cos_angle = (tr-1.0_rp)/2.0_rp
    ! Clamp cos_angle to [-1, 1]
    cos_angle = max(-1.0_rp, min(1.0_rp, cos_angle))

    ! Sin of the angle from the anti-symmetric vector w
    w(1) = 0.5_rp * (rotmat(3,2) - rotmat(2,3))
    w(2) = 0.5_rp * (rotmat(1,3) - rotmat(3,1))
    w(3) = 0.5_rp * (rotmat(2,1) - rotmat(1,2))
    sin_angle = norm(w)

    ! Calculate the angle
    angle = atan2(sin_angle, cos_angle)

    if (cos_angle >= 0.0_rp) then
        ! Angle lies in [0, pi/2]
        if (sin_angle > 0.0_rp) then
            axis = w/sin_angle
        else
            ! Angle = 0; axis non-unique
            angle = 0.0_rp
            axis = [1.0_rp, 0.0_rp, 0.0_rp]
        end if
    else
        ! Find the largest entry on the diagonal of rotmat
        k = maxloc( [rotmat(1,1), rotmat(2,2), rotmat(3,3)], 1)
        select case (k)
        case (1)
            axis(1) = sqrt( (rotmat(1,1)-cos_angle) / (1.0_rp-cos_angle) )
            s = 1.0_rp/(2.0_rp*(1.0_rp-cos_angle)*axis(1))
            axis(2) = s*(rotmat(2,1)+rotmat(1,2))
            axis(3) = s*(rotmat(3,1)+rotmat(1,3))
        case (2)
            axis(2) = sqrt( (rotmat(2,2)-cos_angle) / (1.0_rp-cos_angle) )
            s = 1.0_rp/(2.0_rp*(1.0_rp-cos_angle)*axis(2))
            axis(1) = s*(rotmat(1,2)+rotmat(2,1))
            axis(3) = s*(rotmat(3,2)+rotmat(2,3))
        case (3)
            axis(3) = sqrt( (rotmat(3,3)-cos_angle) / (1.0_rp-cos_angle) )
            s = 1.0_rp/(2.0_rp*(1.0_rp-cos_angle)*axis(3))
            axis(1) = s*(rotmat(1,3)+rotmat(3,1))
            axis(2) = s*(rotmat(2,3)+rotmat(3,2))
        end select
    end if

    ! Resolve axis sign
    if (dot_product(axis, w) < 0.0_rp) axis = -axis

    call aa_fix(axis, angle, .true., estat)

    end subroutine

!************************************************************************

subroutine aa_rotate_vectors (v, axis, angle, vrot)
    !! Rotates vectors about `axis` by `angle`.
    !!
    !! The components of `axis`, `v`, and the rotated vectors are all
    !! expressed in the same fixed basis.
    real(rp), dimension(:,:), intent(in) :: v
        !! `(3,n)`. Vectors to rotate stored as columns.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. Positive angles follow the right-hand
        !! rule.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! `(3,n)`. Rotated vectors.
    real(rp), dimension(3,3) :: rotmat

    call aa_rotmat(axis, angle, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine aa_shiftmat (axis, angle, forward, shiftmat)
    !! Returns the shifter matrix corresponding to an *axis-angle*
    !! representation.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. Positive angles follow the right-hand
        !! rule.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`, 
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: shiftmat
        !! Shifter matrix.

    !Get shiftmat assuming forward == .true.
    call aa_rotmat(-axis, angle, shiftmat) 

    if (.not. forward) then
        shiftmat = transpose(shiftmat)
    end if

    end subroutine

!************************************************************************

subroutine aa_shift_vectors (v, axis, angle, forward, vshift)
    !! Transform vector components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by `axis` and `angle`,
    !! with `axis` expressed in A and positive angles following the right-hand
    !! rule. This orientation is used for both transformation directions.
    !! The vectors remain unchanged; only their components are expressed in a
    !! different frame.
    real(rp), dimension(:,:), intent(in) :: v
        !! `(3,n)`. Components of a vectors stored as columns. Expressed in
        !! A if `forward` is `.true.`, or in B otherwise.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis expressed in A, regardless of `forward`
    real(rp), intent(in) :: angle
        !! Angle specifying the orientation of B relative to A, in radians.
        !! Positive angles follow the right-hand rule.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! `(3,n)`. Components of the same vectors expressed in B if 
        !! `forward` is `.true.`, or in A otherwise. Has the same shape as `v`.
    real(rp), dimension(3,3) :: shiftmat

    call aa_shiftmat(axis, angle, forward, shiftmat)
    vshift = matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine aa_shift_tensor2 (a, axis, angle, forward, ashift)
    !! Transforms second-order tensor components between frames A and B.
    !!  
    !! The orientation of B relative to A is specified by `axis` and `angle`,
    !! with `axis` expressed in A and positive angles following the right-hand
    !! rule. This orientation is used for both transformation directions.
    !! The tensor remains unchanged; only its components are expressed in a
    !! different frame.
    real(rp), dimension(3,3), intent(in) :: a
        !! Components of a second-order tensor, expressed in A if `forward` is
        !! `.true.`, or in B otherwise.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis expressed in A, regardless of `forward`
    real(rp), intent(in) :: angle
        !! Angle specifying the orientation of B relative to A, in radians.
        !! Positive angles follow the right-hand rule.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Components of the second-order tensor expressed in B if `forward` is
        !! `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call aa_shiftmat(axis, angle, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine aa_shift_tensor3 (a, axis, angle, forward, ashift)
    !! Transforms third-order tensor components between frames A and B.
    !!  
    !! The orientation of B relative to A is specified by `axis` and `angle`,
    !! with `axis` expressed in A and positive angles following the right-hand
    !! rule. This orientation is used for both transformation directions.
    !! The tensor remains unchanged; only its components are expressed in a
    !! different frame.
    real(rp), dimension(3,3,3), intent(in) :: a
        !! Components of a third-order tensor, expressed in A if `forward` is
        !! `.true.`, or in B otherwise.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis expressed in A, regardless of `forward`
    real(rp), intent(in) :: angle
        !! Angle specifying the orientation of B relative to A, in radians.
        !! Positive angles follow the right-hand rule.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Components of the third-order tensor expressed in B if `forward` is
        !! `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call aa_shiftmat(axis, angle, forward, shiftmat)
    call shift_tensor3(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine aa_to_quat(axis, angle, q)
    !! Converts an *axis-angle* representation to a unit quaternion.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian.
    real(rp), dimension(4), intent(out) :: q
        !! Unit quaternion in the canonical form.

    q(1) = cos(0.5_rp*angle)
    q(2:4) = sin(0.5_rp*angle)*axis

    !Normalize to ensure unit quaternion and canonicalize 
    call quat_normalized(q, .true.)

    end subroutine

!************************************************************************

subroutine aa_to_dcm (axis, angle, dcm)
    !! Converts an *axis-angle* representation to a direction cosine matrix.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along  the axis.
    real(rp), intent(in)   :: angle
        !! Angle of rotation in radian.
    real(rp), dimension(3,3), intent(out) :: dcm
        !! Direction cosine matrix.

    call aa_shiftmat(axis, angle, .true., dcm)

    end subroutine

!************************************************************************

subroutine aa_to_euler (axis, angle, seq, world, tol, euler)
    !! Coverts an *axis-angle* representation to *Euler angles*.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle in radian.
    character(len=3), intent(in) :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), intent(in) :: tol
        !! Tolerance (must be >= 0) for setting one outer angle to zero.
    real(rp), dimension(3), intent(out) :: euler
        !! Euler angles.
    real(rp), dimension(3,3) :: rotmat
    
    call aa_rotmat(axis, angle, rotmat)
    call factor_rotmat(rotmat, seq, world, tol, euler)

    end subroutine

!*****************************************************************************!
!                           DIRECTION COSINE MATRIX                           !
!*****************************************************************************!

logical function mat_is_dcm (mat, tol)
    !! Checks if `mat` is a direction cosine matrix or not.
    real(rp), dimension(3,3), intent(in) :: mat
        !! Matrix to check.
    real(rp), intent(in) :: tol
        !! Absolute tolerance for comparing determinant and orthogonality.

    mat_is_dcm = mat_is_rotmat(mat, tol)

    end function

!************************************************************************

subroutine dcm_from_frames (A, B, dcm)
    !! Returns the direction cosine matrix of frame B with respect to frame A,
    !! where the bases of the frames are resolved in an underlying fixed frame F.
    !! If F coincides with A, then `A` reduces to an identity matrix.
    real(rp), dimension(3,3), intent(in)  :: A
       !! The columns of A represent the orthonormal basis vectors of a
       !! right-handed frame A resolved in F.
    real(rp), dimension(3,3), intent(in)  :: B
       !! The columns of B represent the orthonormal basis vectors of a
       !! right-handed frame B resolved in F.
    real(rp), dimension(3,3), intent(out) :: dcm
        !! Direction cosine matrix of frame B with respect to frame A

    dcm =  matmul(transpose(B), A)

    end subroutine

!************************************************************************

subroutine dcm_rotmat (dcm, rotmat)
    !! Returns the rotation matrix corresponding to a direction cosine matrix.
    !!
    !! Let D be the direction cosine matrix of frame B relative to frame A.  The
    !! rotation matrix R corresponding to D rotates vectors resolved in A such
    !! that the bases of A align with those of B, when all the bases are resolved
    !! in A.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(3,3), intent(out) :: rotmat
        !! Rotation matrix.

    rotmat = transpose(dcm)

    end subroutine

!************************************************************************

subroutine dcm_shiftmat (dcm, forward, shiftmat)
    !! Returns the shifter matrix corresponding to a direction cosine matrix.
    !!
    !! Let D be the direction cosine matrix of frame B relative to frame A. The
    !! shifter matrix S corresponding to D transforms components of vectors (or
    !! tensors) resolved in A to those in B (or vice-versa).
    real(rp), dimension(3,3), intent(in)  :: dcm
        !! Direction cosine matrix.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: shiftmat
        !! Shifter matrix.

    if (forward) then
        shiftmat = dcm
    else
        shiftmat = transpose(dcm)
    end if

    end subroutine

!************************************************************************

subroutine dcm_rotate_vectors (v, dcm, vrot)
    !! Rotates vectors by a direction cosine matrix.
    !!
    !! Let D be the direction cosine matrix of frame B relative to frame A.
    !! Vectors resolved in A are rotated such that the bases of A align with those
    !! of B, when all the bases are resolved in A.
    real(rp), dimension(:,:), intent(in) :: v
        !! `(3,n)`. Components of a single vector or n vectors stored as columns.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! `(3,n)`. Rotated vectors.
    real(rp) :: rotmat(3,3)

    call dcm_rotmat(dcm, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine dcm_shift_vectors (v, dcm, forward, vshift)
    !! Transform vector components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by `dcm`. This orientation
    !! is used for both transformation directions. The vectors remain unchanged;
    !! only their components are expressed in a different frame.
    real(rp), dimension(:,:), intent(in) :: v
        !! `(3,n)`. Components of vectors stored as columns. Expressed
        !! in A if `forward` is `.true.`, or in B otherwise.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix of frame B with respect to frame A.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! `(3,n)`. Components of the same vectors expressed in B if 
        !! `forward` is `.true.`, or in A otherwise. Has the same shape as `v`.
    real(rp), dimension(3,3) :: shiftmat

    call dcm_shiftmat(dcm, forward, shiftmat)
    vshift = matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine dcm_shift_tensor2 (a, dcm, forward, ashift)
    !! Transform second-order tensor components between frames A and B.
    !!
    !! The orientation of B relative to A is specified by `dcm`. This orientation
    !! is used for both transformation directions. The tensor remains unchanged;
    !! only its components are expressed in a different frame.
    real(rp), dimension(3,3), intent(in) :: a
        !! A second-order tensor, resolved in A if `forward` is `.true.`,
        !! or in B otherwise.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix of frame B with respect to frame A.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Components of the same tensor expressed in B if 
        !! `forward` is `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call dcm_shiftmat(dcm, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine dcm_shift_tensor3 (a, dcm, forward, ashift)
    !! Transform third-order tensor components between frames A and B.
    !!
    !! The orientation of B relative to A is specified by `dcm`. This orientation
    !! is used for both transformation directions. The tensor remains unchanged;
    !! only its components are expressed in a different frame.
    real(rp), dimension(3,3,3), intent(in) :: a
        !! A third-order tensor, resolved in A if `forward` is `.true.`,
        !! or in B otherwise.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix of frame B with respect to frame A.
    logical, intent(in) :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Components of the same tensor expressed in B if 
        !! `forward` is `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call dcm_shiftmat(dcm, forward, shiftmat)
    call shift_tensor3(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine dcm_to_quat (dcm, q)
    !! Converts a direction cosine matrix to a unit quaternion.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(4), intent(out)  :: q
        !! Unit quaternion.
    real(rp), dimension(3,3) :: rotmat
    real(rp), dimension(3) :: axis(3)
    real(rp) :: angle

    call dcm_rotmat(dcm, rotmat)
    call aa_from_rotmat(rotmat, axis, angle)
    call aa_to_quat(axis, angle, q)

    end subroutine

!************************************************************************

subroutine dcm_to_aa (dcm, axis, angle)
    !! Converts a direction cosine matrix to an *axis-angle* representation.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(3), intent(out)  :: axis
        !! Unit vector along the axis.
    real(rp), intent(out) :: angle
        !! Angle in radian (lies in `[0, pi]`).
    real(rp), dimension(3,3) :: rotmat

    call dcm_rotmat(dcm, rotmat)
    call aa_from_rotmat(rotmat, axis, angle)

    end subroutine

!************************************************************************

subroutine dcm_to_euler (dcm, seq, world, tol, euler)
    !! Converts a direction cosine matrix to Euler angles.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    character(3), intent(in) :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), intent(in) :: tol
        !! Tolerance (must be >= 0) for setting one outer angle to zero.
    real(rp), dimension(3), intent(out)  :: euler
        !! Euler angles
    real(rp), dimension(3,3) :: rotmat

    call dcm_rotmat (dcm, rotmat)
    call factor_rotmat(rotmat, seq, world, tol, euler)

    end subroutine

!*****************************************************************************!
!                                 EULER ANGLES                                !
!*****************************************************************************!

subroutine euler_rotmat (euler, seq, world, rotmat)
    !! Returns the rotation matrix for a set of Euler angles.
    real(rp), dimension(3), intent(in) :: euler
    character(len=3), intent(in) :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3,3), intent(out) :: rotmat
        !! Rotation matrix.
    real(rp), dimension(3) :: euler_

    if (.not. world) then 
        euler_ = -euler
    else
        euler_ = euler
    end if

    select case (seq)
    case ('XYZ')
        call rotmat_xyz(euler_, rotmat)
    case ('XZY')
        call rotmat_xzy(euler_, rotmat)
    case ('YXZ')
        call rotmat_yxz(euler_, rotmat)
    case ('YZX')
        call rotmat_yzx(euler_, rotmat)
    case ('ZXY')
        call rotmat_zxy(euler_, rotmat)
    case ('ZYX')
        call rotmat_zyx(euler_, rotmat)
    case default
        write(*,*) 'seq = ', seq
        error stop 'seq out of allowed values.'
    end select

    if (.not. world) rotmat = transpose(rotmat)

    end subroutine

!*******************************************************************************

subroutine euler_shiftmat (euler, seq, world, forward, shiftmat)
    !! Returns the shifter matrix for a set of Euler angles.
    !!
    !! The orientation of frame B relative to frame A is specified by a set of
    !! Euler angles. The shifter matrix S transforms components of vectors (or
    !! tensors) resolved in A to those in B (or vice-versa).

    real(rp), dimension(3), intent(in) :: euler
        !! Euler angles.
    character(3),           intent(in) :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: shiftmat
        !! Shifter matrix.
    real(rp), dimension(3,3) :: rotmat

    call euler_rotmat(euler, seq, world, rotmat)

    if (forward) then
        shiftmat = transpose(rotmat)
    else
        shiftmat = rotmat
    end if

    end subroutine

!************************************************************************

subroutine euler_rotate_vectors (v, euler, seq, world, vrot)
    !! Rotates vectors with a set of Euler angles. The components of `v`
    !! and the rotated vectors are all expressed in the same fixed basis.
    real(rp), dimension(:,:), intent(in)  :: v
        !! `(3,n)`. Components of vectors stored as columns.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! `(3,n)`. Rotated vectors.
    real(rp), dimension(3,3) :: rotmat

    call euler_rotmat (euler, seq, world, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine euler_shift_vectors (v, euler, seq, world, forward, vshift)
    !! Transform vector components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by Euler angles.
    !! This orientation is used for both transformation directions. The vectors
    !! remain unchanged; only their components are expressed in a different
    !! frame.
    real(rp), dimension(:,:), intent(in)  :: v
        !! `(3,n)`. Components of vectors stored as columns. Expressed
        !! in A if `forward` is `.true.`, or in B otherwise.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! `(3,n)`. Components of the same vectors expressed in B if 
        !! `forward` is `.true.`, or in A otherwise. Has the same shape as `v`.
    real(rp), dimension(3,3) :: shiftmat

    call euler_shiftmat (euler, seq, world, forward, shiftmat)
    vshift =  matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine euler_shift_tensor2 (a, euler, seq, world, forward, ashift)
    !! Transform second-order tensor components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by Euler angles.
    !! This orientation is used for both transformation directions. The vectors
    !! remain unchanged; only their components are expressed in a different
    !! frame.
    real(rp), dimension(3,3), intent(in)  :: a
        !! A second-order tensor, resolved in A if `forward` is `.true.`,
        !! or in B otherwise.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Components of the same tensor expressed in B if 
        !! `forward` is `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call euler_shiftmat (euler, seq, world, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine euler_shift_tensor3 (a, euler, seq, world, forward, ashift)
    !! Transform third-order tensor components between frames A and B.
    !! 
    !! The orientation of B relative to A is specified by Euler angles.
    !! This orientation is used for both transformation directions. The vectors
    !! remain unchanged; only their components are expressed in a different
    !! frame.
    real(rp), dimension(3,3,3), intent(in)  :: a
        !! A third-order tensor, resolved in A if `forward` is `.true.`,
        !! or in B otherwise.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, transform components from A to B. If `.false.`,
        !! transform components from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Components of the same tensor expressed in B if 
        !! `forward` is `.true.`, or in A otherwise.
    real(rp), dimension(3,3) :: shiftmat

    call euler_shiftmat (euler, seq, world, forward, shiftmat)
    call shift_tensor3(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine euler_to_quat (euler, seq, world, q)
    !! Convert Euler angles to a unit quaternion.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(4),     intent(out) :: q
        !! Unit quaternion in canonical form.
    real(rp), dimension(3) :: axis
    real(rp) :: angle

    call euler_to_aa(euler, seq, world, axis, angle)
    call aa_to_quat(axis, angle, q)

    end subroutine

!************************************************************************

subroutine euler_to_aa (euler, seq, world, axis, angle)
    !! Convert Euler angles to an *axis-angle* representation.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3),     intent(out) :: axis
        !! Unit vector along the axis of rotation.
    real(rp),                   intent(out) :: angle
        !! Angle of rotation in radian (lies in `[0, pi]`).
    real(rp), dimension(3,3) :: rotmat

    call euler_rotmat (euler, seq, world, rotmat)
    call aa_from_rotmat(rotmat, axis, angle)

    end subroutine

!************************************************************************

subroutine euler_to_dcm (euler, seq, world, dcm)
    !! Convert Euler angles to a direction cosine matrix.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3,3),   intent(out) :: dcm
        !! Direction cosine matrix.

    call euler_shiftmat (euler, seq, world, .true., dcm)

    end subroutine

!************************************************************************

subroutine euler_to_euler (euler, seq, world, to_seq, to_world, tol, to_euler)
    !! Convert one set of Euler angles to another.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3), intent(in)  :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Sequence for `euler`.
    logical, intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    character(3), intent(in) :: to_seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Convert `euler` with sequence `seq` to the sequence `to_seq`.
    logical, intent(in) :: to_world
        !! Whether the converted euler angles are with respect to the *world* frame or not.
    real(rp), intent(in) :: tol
        !! Tolerance (must be >= 0) for setting one outer angle to zero.
    real(rp), dimension(3),     intent(out) :: to_euler
        !! Converted Euler angles.
    real(rp), dimension(3,3) :: rotmat

    call euler_rotmat (euler, seq, world, rotmat)
    call factor_rotmat(rotmat, to_seq, to_world, tol, to_euler)

    end subroutine

!*****************************************************************************!
!                           INTERNAL ROUTINES                                 !
!*****************************************************************************!

subroutine factor_rotmat (rotmat, seq, world, tol, euler)
    !! Factorize a rotation matrix to obtain the three Euler angles.
    real(rp), dimension(3,3), intent(in) :: rotmat
        !! Rotation matrix.
    character(3), intent(in) :: seq
        !! `'XYZ'` | `'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), intent(in) :: tol
        !! Tolerance (must be >= 0) for setting one outer angle to zero.
    real(rp), dimension(3), intent(out) :: euler
        !! Euler angles.
    real(rp), dimension(3,3) :: rotmat_

    if (.not. world) then
        rotmat_ = transpose(rotmat)
    else
        rotmat_ = rotmat
    end if

    select case (seq)
    case ('XYZ')
        call factor_rotmat_xyz(rotmat_, tol, euler)
    case ('XZY')
        call factor_rotmat_xzy(rotmat_, tol, euler)
    case ('YXZ')
        call factor_rotmat_yxz(rotmat_, tol, euler)
    case ('YZX')
        call factor_rotmat_yzx(rotmat_, tol, euler)
    case ('ZXY')
        call factor_rotmat_zxy(rotmat_, tol, euler)
    case ('ZYX')
        call factor_rotmat_zyx(rotmat_, tol, euler)
    case default
        write(*,*) 'seq = ', seq
        error stop 'seq out of allowed values.'
    end select

    if (.not. world) euler = -euler

    end subroutine

!************************************************************************

subroutine shift_tensor2(a, shiftmat, ashift)
    !! Shifts a second order tensor, i.e. `(3,3)` matrix using a shifter matrix

    real(rp), dimension(3,3), intent(in) :: a
    real(rp), dimension(3,3), intent(in) :: shiftmat
    real(rp), dimension(3,3), intent(out) :: ashift

    ashift = matmul(shiftmat, matmul(a, transpose(shiftmat)))

    end subroutine

!************************************************************************

subroutine shift_tensor3(a, shiftmat, ashift)
    !! Shifts a third order tensor, i.e. `(3,3,3)` matrix using a shifter matrix

    real(rp), dimension(3,3,3), intent(in) :: a
    real(rp), dimension(3,3), intent(in) :: shiftmat
    real(rp), dimension(3,3,3), intent(out) :: ashift
    integer :: i, j, k, p, q, r

    ashift = 0.0_rp

    do r = 1, 3
        do q = 1, 3
            do p = 1, 3
                do k = 1, 3
                    do j = 1, 3
                        do i = 1, 3
                            ashift(i,j,k) = ashift(i,j,k) + &
                                shiftmat(i,p)*shiftmat(j,q)*shiftmat(k,r)*a(p,q,r)
                        end do
                    end do
                end do
            end do
        end do
    end do


    end subroutine

!*******************************************************************************

subroutine rotmat_xyz(euler, rotmat)

    real(rp), dimension(3), intent(in) :: euler
    real(rp), dimension(3,3), intent(out) :: rotmat
    real(rp) :: phi, theta, psi
    real(rp) :: sin_phi, sin_theta, sin_psi
    real(rp) :: cos_phi, cos_theta, cos_psi

    rotmat = 0.0_rp
    phi = euler(1); theta = euler(2); psi = euler(3)

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi)
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi)

    rotmat(1,1) = cos_theta*cos_psi
    rotmat(2,1) = cos_theta*sin_psi
    rotmat(3,1) = -sin_theta

    rotmat(1,2) = sin_phi*sin_theta*cos_psi - cos_phi*sin_psi
    rotmat(2,2) = sin_psi*sin_theta*sin_phi + cos_phi*cos_psi
    rotmat(3,2) = sin_phi*cos_theta

    rotmat(1,3) = cos_phi*sin_theta*cos_psi + sin_phi*sin_psi
    rotmat(2,3) = cos_phi*sin_theta*sin_psi - sin_phi*cos_psi
    rotmat(3,3) = cos_phi*cos_theta

    end subroutine

!*******************************************************************************

subroutine rotmat_xzy(euler, rotmat)

    real(rp), dimension(3), intent(in) :: euler
    real(rp), dimension(3,3), intent(out) :: rotmat
    real(rp) :: phi, theta, psi
    real(rp) :: sin_phi, sin_theta, sin_psi
    real(rp) :: cos_phi, cos_theta, cos_psi

    rotmat = 0.0_rp
    phi = euler(1); theta = euler(2); psi = euler(3)

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi)
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi)

    rotmat(1,1) = cos_theta*cos_psi
    rotmat(2,1) = sin_psi
    rotmat(3,1) = -sin_theta*cos_psi

    rotmat(1,2) = sin_phi*sin_theta - cos_phi*cos_theta*sin_psi
    rotmat(2,2) = cos_phi*cos_psi
    rotmat(3,2) = sin_phi*cos_theta + cos_phi*sin_theta*sin_psi

    rotmat(1,3) = cos_phi*sin_theta + sin_phi*cos_theta*sin_psi
    rotmat(2,3) = -sin_phi*cos_psi
    rotmat(3,3) = cos_phi*cos_theta - sin_phi*sin_theta*sin_psi

    end subroutine

!*******************************************************************************

subroutine rotmat_yxz(euler, rotmat)

    real(rp), dimension(3), intent(in) :: euler
    real(rp), dimension(3,3), intent(out) :: rotmat
    real(rp) :: phi, theta, psi
    real(rp) :: sin_phi, sin_theta, sin_psi
    real(rp) :: cos_phi, cos_theta, cos_psi

    rotmat = 0.0_rp
    phi = euler(1); theta = euler(2); psi = euler(3)

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi)
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi)

    rotmat(1,1) = cos_theta*cos_psi - sin_phi*sin_theta*sin_psi
    rotmat(2,1) = sin_phi*sin_theta*cos_psi + cos_theta*sin_psi
    rotmat(3,1) = -cos_phi*sin_theta

    rotmat(1,2) = -cos_phi*sin_psi
    rotmat(2,2) = cos_phi*cos_psi
    rotmat(3,2) = sin_phi

    rotmat(1,3) = sin_theta*cos_psi + sin_phi*cos_theta*sin_psi
    rotmat(2,3) = sin_theta*sin_psi - sin_phi*cos_theta*cos_psi
    rotmat(3,3) = cos_phi*cos_theta

    end subroutine

!*******************************************************************************

subroutine rotmat_yzx(euler, rotmat)

    real(rp), dimension(3), intent(in) :: euler
    real(rp), dimension(3,3), intent(out) :: rotmat
    real(rp) :: phi, theta, psi
    real(rp) :: sin_phi, sin_theta, sin_psi
    real(rp) :: cos_phi, cos_theta, cos_psi

    rotmat = 0.0_rp
    phi = euler(1); theta = euler(2); psi = euler(3)

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi)
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi)

    rotmat(1,1) = cos_theta*cos_psi
    rotmat(2,1) = sin_phi*sin_theta + cos_phi*cos_theta*sin_psi
    rotmat(3,1) = sin_phi*cos_theta*sin_psi - cos_phi*sin_theta

    rotmat(1,2) = -sin_psi
    rotmat(2,2) = cos_phi*cos_psi
    rotmat(3,2) = sin_phi*cos_psi

    rotmat(1,3) = sin_theta*cos_psi
    rotmat(2,3) = cos_phi*sin_theta*sin_psi - sin_phi*cos_theta
    rotmat(3,3) = sin_phi*sin_theta*sin_psi + cos_phi*cos_theta

    end subroutine

!*******************************************************************************

subroutine rotmat_zxy(euler, rotmat)

    real(rp), dimension(3), intent(in) :: euler
    real(rp), dimension(3,3), intent(out) :: rotmat
    real(rp) :: phi, theta, psi
    real(rp) :: sin_phi, sin_theta, sin_psi
    real(rp) :: cos_phi, cos_theta, cos_psi

    rotmat = 0.0_rp
    phi = euler(1); theta = euler(2); psi = euler(3)

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi)
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi)

     rotmat(1,1) = cos_theta*cos_psi + sin_phi*sin_theta*sin_psi
     rotmat(2,1) = cos_phi*sin_psi
     rotmat(3,1) = sin_phi*cos_theta*sin_psi - sin_theta*cos_psi

     rotmat(1,2) = sin_phi*sin_theta*cos_psi - cos_theta*sin_psi
     rotmat(2,2) = cos_phi*cos_psi
     rotmat(3,2) = sin_phi*cos_theta*cos_psi + sin_theta*sin_psi

     rotmat(1,3) = cos_phi*sin_theta
     rotmat(2,3) = -sin_phi
     rotmat(3,3) = cos_phi*cos_theta

    end subroutine

!*******************************************************************************

subroutine rotmat_zyx (euler, rotmat)

    real(rp), dimension(3), intent(in) :: euler
    real(rp), dimension(3,3), intent(out) :: rotmat
    real(rp)              :: phi, theta, psi
    real(rp)              :: sin_phi, sin_theta, sin_psi
    real(rp)              :: cos_phi, cos_theta, cos_psi

    rotmat = 0.0_rp
    phi = euler(1); theta = euler(2); psi = euler(3)

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi)
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi)

    rotmat(1,1) = cos_theta*cos_psi
    rotmat(2,1) = sin_phi*sin_theta*cos_psi + cos_phi*sin_psi
    rotmat(3,1) = sin_phi*sin_psi - cos_phi*sin_theta*cos_psi

    rotmat(1,2) = -cos_theta*sin_psi
    rotmat(2,2) = cos_phi*cos_psi - sin_phi*sin_theta*sin_psi
    rotmat(3,2) = sin_phi*cos_psi + cos_phi*sin_theta*sin_psi

    rotmat(1,3) = sin_theta
    rotmat(2,3) = -sin_phi*cos_theta
    rotmat(3,3) = cos_phi*cos_theta

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_xyz (rotmat, tol, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), intent(in) :: tol
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi
    real(rp) :: c, s, sphi, cphi

    ! Recover the middle angle theta
    s = -rotmat(3,1)
    c = hypot(rotmat(1,1), rotmat(2,1)) ! = abs(cos theta))
    theta = atan2(s,c) !Returns in [-pi/2, pi/2]
    !Calculate phi (accounting for gimbal lock)
    if (c > tol) then
        phi = atan2(rotmat(3,2), rotmat(3,3))
    else
        phi = 0.0_rp
    end if
    !Calculate psi using phi
    sphi = sin(phi)
    cphi = cos(phi)
    psi = atan2(sphi*rotmat(1,3)-cphi*rotmat(1,2), &
                cphi*rotmat(2,2)-sphi*rotmat(2,3))

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_xzy (rotmat, tol, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), intent(in) :: tol
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi
    real(rp) :: c, s, sphi, cphi

    !Recover angle psi
    s = rotmat(2,1)
    c = hypot(rotmat(1,1), rotmat(3,1)) ! = abs(cos psi))
    psi = atan2(s,c) !Returns in [-pi/2, pi/2]
    !Calculate phi (accounting for gimbal lock)
    if (c > tol) then
        phi = atan2(-rotmat(2,3), rotmat(2,2))
    else
        phi = 0.0_rp
    end if
    !Calculate theta using phi
    sphi = sin(phi)
    cphi = cos(phi)
    theta = atan2(sphi*rotmat(1,2)+cphi*rotmat(1,3), &
                  sphi*rotmat(3,2)+cphi*rotmat(3,3))
    
    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_yxz (rotmat, tol, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), intent(in) :: tol
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi
    real(rp) :: c, s, stheta, ctheta

    !Recover middle angle phi
    s = rotmat(3,2)
    c = hypot(rotmat(1,2), rotmat(2,2)) ! = abs(cos phi))
    phi = atan2(s,c) !Returns in [-pi/2, pi/2]
    !Calculate theta (accounting for gimbal lock)
    if (c > tol) then
        theta = atan2(-rotmat(3,1), rotmat(3,3))
    else
        theta = 0.0_rp
    end if
    !Calculate psi using theta
    stheta = sin(theta)
    ctheta = cos(theta)
    psi = atan2(ctheta*rotmat(2,1)+stheta*rotmat(2,3), &
                ctheta*rotmat(1,1)+stheta*rotmat(1,3))

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_yzx (rotmat, tol, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), intent(in) :: tol
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi
    real(rp) :: c, s, stheta, ctheta

    !Recover middle angle psi
    s = -rotmat(1,2)
    c = hypot(rotmat(2,2), rotmat(3,2)) ! = abs(cos psi)
    psi = atan2(s,c) !Returns in [-pi/2, pi/2]
    !Calculate theta (accounting for gimbal lock)
    if (c > tol) then
        theta = atan2(rotmat(1,3), rotmat(1,1))
    else
        theta = 0.0_rp
    end if
    !Calculate phi using theta
    stheta = sin(theta)
    ctheta = cos(theta)
    phi = atan2(stheta*rotmat(2,1)-ctheta*rotmat(2,3), &
                ctheta*rotmat(3,3)-stheta*rotmat(3,1))

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_zxy (rotmat, tol, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), intent(in) :: tol
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi
    real(rp) :: c, s, spsi, cpsi

    !Recover middle angle phi
    s = -rotmat(2,3)
    c = hypot(rotmat(1,3), rotmat(3,3)) ! = abs(cos phi))
    phi = atan2(s,c) !Returns in [-pi/2, pi/2]
    !Calculate psi (accounting for gimbal lock)
    if (c > tol) then
        psi = atan2(rotmat(2,1), rotmat(2,2))
    else
        psi = 0.0_rp
    end if
    !Calculate theta using psi
    spsi = sin(psi)
    cpsi = cos(psi)
    theta = atan2(spsi*rotmat(3,2)-cpsi*rotmat(3,1), &
                  cpsi*rotmat(1,1)-spsi*rotmat(1,2))

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_zyx (rotmat, tol, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), intent(in) :: tol
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi
    real(rp) :: c, s, spsi, cpsi

    ! Recover the middle angle theta
    s = rotmat(1,3)
    c = hypot(rotmat(2,3), rotmat(3,3)) ! = abs(cos theta))
    theta = atan2(s,c) !Returns in [-pi/2, pi/2]
    !Calculate psi (accounting for gimbal lock)
    if (c > tol) then
        psi = atan2(-rotmat(1,2), rotmat(1,1))
    else
        psi = 0.0_rp
    end if
    !Calculate phi using psi
    spsi = sin(psi)
    cpsi = cos(psi)
    phi = atan2(spsi*rotmat(3,1)+cpsi*rotmat(3,2), &
                spsi*rotmat(2,1)+cpsi*rotmat(2,2))

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!******************************************************************************

end module rotlib_m
