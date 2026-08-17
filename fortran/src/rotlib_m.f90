module rotlib_m
!!
!!  Euler angle ranges:
!!
!!  * XYZ, ZYX: phi in [-pi,   pi],   theta in [-pi/2, pi/2], psi in [-pi,   pi]  
!!  * XZY, YZX: phi in [-pi,   pi],   theta in [-pi,   pi],   psi in [-pi/2, pi/2]
!!  * ZXY, YXZ: phi in [-pi/2, pi/2], theta in [-pi,   pi],   psi in [-pi,   pi]  
!!
!!  Euler angle sequence: XYZ (world). First rotation about X, second rotation
!!  about Y, and the third rotation about Z axis of the world(i.e. fixed) frame.
!!  This is the same as the sequence used in *Blender*.
!!
!!  In contrast, the XYZ sequence is understood in the Aerospace community as:
!!  XYZ: First rotation about Z-axis, second rotation about Y-axis, and the third
!!  rotation about X-axis of the body frame.
!!
!!  *Reference*: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>

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
    !! Rotates vectors such that a set of mutually orthogonal unit vectors
    !! aligns with another set of mutually orthogonal unit vectors.
    !!       
    !! The vectors `v` are rotated using a rotation matrix that, when applied
    !! to the unit vectors `old`, aligns them in the direction of the unit
    !! vectors `new`.
    real(rp), dimension(:,:), intent(in) :: v
        !! *(3,n)*. Vectors to align.
    real(rp), dimension(:,:), intent(in) :: old
        !! *(3,m)*. A set of at most three mutually orthogonal unit vectors.
        !!  `m` must be one of 1, 2 or 3.
    real(rp), dimension(:,:), intent(in) :: new
        !! *(3,m)*. Another set of mutually orthogonal unit vectors of the same
        !! shape as `old`. The angles between the vectors must be the same as
        !! those of `old`.
    real(rp), dimension(:,:), intent(out) :: valigned
        !! *(3,n)*. Aligned vectors, of the same shape as `v`.

    real(rp) :: axis(3), z_old(3), z_new(3), axes_old(3,3), &
                axes_new(3,3), dcm(3,3)
    real(rp) :: cos_angle, angle
    integer :: n

    if ( any(shape(old) /= shape(new)) ) then
        write(*,*) 'Shape mismatch between `old` and `new`.'
        write(*,*) 'Shape of `old` = ', shape(old)
        write(*,*) 'Shape of `new` = ', shape(new)
        stop 'Shape mismatch.'
    end if
    n = size(old, dim=2)
    if (n == 1) then
        !Angle between old and new
        cos_angle = dot_product(old(:,1), new(:,1))
        if ( (cos_angle > 1.0_rp) .or. (cos_angle < -1.0_rp) ) then
            if ( isclose(abs(cos_angle), 1.0_rp, 1e-9_rp, 1e-15_rp) ) then
                cos_angle = sign(1.0_rp, cos_angle)
            else
                stop 'Cos of angle out-of-range.'
            end if
        end if

        if ( isclose(cos_angle, 1.0_rp, 1e-9_rp, 1e-15_rp) ) then
            !Angle = zero: No rotation needed
            valigned = v
        else if ( isclose(cos_angle, -1.0_rp, 1e-9_rp, 1e-15_rp) ) then
            !Angle = pi: Rotation axis not unique. We flip the 
            ! vectors `v` to align.
            valigned = -v
        else
            angle = acos(cos_angle)
            !Axis of rotation 
            call cross(old, new, axis)
            axis = axis/norm2(axis)
            call aa_rotate_vectors(v, axis, angle, valigned)
        end if
    else if (n == 2) then
        !The third vector orthogonal to `old`
        call cross( old(:,1), old(:,2), z_old )
        z_old = z_old/norm2(z_old)
        !The third vector orthogonal to `new`
        call cross( new(:,1), new(:,2), z_new )
        z_new = z_new/norm2(z_new)
        !The `old` and new `frames`
        axes_old(:,1:2) = old; axes_old(:,3) = z_old
        axes_new(:,1:2) = new; axes_new(:,3) = z_new
        call dcm_from_axes(axes_old, axes_new, dcm)
        call dcm_rotate_vectors(v, dcm, valigned)
    else if (n == 3) then
        call dcm_from_axes(old, new, dcm)
        call dcm_rotate_vectors(v, dcm, valigned)
    end if

    end subroutine

!************************************************************************

logical function mat_is_rotmat(mat)
    !! Returns `.true.` if a (3,3) matrix is a rotation matrix.

    real(rp), dimension(3,3), intent(in) :: mat
        !! (3,3). Matrix to check.
    real(rp), dimension(3,3) :: mat_trans
    real(rp), dimension(3,3) :: mat_prod
    logical :: det_is_one
    logical :: is_orthogonal
    real(rp), dimension(3,3) :: eye

    call identity(eye)
    mat_prod = matmul(mat, transpose(mat))

    det_is_one = isclose(det(mat), 1.0_rp)
    is_orthogonal = allclose(mat_prod, eye)
    mat_is_rotmat = (is_orthogonal .and. det_is_one)

    end function

!******************************************************************************!
!                                  QUATERNION                                  !
!******************************************************************************!

subroutine quat_rand(q)
    !! Returns a random unit quaternion.
    real(rp), dimension(4), intent(out) :: q
        !! Unit quaternion.
    real(rp), dimension(3) :: axis
    real(rp) :: angle

    call aa_rand(axis, angle)
    call aa_to_quat(axis, angle, q)

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

subroutine quat_inverted(q)
    !! Inverts a quaternion in-place.

    real(rp), dimension(4), intent(in out) :: q
        !! On entry, a quaternion. On return, the inverted quaternion.

    call quat_conjugated(q)
    call quat_normalized(q)

    end subroutine

!************************************************************************

subroutine quat_normalized(q)
    !! Normalizes a quaternion in-place.

    real(rp), dimension(4), intent(in out) :: q
        !! On entry, a quaternion. On return, the normalized quaternion.

    q = q /norm2(q)

    end subroutine

!************************************************************************

logical function quat_is_normalized(q)
    !! Returns `.true.` if quaternion is normalized, i.e. if it is a unit
    !! quaternion

    real(rp), dimension(4), intent(in) :: q
        !! A quaternion
    real(rp) :: norm

    norm = norm2(q)
    if (isclose(norm, 1.0_rp, rel_tol=1e-14_rp)) then
        quat_is_normalized =  .true.
    else
        quat_is_normalized =  .false.
    end if

    end function

!************************************************************************

subroutine quat_prod (p, q, pq, normalize)
    !! Calculates the product of two quaternions

    real(rp), dimension(4), intent(in)  :: p
        !! Quaternion.
    real(rp), dimension(4), intent(in)  :: q
        !! Quaternion.
    real(rp), dimension(4), intent(out) :: pq
        !! Product of `p` and `q`. This is a unit quaternion if `normalize`
        !! is `.true.`.
    logical, intent(in), optional :: normalize
        !! Whether to normalize the product, defaults to `.true.`.
    logical :: normalize_ = .true.
    real(rp), dimension(4,4) :: prod_mat

    if (present(normalize)) normalize_ = normalize

    pq(1) = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)
    pq(2) = p(2)*q(1) + p(1)*q(2) - p(4)*q(3) + p(3)*q(4)
    pq(3) = p(3)*q(1) + p(4)*q(2) + p(1)*q(3) - p(2)*q(4)
    pq(4) = p(4)*q(1) - p(3)*q(2) + p(2)*q(3) + p(1)*q(4)

    if (normalize_) call quat_normalized(pq)

    end subroutine

!************************************************************************

function quat_angle_between(p, q) result(angle)
    !! Returns the angle between two _unit_ quaternions

    real(rp), dimension(4), intent(in) :: p
        !! Unit quaternion.
    real(rp), dimension(4), intent(in) :: q
        !! Unit quaternion.
    real(rp) :: angle, cos_angle
        !! Angle in radian.

    cos_angle = dot_product(p, q)
    if ( (cos_angle > 1.0_rp) .or. (cos_angle < -1.0_rp) ) then
        if ( isclose(abs(cos_angle), 1.0_rp, 1e-9_rp, 1e-15_rp) ) then
            cos_angle = sign(1.0_rp, cos_angle)
        else
            stop 'Cos of angle out-of-range.'
        end if
    end if
    angle = acos(cos_angle)

    end function

!************************************************************************

subroutine quat_interpolate(q1, q2, t, q)
    !! Interpolates between two _unit_ quaternions.

    real(rp), dimension(4), intent(in) :: q1
        !! Unit quaternion.
    real(rp), dimension(4), intent(in) :: q2
        !! Unit quaternion.
    real(rp),               intent(in)  :: t
        !! Interpolation factor, *0 <= t <= 1*.
    real(rp), dimension(4), intent(out) :: q
        !! Interpolated unit quaternion.
    real(rp) :: theta

    theta = quat_angle_between(q1, q2)

    q = (q1*sin((1.0_rp-t)*theta) + q2*sin(t*theta))/sin(theta)

    end subroutine

!************************************************************************

subroutine quat_deriv_to_angvel_mat(qdot, mat)
    !! Returns the matrix mapping the derivative of a unit quaternion to
    !! angular velocity.
    real(rp), dimension(4),   intent(in)  :: qdot
        !! Derivative of a unit quaternion.
    real(rp), dimension(3,4), intent(out) :: mat
        !! Angular velocity matrix.

    mat = reshape([-qdot(2), -qdot(3), -qdot(4), &
                    qdot(1),  qdot(4), -qdot(3), &
                   -qdot(4),  qdot(1),  qdot(2), &
                    qdot(3), -qdot(2),  qdot(1) ], [3,4])
    mat = 2*mat

    end subroutine

!************************************************************************

subroutine quat_deriv_to_angvel (q, qdot, angvel)
    !! Calculates the angular velocity from a unit quaternion and its time
    !! derivative.
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(4), intent(in)  :: qdot
        !! Derivative of `q`.
    real(rp), dimension(3), intent(out) :: angvel
        !! Angular velocity.
    real(rp) :: mat(3,4)

    call quat_deriv_to_angvel_mat(q, mat)
    angvel = matmul(mat, qdot)

    end subroutine

!************************************************************************

subroutine quat_deriv_from_angvel_mat (q, mat)
    !! Returns the matrix mapping angular velocity to time derivative of a unit
    !! quaternion.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(4,3), intent(out) :: mat
        !! Quaternion derivative matrix.

    mat = reshape([-q(2),  q(1), -q(4),  q(3), &
                   -q(3),  q(4),  q(1), -q(2), &
                   -q(4), -q(3),  q(2),  q(1)], [4,3])
    mat = mat/2

    end subroutine

!************************************************************************

subroutine quat_deriv_from_angvel(q, angvel, qdot)
    !! Calculates the time derivative of a unit quaternion from angular velocity.
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(3), intent(in)  :: angvel
        !! Angular velocity.
    real(rp), dimension(4), intent(out) :: qdot
        !! Time derivative of `q`.
    real(rp) :: mat(4,3)

    call quat_deriv_from_angvel_mat(q, mat)
    qdot = matmul(mat, angvel)

    end subroutine

!************************************************************************

subroutine quat_rotmat (q, rotmat)
    !! Returns the rotation matrix corresponding to a unit quaternion.
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

    rotmat(1,1) = 2*(q1sq + q2sq) - 1.0_rp
    rotmat(2,1) = 2*(q2q3 + q1q4)
    rotmat(3,1) = 2*(q2q4 - q1q3)

    rotmat(1,2) = 2*(q2q3 - q1q4)
    rotmat(2,2) = 2*(q1sq + q3sq) - 1.0_rp
    rotmat(3,2) = 2*(q3q4 + q1q2)

    rotmat(1,3) = 2*(q2q4 + q1q3)
    rotmat(2,3) = 2*(q3q4 - q1q2)
    rotmat(3,3) = 2*(q1sq + q4sq) - 1.0_rp

    end subroutine

!************************************************************************

subroutine quat_shiftmat(q, forward, shiftmat)
    !! Returns the shifter matrix corresponding to a unit quaternion.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    logical,                  intent(in)  :: forward
        !! Whether to shift forward, i.e., along the orientation or
        !! shift reverse.
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
    !! Rotates vectors by a unit quaternion.
    real(rp), dimension(:,:), intent(in)  :: v
        !! *(3,n)*. Vectors to rotate.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! *(3,n)* Rotated vectors.
    real(rp), dimension(3,3) :: rotmat

    call quat_rotmat(q, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine quat_shift_vectors(v, q, forward, vshift)
    !! Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by a unit quaternion, shifts vectors from A
    !! to B or B to A.
    real(rp), dimension(:,:), intent(in)  :: v
        !! *(3,n)*. Vectors to shift.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    logical,                  intent(in)  :: forward
        !! If `.true.`, shift from A to B. If `.false.`, shift from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! *(3,n)* Shifted vectors.
    real(rp), dimension(3,3) :: shiftmat

    call quat_shiftmat(q, forward, shiftmat)
    vshift = matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine quat_shift_tensor2(a, q, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B with respect
    !! to frame A is given by a unit quaternion, shifts second-order tensors from A
    !! to B or B to A.
    real(rp), dimension(3,3), intent(in)  :: a
        !! A second-order tensor.
    real(rp), dimension(4),   intent(in)  :: q
        !! Unit quaternion.
    logical,                  intent(in)  :: forward
        !! If `.true.`, shift from A to B. If `.false.`, shift from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Shifted second-order tensor.
    real(rp), dimension(3,3) :: shiftmat

    call quat_shiftmat(q, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine quat_shift_tensor3(a, q, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B with respect
    !! to frame A is given by a unit quaternion, shifts third-order tensors from A
    !! to B or B to A.
    real(rp), dimension(3,3,3), intent(in)  :: a
        !! A third-order tensor.
    real(rp), dimension(4),     intent(in)  :: q
        !! Unit quaternion.
    logical,                    intent(in)  :: forward
        !! If `.true.`, shift from A to B. If `.false.`, shift from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Shifted third-order tensor.
    real(rp), dimension(3,3) :: shiftmat

    call quat_shiftmat(q, forward, shiftmat)
    call shift_tensor3(a, shiftmat, ashift)

    end subroutine 

!************************************************************************

subroutine quat_to_aa (q, axis, angle)
    !! Converts a unit quaternion to an *axis-angle* representation.
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    real(rp), dimension(3), intent(out) :: axis
        !! Unit vector along the axis of rotation.
    real(rp),               intent(out) :: angle
        !! Angle in radian.
    real(rp), dimension(3,3) :: rotmat(3,3)
    real(rp) :: sine

    angle = 2*acos(q(1))
    sine = sqrt(1.0_rp-q(1)*q(1))

    if (angle > 0.0_rp) then
        if (angle < math_pi) then
            axis = q(2:4)/sine
        else
            call quat_rotmat(q, rotmat)
            call aa_from_rotmat(rotmat, axis, angle)
        end if
    else
        axis = [1.0_rp, 0.0_rp, 0.0_rp]
    end if

    call aa_fix(axis, angle, .true.)

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

subroutine quat_to_euler (q, seq, world, euler)
    !! Converts a unit quaternion to an Euler angle sequence.
    real(rp), dimension(4), intent(in)  :: q
        !! Unit quaternion.
    character(3),           intent(in)  :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3), intent(out) :: euler
        !! Euler angles.
    real(rp), dimension(3,3) :: rotmat

    call quat_rotmat(q, rotmat)
    call factor_rotmat(rotmat, seq, world, euler)

    end subroutine

!*****************************************************************************!
!                                  AXIS ANGLE                                 !
!*****************************************************************************!

subroutine aa_fix (axis, angle, normalize)
    !! Modifies `axis` and `angle` to ensure a right handed rotation with
    !! `angle` \(\in [0, \pi)\).
    real(rp), dimension(3), intent(in out) :: axis
        !! Axis of rotation. If this is not a unit vector, set `normalize` to
        !! ``True``.
    real(rp), intent(in out) :: angle
        !! Angle in radian.
    logical, intent(in), optional :: normalize
        !! Whether to normalize the axis to a unit vector.
    logical :: normalize_
    real(rp) :: norm

    normalize_ = .true.
    if (present(normalize)) normalize_ = normalize
    if (normalize_) call unitize(axis)
    angle = mod(angle, 2*math_pi)
    if (angle < 0.0_rp) then
        angle = -angle
        axis = -axis
    end if
    if (angle > math_pi) then
        angle = 2*math_pi - angle
        axis = -axis
    end if

    end subroutine

!************************************************************************

subroutine aa_rand (axis, angle)
    !! Generates a random orientation in *axis-angle* representation.
    !! The axis is a random vector drawn from a uniform distribution on the
    !! surface of a unit sphere. The current implementation in based on the
    !! algorithm from Allen & Tildesley p. 349.
    real(rp), dimension(3), intent(out) :: axis
    real(rp), intent(out) :: angle
    real(rp) :: zeta1, zeta2, zetasq
    real(rp) :: rt

    axis = 0.0_rp
   !Generate angle: A uniform random number from [0, 2*pi)
    call random_number(angle)
    angle = 2*math_pi*angle

    do
        !Generate two uniform random numbers from [-1, 1)
        call random_number(zeta1)
        call random_number(zeta2)
        zeta1 = 2*zeta1 - 1.0_rp
        zeta2 = 2*zeta2 - 1.0_rp
        zetasq = zeta1**2 + zeta2**2
        if (zetasq <= 1._rp) exit
    end do

    rt = sqrt(1.0_rp-zetasq)
    axis(1) = 2*zeta1*rt
    axis(2) = 2*zeta2*rt
    axis(3) = 1.0_rp - 2*zetasq

    call aa_fix(axis, angle, .true.)

    end subroutine

!************************************************************************

subroutine aa_rotmat (axis, angle, rotmat)
    !! Returns the rotation matrix corresponding to an *axis-angle*
    !! representation.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. `0 <= angle < pi`.
    real(rp), dimension(3,3), intent(out) :: rotmat
        !! Rotation matrix.
    real(rp) :: sine, cosine, icos

    sine = sin(angle)
    cosine = cos(angle)
    icos = 1.0_rp - cosine

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
    real(rp), dimension(3,3), intent(in) :: rotmat
        !! Rotation matrix (must be orthonormal).
    real(rp), dimension(3), intent(out) :: axis
        !! Unit vector along the axis.
    real(rp), intent(out) :: angle
        !! Angle of rotation in radian.
    real(rp) :: tr
    real(rp) :: s
    integer  :: k

    tr = rotmat(1,1) + rotmat(2,2) + rotmat(3,3)
    angle = acos((tr-1.0_rp)/2)

    if (angle > 0.0_rp) then
        if (angle < math_pi) then
            axis(1) = rotmat(3,2) - rotmat(2,3)
            axis(2) = rotmat(1,3) - rotmat(3,1)
            axis(3) = rotmat(2,1) - rotmat(1,2)
        else
            !Find the largest entry in the diagonal of rotmat
            k = maxloc([rotmat(1,1), rotmat(2,2), rotmat(3,3)], 1)

            select case (k)
                case (1)
                    axis(1) = sqrt(rotmat(1,1)-rotmat(2,2)-rotmat(3,3)+1.0_rp)/2
                    s = 1.0_rp/(2*axis(1))
                    axis(2) = s*rotmat(1,2)
                    axis(3) = s*rotmat(1,3)
                case (2)
                    axis(2) = sqrt(rotmat(2,2)-rotmat(1,1)-rotmat(3,3)+1.0_rp)/2
                    s = 1.0_rp/(2*axis(2))
                    axis(1) = s*rotmat(1,2)
                    axis(3) = s*rotmat(2,3)
                case (3)
                    axis(3) = sqrt(rotmat(3,3)-rotmat(1,1)-rotmat(2,2)+1.0_rp)/2
                    s = 1.0_rp/(2*axis(3))
                    axis(1) = s*rotmat(1,3)
                    axis(2) = s*rotmat(2,3)
            end select
        end if
    else
        axis = [1.0_rp, 0.0_rp, 0.0_rp]
    end if

    !Ensuring that angle lies in [0, pi)
    call aa_fix(axis, angle, .true.)

    end subroutine

!************************************************************************

subroutine aa_rotate_vectors (v, axis, angle, vrot)
    !! Rotates vectors about `axis` by `angle`.
    real(rp), dimension(:,:), intent(in) :: v
        !! *(3,n)*. Vectors to rotate.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. `0 <= angle < pi`.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! *(3,n)*. Rotated vectors.
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
        !! Angle of rotation in radian. `0 <= angle < pi`.
    logical, intent(in) :: forward
        !! Whether to shift forward, i.e., along the orientation or shift reverse.
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
    !! Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by an *axis-angle* representation, shifts
    !! vectors from A to B or B to A.
    real(rp), dimension(:,:), intent(in) :: v
        !! *(3,n)*. Vectors to shift.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. `0 <= angle < pi`.
    logical, intent(in) :: forward
        !! If ``True``, shift from A to B. If ``False``, shift from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! *(3,n)*. Shifted vectors.
    real(rp), dimension(3,3) :: shiftmat

    call aa_shiftmat(axis, angle, forward, shiftmat)
    vshift = matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine aa_shift_tensor2 (a, axis, angle, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by an axis-angle representation, shifts
    !! second-order tensors from A to B or B to A.
    real(rp), dimension(3,3), intent(in) :: a
        !! Second-order tensor
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. `0 <= angle < pi`.
    logical, intent(in) :: forward
        !! If ``True``, shift from A to B. If ``False``, shift from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Shifted second-order tensor.
    real(rp), dimension(3,3) :: shiftmat

    call aa_shiftmat(axis, angle, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine aa_shift_tensor3 (a, axis, angle, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by an axis-angle representation, shifts
    !! third-order tensors from A to B or B to A.
    real(rp), dimension(3,3,3), intent(in) :: a
        !! Third-order tensor
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle of rotation in radian. `0 <= angle < pi`.
    logical, intent(in), optional :: forward
        !! If ``True``, shift from A to B. If ``False``, shift from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Shifted third-order tensor.
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
        !! Angle of rotation in radian. `0 <= angle < pi`.
    real(rp), dimension(4), intent(out) :: q
        !! Unit quaternion.

    q(1) = cos(angle/2)
    q(2:4) = sin(angle/2)*axis

    !Normalize to ensure unit quaternion
    call quat_normalized(q)

    end subroutine

!************************************************************************

subroutine aa_to_dcm (axis, angle, dcm)
    !! Converts an *axis-angle* representation to a direction cosine matrix.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along  the axis.
    real(rp), intent(in)   :: angle
        !! Angle of rotation in radian. `0 <= angle < pi`.
    real(rp), dimension(3,3), intent(out) :: dcm
        !! Direction cosine matrix.

    call aa_shiftmat(axis, angle, .true., dcm)

    end subroutine

!************************************************************************

subroutine aa_to_euler (axis, angle, seq, world, euler)
    !! Coverts an *axis-angle* representation to *Euler angles*.
    real(rp), dimension(3), intent(in) :: axis
        !! Unit vector along the axis.
    real(rp), intent(in) :: angle
        !! Angle in radian. `0 <= angle < pi`.
    character(len=3), intent(in) :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3), intent(out) :: euler
        !! Euler angles.
    real(rp), dimension(3,3) :: rotmat
    
    call aa_rotmat(axis, angle, rotmat)
    call factor_rotmat(rotmat, seq, world, euler)

    end subroutine

!*****************************************************************************!
!                           DIRECTION COSINE MATRIX                           !
!*****************************************************************************!

logical function mat_is_dcm (mat)
    !! Returns *true* if a (3,3) matrix is a direction cosine matrix.
    real(rp), dimension(3,3), intent(in) :: mat
        !! Matrix to check.

    mat_is_dcm = mat_is_rotmat(mat)

    end function

!************************************************************************

subroutine dcm_from_axes (A, B, dcm)
   !!  Returns the direction cosine matrix of axes(i.e. frame) B with
   !!  respect to axes(i.e. frame) A.
    real(rp), dimension(3,3), intent(in)  :: A
       !! The columns of A represent the orthonormal basis vectors of frame A.
    real(rp), dimension(3,3), intent(in)  :: B
       !! The columns of B represent the orthonormal basis vectors of frame B.
    real(rp), dimension(3,3), intent(out) :: dcm
        !! Direction cosine matrix of frame B with respect to frame A

    dcm =  matmul(transpose(A), B)

    end subroutine

!************************************************************************

subroutine dcm_rotmat (dcm, rotmat)
    !! Returns the rotation matrix corresponding to a direction cosine matrix.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(3,3), intent(out) :: rotmat
        !! Rotation matrix.

    rotmat = transpose(dcm)

    end subroutine

!************************************************************************

subroutine dcm_shiftmat (dcm, forward, shiftmat)
    !! Returns the shifter matrix corresponding to a direction cosine matrix.
    real(rp), dimension(3,3), intent(in)  :: dcm
        !! Direction cosine matrix.
    logical, intent(in) :: forward
        !! Whether to shift forward, i.e., along the orientation or shift reverse.
    real(rp), dimension(3,3), intent(out) :: shiftmat
        !! Shifter matrix.

    if (forward) then
        shiftmat = dcm
    else
        shiftmat = transpose(shiftmat)
    end if

    end subroutine

!************************************************************************

subroutine dcm_rotate_vectors (v, dcm, vrot)
    !! Rotates vectors by a direction cosine matrix.
    real(rp), dimension(:,:), intent(in) :: v
        !! *(3,n)*. Vectors to rotate.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! *(3,n)*. Rotated vectors.
    real(rp) :: rotmat(3,3)

    call dcm_rotmat(dcm, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine dcm_shift_vectors (v, dcm, forward, vshift)
    !!  Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by a direction cosine matrix, shifts
    !! vectors from A to B or B to A.
    real(rp), dimension(:,:), intent(in) :: v
        !! *(3,n)*. Vectors to shift.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    logical, intent(in) :: forward
        !! If ``True``, shift from A to B. If ``False``, shift from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! *(3,n)*. Shifted vectors.
    real(rp), dimension(3,3) :: shiftmat

    call dcm_shiftmat(dcm, forward, shiftmat)
    vshift = matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine dcm_shift_tensor2 (a, dcm, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by a direction cosine matrix, shifts
    !! second-order tesnors from A to B or B to A.
    real(rp), dimension(3,3), intent(in) :: a
        !! A second-order tensor.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    logical, intent(in) :: forward
        !! If ``True``, shift from A to B. If ``False``, shift from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! *(3,n)*. Shifted second-order tensor.
    real(rp), dimension(3,3) :: shiftmat

    call dcm_shiftmat(dcm, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine dcm_shift_tensor3 (a, dcm, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B with
    !! respect to frame A is given by a direction cosine matrix, shifts
    !! third-order tesnors from A to B or B to A.
    real(rp), dimension(3,3,3), intent(in) :: a
        !! A second-order tensor.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    logical, intent(in) :: forward
        !! If ``True``, shift from A to B. If ``False``, shift from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! *(3,n)*. Shifted second-order tensor.
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
        !! Angle in radian.
    real(rp), dimension(3,3) :: rotmat

    call dcm_rotmat(dcm, rotmat)
    call aa_from_rotmat(rotmat, axis, angle)

    end subroutine

!************************************************************************

subroutine dcm_to_euler (dcm, seq, world, euler)
    !! Converts a direction cosine matrix to Euler angles.
    real(rp), dimension(3,3), intent(in) :: dcm
        !! Direction cosine matrix.
    character(3), intent(in) :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3), intent(out)  :: euler
        !! Euler angles
    real(rp), dimension(3,3) :: rotmat

    call dcm_rotmat (dcm, rotmat)
    call factor_rotmat(rotmat, seq, world, euler)

    end subroutine

!*****************************************************************************!
!                                 EULER ANGLES                                !
!*****************************************************************************!

subroutine euler_rotmat (euler, seq, world, rotmat)
    !!  Returns the rotation matrix for a set of Euler angles.
    real(rp), dimension(3), intent(in) :: euler
    character(len=3), intent(in) :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
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
        stop 'unknown value of seq'
    end select

    if (.not. world) rotmat = transpose(rotmat)

    end subroutine

!*******************************************************************************

subroutine euler_shiftmat (euler, seq, world, forward, shiftmat)
    !! Returns the shifter matrix for a set of Euler angles.
    real(rp), dimension(3), intent(in) :: euler
        !! Euler angles.
    character(3),           intent(in) :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! Whether to shift forward, i.e., along the orientation or shift reverse.
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
    !! Rotates vectors with a set of Euler angles.
    real(rp), dimension(:,:), intent(in)  :: v
        !! *(3,n)*. Vectors to rotate.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(:,:), intent(out) :: vrot
        !! *(3,n)*. Rotated vectors.
    real(rp), dimension(3,3) :: rotmat

    call euler_rotmat (euler, seq, world, rotmat)
    vrot = matmul(rotmat, v)

    end subroutine

!************************************************************************

subroutine euler_shift_vectors (v, euler, seq, world, forward, vshift)
    !! Given two frames A and B such that the orientation of frame B
    !! with respect to frame A is given by three Euler angles, shifts
    !! vectors from A to B or B to A.
    real(rp), dimension(:,:), intent(in)  :: v
        !! *(3,n)*. Vectors to rotate.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, shift from A to B. If `.false.`, shift from B to A.
    real(rp), dimension(:,:), intent(out) :: vshift
        !! *(3,n)*. Shifted vectors.
    real(rp), dimension(3,3) :: shiftmat

    call euler_shiftmat (euler, seq, world, forward, shiftmat)
    vshift =  matmul(shiftmat, v)

    end subroutine

!************************************************************************

subroutine euler_shift_tensor2 (a, euler, seq, world, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B
    !! with respect to frame A is given by three Euler angles, shifts
    !! second-order tensors from A to B or B to A.
    real(rp), dimension(3,3), intent(in)  :: a
        !! Second-order tensor.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, shift from A to B. If `.false.`, shift from B to A.
    real(rp), dimension(3,3), intent(out) :: ashift
        !! Shifted second-order tensor.
    real(rp), dimension(3,3) :: shiftmat

    call euler_shiftmat (euler, seq, world, forward, shiftmat)
    call shift_tensor2(a, shiftmat, ashift)

    end subroutine

!************************************************************************

subroutine euler_shift_tensor3 (a, euler, seq, world, forward, ashift)
    !! Given two frames A and B such that the orientation of frame B
    !! with respect to frame A is given by three Euler angles, shifts
    !! third-order tensors from A to B or B to A.
    real(rp), dimension(3,3,3), intent(in)  :: a
        !! Third-order tensor.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    logical,                  intent(in)  :: forward
        !! If `.true.`, shift from A to B. If `.false.`, shift from B to A.
    real(rp), dimension(3,3,3), intent(out) :: ashift
        !! Shifted third-order tensor.
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
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(4),     intent(out) :: q
        !! Unit quaternion.
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
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3),     intent(out) :: axis
        !! Unit vector along the axis of rotation.
    real(rp),                   intent(out) :: angle
        !! Angle of rotation in radian.
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
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3,3),   intent(out) :: dcm
        !! Direction cosine matrix.
    real(rp), dimension(3,3) :: rotmat

    call euler_shiftmat (euler, seq, world, .true., dcm)

    end subroutine

!************************************************************************

subroutine euler_to_euler (euler, seq, world, to_seq, to_world, to_euler)
    !! Convert one set of Euler angles to another.
    real(rp), dimension(3),   intent(in)  :: euler
        !! Euler angles.
    character(3),             intent(in)  :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Sequence for `euler`.
    logical,                  intent(in)  :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
    character(3),               intent(out) :: to_seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Convert `euler` with sequence `seq` to the sequence `to_seq`.
    logical,                    intent(out) :: to_world
        !! Whether the converted euler angles are with respect to the *world* frame or not.
    real(rp), dimension(3),     intent(out) :: to_euler
        !! Converted Euler angles.
    real(rp), dimension(3,3) :: rotmat

    call euler_rotmat (euler, seq, world, rotmat)
    call factor_rotmat(rotmat, to_seq, to_world, to_euler)

    end subroutine

!*****************************************************************************!
!                           INTERNAL ROUTINES                                 !
!*****************************************************************************!

subroutine factor_rotmat (rotmat, seq, world, euler)
    !! Factorize a rotation matrix to obtain the three Euler angles.
    real(rp), dimension(3,3), intent(in) :: rotmat
        !! Rotation matrix.
    character(3), intent(in) :: seq
        !! `'XYZ'`|`'XZY'` | `'YXZ'` | `'YZX'` | `'ZXY'` | `'ZYX'`.
        !! Euler angle sequence.
    logical, intent(in) :: world
        !! Whether the Euler angles are with respect to the *world* frame or not.
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
        call factor_rotmat_xyz(rotmat_, euler)
    case ('XZY')
        call factor_rotmat_xzy(rotmat_, euler)
    case ('YXZ')
        call factor_rotmat_yxz(rotmat_, euler)
    case ('YZX')
        call factor_rotmat_yzx(rotmat_, euler)
    case ('ZXY')
        call factor_rotmat_zxy(rotmat_, euler)
    case ('ZYX')
        call factor_rotmat_zyx(rotmat_, euler)
    case default
        stop 'Unknown value of seq'
    end select

    if (.not. world) euler = -euler

    end subroutine

!************************************************************************

subroutine shift_tensor2(a, shiftmat, ashift)
    !! Shifts a second order tensor, i.e. (3,3) matrix using a shifter matrix

    real(rp), dimension(3,3), intent(in) :: a
    real(rp), dimension(3,3), intent(in) :: shiftmat
    real(rp), dimension(3,3), intent(out) :: ashift

    ashift = matmul(shiftmat, matmul(a, transpose(shiftmat)))

    end subroutine

!************************************************************************

subroutine shift_tensor3(a, shiftmat, ashift)
    !! Shifts a third order tensor, i.e. (3,3,3) matrix using a shifter matrix

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

subroutine factor_rotmat_xyz (rotmat, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi

    if (rotmat(3,1) < 1.0_rp) then
        if (rotmat(3,1) > -1.0_rp) then
            theta = asin(-rotmat(3,1))
            psi   = atan2(rotmat(2,1), rotmat(1,1))
            phi   = atan2(rotmat(3,2), rotmat(3,3))
        else
            !Not unique: phi - psi = atan2(-rotmat(2,3), rotmat(2,2))
            theta = math_pi_2
            psi = -atan2(-rotmat(2,3), rotmat(2,2))
            phi = 0.0_rp
        end if
    else
        !Not unique: phi + psi = atan2(-rotmat(2,3), rotmat(2,2))
        phi = 0.0_rp
        theta = -math_pi_2
        psi = atan2(-rotmat(2,3), rotmat(2,2))
    end if

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_xzy (rotmat, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi

    if (rotmat(2,1) < 1.0_rp) then
        if (rotmat(2,1) > -1.0_rp) then
            phi   = atan2(-rotmat(2,3), rotmat(2,2))
            theta = atan2(-rotmat(3,1), rotmat(1,1))
            psi   = asin(rotmat(2,1))
        else
            !Not unique: phi - theta = atan2(rotmat(3,2), rotmat(3,3))
            phi = 0.0_rp
            theta = -atan2(rotmat(3,2), rotmat(3,3))
            psi = -math_pi_2
        end if
    else
        !Not unique: phi + theta = atan2(rotmat(3,2), rotmat(3,3))
        phi = 0.0_rp
        theta = atan2(rotmat(3,2), rotmat(3,3))
        psi = math_pi_2
    end if

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_yxz (rotmat, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi

    if (rotmat(3,2) < 1.0_rp) then
        if (rotmat(3,2) > -1.0_rp) then
            phi   = asin(rotmat(3,2))
            theta = atan2(-rotmat(3,1), rotmat(3,3))
            psi   = atan2(-rotmat(1,2), rotmat(2,2))
        else
            !Not unique: theta - psi = atan2(rotmat(1,3), rotmat(1,1))
            phi = -math_pi_2
            theta = 0.0_rp
            psi = -atan2(rotmat(1,3), rotmat(1,1))
        end if
    else
        !Not unique: theta + psi = atan2(rotmat(1,3), rotmat(1,1))
        phi = math_pi_2
        theta = 0.0_rp
        psi = atan2(rotmat(1,3), rotmat(1,1))
    end if

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_yzx (rotmat, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi

    if (rotmat(1,2) < 1.0_rp) then
        if (rotmat(1,2) > -1.0_rp) then
            phi   = atan2(rotmat(3,2), rotmat(2,2))
            theta = atan2(rotmat(1,3), rotmat(1,1))
            psi   = asin(-rotmat(1,2))
        else
            !Not unique: theta - phi = atan2(-rotmat(3,1), rotmat(3,3))
            phi = -atan2(-rotmat(3,1), rotmat(3,3))
            theta = 0.0_rp
            psi = math_pi_2
        end if
    else
        !Not unique: theta + phi = atan2(-rotmat(3,1), rotmat(3,3))
        phi = atan2(-rotmat(3,1), rotmat(3,3))
        theta = 0.0_rp
        psi = -math_pi_2
    end if

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_zxy (rotmat, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi

    if (rotmat(2,3) < 1.0_rp) then
        if (rotmat(2,3) > -1.0_rp) then
            phi   = asin(-rotmat(2,3))
            theta = atan2(rotmat(1,3), rotmat(3,3))
            psi   = atan2(rotmat(2,1), rotmat(2,2))
        else
            !Not unique: psi - theta = atan2(-rotmat(1,2), rotmat(1,1))
            phi = math_pi_2
            theta = -atan2(-rotmat(1,2), rotmat(1,1))
            psi = 0.0_rp
        end if
    else
        !Not unique: psi + theta = atan2(-rotmat(1,2), rotmat(1,1))
        phi = -math_pi_2
        theta = atan2(-rotmat(1,2), rotmat(1,1))
        psi = 0.0_rp
    end if

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!*******************************************************************************

subroutine factor_rotmat_zyx (rotmat, euler)

    real(rp), dimension(3,3), intent(in) :: rotmat
    real(rp), dimension(3), intent(out) :: euler
    real(rp) :: phi, theta, psi

    if (rotmat(1,3) < 1.0_rp) then
        if (rotmat(1,3) > -1.0_rp) then
            phi   = atan2(-rotmat(2,3), rotmat(3,3))
            theta = asin(rotmat(1,3))
            psi   = atan2(-rotmat(1,2), rotmat(1,1))
        else 
            !Not unique: psi - phi = atan2(rotmat(2,1), rotmat(2,2))
            phi = -atan2(rotmat(2,1), rotmat(2,2))
            theta = -math_pi_2
            psi = 0.0_rp
        end if
    else
        !Not unique: psi + phi = atan2(rotmat(2,1), rotmat(2,2))
        phi = atan2(rotmat(2,1), rotmat(2,2))
        theta = math_pi_2
        psi = 0.0_rp
    end if

    euler(1) = phi; euler(2) = theta; euler(3) = psi

    end subroutine

!******************************************************************************

end module rotlib_m
