program test_main

use m_precision
use m_utils_math
use m_rotlib

implicit none

! Reset random number generator with an arbitrary seed
call random_seed()

write(*,*) "TESTING ROTLIB"
call test_quat_to_aa_roundtrip()
call test_quat_to_dcm_roundtrip()
call test_quat_to_euler_roundtrip()
call test_aa_to_dcm_roundtrip()
call test_aa_to_euler_roundtrip()
call test_dcm_to_euler_roundtrip()
call test_euler_to_euler_roundtrip()
call test_angvel_to_qdot_roundtrip()

contains

!*******************************************************************************

subroutine test_quat_to_aa_roundtrip()
    real(rp) :: angle
    real(rp) :: qf(4), qb(4), axis(3)

    write(*,*) "Testing quaternion -> axis angle (roundtrip)"
    !Generate a random quaternion
    call quat_rand(qf)
    write(*,'(a,1x,4(g0.8,1x))') " Quaternion", qf
    !Convert qf to axis angle
    call quat_to_aa(qf, axis, angle)
    write(*,'(a,1x,3(g0.8,1x),1x,a,1x,g0.8)') " Axis", axis, "Angle", angle
    !Convert axis angle back to quaternion
    call aa_to_quat(axis, angle, qb)
    write(*,'(a,1x,4(g0.8,1x))') " Quaternion from axis angle", qb
    if ( allclose(qf, qb, 1e-8_rp, 1e-14_rp) ) then
        write(*,*) "Passed"
    else
        write(*,*) "Failed"
    end if
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

subroutine test_quat_to_dcm_roundtrip()
    real(rp) :: qf(4), qb(4), dcm(3,3)

    write(*,*)
    write(*,*) "Testing quaternion -> dcm (roundtrip)"
    !Generate a random quaternion
    call quat_rand(qf)
    write(*,'(a,1x,4(g0.8,1x))') " Quaternion", qf
    !Convert qf to dcm
    call quat_to_dcm(qf, dcm)
    write(*,*) " Direction cosine matrix"
    write(*,'(1x,3(g0.8,1x))') dcm(1,:)
    write(*,'(1x,3(g0.8,1x))') dcm(2,:)
    write(*,'(1x,3(g0.8,1x))') dcm(3,:)
    !Convert dcm back to quaternion
    call dcm_to_quat(dcm, qb)
    write(*,'(a,1x,4(g0.8,1x))') " Quaternion from dcm", qb
    if ( allclose(qf, qb, 1e-8_rp, 1e-14_rp) ) then
        write(*,*) "Passed"
    else
        write(*,*) "Failed"
    end if
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

subroutine test_quat_to_euler_roundtrip()
    real(rp) :: qf(4), qb(4), euler(3)
    character(3) :: seq(6)
    logical :: world(2)
    integer :: i, j

    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [.true., .false.]

    write(*,*)
    write(*,*) "Testing quaternion -> euler (roundtrip)"
    !Generate a random quaternion
    call quat_rand(qf)
    write(*,'(a,1x,4(g0.8,1x))') " Quaternion", qf
    !Convert qf to euler
    do i = 1, 6
        do j = 1,2
            call quat_to_euler(qf, seq(i), world(j), euler)
            write(*,'(a,1x,a,1x,a,1l)') " Euler angles", seq(i), &
                "world=", world(j)
            write(*,'(1x,3(g0.8,1x))') euler
            !Convert euler back to quaternion
            call euler_to_quat(euler, seq(i), world(j), qb)
            write(*,'(a,1x,4(g0.8,1x))') " Quaternion from euler", qb
            if ( allclose(qf, qb, 1e-8_rp, 1e-14_rp) ) then
                write(*,*) "Passed"
            else
                write(*,*) "Failed"
            end if
            write(*,*)
        end do
    end do
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

subroutine test_aa_to_dcm_roundtrip()
    real(rp) :: axisf(3), axisb(3), dcm(3,3)
    real(rp) :: anglef, angleb

    write(*,*)
    write(*,*) "Testing axis angle -> dcm (roundtrip)"
    !Generate a random axis angle
    call aa_rand(axisf, anglef)
    write(*,'(a,1x,3(g0.8,1x),1x,a,1x,g0.8)') " Axis", axisf, "Angle", anglef
    !Convert axis angle to dcm
    call aa_to_dcm(axisf, anglef, dcm)
    write(*,*) " Direction cosine matrix"
    write(*,'(1x,3(g0.8,1x))') dcm(1,:)
    write(*,'(1x,3(g0.8,1x))') dcm(2,:)
    write(*,'(1x,3(g0.8,1x))') dcm(3,:)
    !Convert dcm back to axis angle
    call dcm_to_aa(dcm, axisb, angleb)
    write(*,'(a,1x,3(g0.8,1x),1x,a,1x,g0.8)') " Axis", axisb, "Angle", angleb
    if ( allclose(axisf, axisb, 1e-8_rp, 1e-14_rp) .and. &
        isclose(anglef, angleb, 1e-8_rp, 1e-14_rp) ) then
        write(*,*) "Passed"
    else
        write(*,*) "Failed"
    end if
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

subroutine test_aa_to_euler_roundtrip()
    real(rp) :: axisf(3), axisb(3), euler(3)
    real(rp) :: anglef, angleb
    character(3) :: seq(6)
    logical :: world(2)
    integer :: i, j

    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [.true., .false.]

    write(*,*)
    write(*,*) "Testing axis angle -> euler (roundtrip)"
    !Generate a random axis angle
    call aa_rand(axisf, anglef)
    write(*,'(a,1x,3(g0.8,1x),1x,a,1x,g0.8)') " Axis", axisf, "Angle", anglef
    !Convert axis angle to euler
    do i = 1, 6
        do j = 1,2
            call aa_to_euler(axisf, anglef, seq(i), world(j), euler)
            write(*,'(a,1x,a,1x,a,1l)') " Euler angles", seq(i), &
                "world=", world(j)
            write(*,'(1x,3(g0.8,1x))') euler
            !Convert euler back to axis angle
            call euler_to_aa(euler, seq(i), world(j), axisb, angleb)
            write(*,'(a,1x,3(g0.8,1x),1x,a,1x,g0.8)') " Axis", axisb,&
                "Angle", angleb
            if ( allclose(axisf, axisb, 1e-8_rp, 1e-14_rp) .and. &
                isclose(anglef, angleb, 1e-8_rp, 1e-14_rp) ) then
                write(*,*) "Passed"
            else
                write(*,*) "Failed"
            end if
            write(*,*)
        end do
    end do
    write(*,*) "----------------------------------------------------------"
    
    end subroutine

!*******************************************************************************

subroutine test_dcm_to_euler_roundtrip()
    real(rp) :: dcmf(3,3), dcmb(3,3), euler(3), axis(3)
    real(rp) :: angle
    character(3) :: seq(6)
    logical :: world(2)
    integer :: i, j

    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [.true., .false.]

    write(*,*)
    write(*,*) "Testing dcm -> euler (roundtrip)"
    !Generate a random dcm
    call aa_rand(axis, angle)
    call aa_to_dcm(axis, angle, dcmf)
    write(*,*) " Direction cosine matrix"
    write(*,'(1x,3(g0.8,1x))') dcmf(1,:)
    write(*,'(1x,3(g0.8,1x))') dcmf(2,:)
    write(*,'(1x,3(g0.8,1x))') dcmf(3,:)
    !Convert dcm to euler
    do i = 1, 6
        do j = 1,2
            call dcm_to_euler(dcmf, seq(i), world(j), euler)
            write(*,'(a,1x,a,1x,a,1l)') " Euler angles", seq(i), &
                "world=", world(j)
            write(*,'(1x,3(g0.8,1x))') euler
            !Convert euler back to dcm
            call euler_to_dcm(euler, seq(i), world(j), dcmb)
            write(*,*) " Direction cosine matrix"
            write(*,'(1x,3(g0.8,1x))') dcmb(1,:)
            write(*,'(1x,3(g0.8,1x))') dcmb(2,:)
            write(*,'(1x,3(g0.8,1x))') dcmb(3,:)
            if ( allclose(dcmf, dcmb, 1e-8_rp, 1e-14_rp) ) then
                write(*,*) "Passed"
            else
                write(*,*) "Failed"
            end if
            write(*,*)
        end do
    end do
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

subroutine test_euler_to_euler_roundtrip()
    real(rp) :: eulerf(3), eulerb(3), euler(3), axis(3)
    real(rp) :: angle
    character(3) :: seq(6)
    logical :: world(2) 
    integer :: i, j, ii, jj

    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [.true., .false.]

    write(*,*)
    write(*,*) "Testing euler -> euler (roundtrip)"
    !Generate a random axis angle and convert to euler
    call aa_rand(axis, angle)

    do i = 1, 5
        do j = 1, 1
            call aa_to_euler(axis, angle, seq(i), world(j), eulerf)
            write(*,'(a,1x,a,1x,a,1l)') " Euler", seq(i), "world=", world(j)
            write(*,'(1x,3(g0.8,1x))') eulerf
            write(*,*)
            !Convert one euler sequence to another
            do ii = i, 6
                do jj = j, 2
                    call euler_to_euler(eulerf, seq(i), world(j), &
                        seq(ii), world(jj), euler)
                    write(*,'(a,1x,a,1x,a,1l)') " Euler", seq(ii), &
                        "world=", world(jj)
                    write(*,'(1x,3(g0.8,1x))') euler
                    !Convert the euler sequence back 
                    call euler_to_euler(euler, seq(ii), world(jj), &
                        seq(i), world(j), eulerb)
                    write(*,'(a,1x,a,1x,a,1l)') " Euler", seq(i), &
                        "world=", world(j)
                    write(*,'(1x,3(g0.8,1x))') eulerb
                    if ( allclose(eulerf, eulerb, 1e-8_rp, 1e-14_rp) ) then
                        write(*,*) "Passed"
                    else
                        write(*,*) "Failed"
                    end if
                    write(*,*)
                end do
            end do
        end do
    end do
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

subroutine test_angvel_to_qdot_roundtrip()
    !Inconvenient to do qdot -> angvel round trip as then for a given q, qdot
    !has to be chosen such that it is orthogonal to q. qdot does not have to be
    !a unit vector.

    real(rp) :: q(4), qdot(4), omegaf(3), omegab(3)

    write(*,*) "Testing angvel -> qdot (roundtrip)"
    !Generate a random quaternion
    call quat_rand(q)
    !Generate a random angular velocity between -1 and 1
    call random_number(omegaf)
    omegaf = 2.0_rp*omegaf - 1.0_rp
    write(*,'(a,1x,3(g0.8,1x))') " Omega", omegaf
    !Convert omegaf to qdot
    call quat_deriv_from_angvel(q, omegaf, qdot)
    write(*,'(a,1x,4(g0.8,1x))') " qdot", qdot
    !Convert qdot back to angular velocity
    call quat_deriv_to_angvel(q, qdot, omegab)
    write(*,'(a,1x,4(g0.8,1x))') " Omega from qdot", omegab
    if ( allclose(omegaf, omegab, 1e-8_rp, 1e-14_rp) ) then
        write(*,*) "Passed"
    else
        write(*,*) "Failed"
    end if
    write(*,*) "----------------------------------------------------------"

    end subroutine

!*******************************************************************************

end program
