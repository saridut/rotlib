!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta                                                 !
!                                                                                !
! Permission is hereby granted, free of charge, to any person obtaining a copy   !
! of this software and associated documentation files (the "Software"), to deal  !
! in the Software without restriction, including without limitation the rights   !
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      !
! copies of the Software, and to permit persons to whom the Software is          !
! furnished to do so, subject to the following conditions:                       !
!                                                                                !
! The above copyright notice and this permission notice shall be included in all !
! copies or substantial portions of the Software.                                !
!                                                                                !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     !
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       !
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    !
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         !
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  !
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  !
! SOFTWARE.                                                                      !
!                                                                                !
!********************************************************************************!

module utils_math_m
    !! Various utility functions, particularly for use with small matrices.

use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
use constants_m

implicit none

interface allclose
    !! Checks if two arrays are elementwise close within tolerance
    module procedure allclose_rank1
    module procedure allclose_rank2
    module procedure allclose_rank3
end interface 

interface swap
    !! Swaps two arrays
    module procedure swap_integer
    module procedure swap_real
    module procedure swap_complex
end interface swap

contains

!******************************************************************************

elemental subroutine swap_integer(a, b)
    integer, intent(in out) :: a
    integer, intent(in out) :: b
    integer :: temp

    temp = a
    a = b
    b = temp
    end subroutine

elemental subroutine swap_real(a, b)
    real(rp), intent(in out) :: a
    real(rp), intent(in out) :: b
    real(rp) :: temp

    temp = a
    a = b
    b = temp
    end subroutine

elemental subroutine swap_complex(a, b)
    complex(rp), intent(in out) :: a
    complex(rp), intent(in out) :: b
    complex(rp) :: temp

    temp = a
    a = b
    b = temp
    end subroutine

!******************************************************************************
 
elemental logical function isclose(a, b, rel_tol, abs_tol)
    !! Checks if two floating point numbers of type double are close within
    !! tolerance.
    !!
    !! Based on python implementation at
    !! <https://github.com/PythonCHB/close_pep/blob/master/is_close.py>.
    !! The `method='weak'` option is used here.
    !""

    real(rp), intent(in)    :: a
    real(rp), intent(in)    :: b
    real(rp), intent(in), optional    :: rel_tol
        !! Relative tolerance, `rel_tol >= 0`, default 1e-10
    real(rp), intent(in), optional    :: abs_tol
        !! Absolute tolerance, `abs_tol >= 0`, default 0.0
    real(rp) :: rel_tol_
    real(rp) :: abs_tol_
    real(rp) :: diff

    rel_tol_ = 1e-10_rp
    abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol
    
    if (a == b) then  ! short-circuit exact equality
        isclose = .true.
        return
    end if
    
    if ((.not. ieee_is_finite(a)) .or. (.not. ieee_is_finite(b))) then
         ! Includes the case of two infinities of opposite sign, or
         ! one infinity and one finite number. Two infinities of opposite sign
         ! would otherwise have an infinite relative tolerance.
         isclose = .false.
        return
    end if
    
    diff = abs(b - a)
    isclose = ( ((diff <= abs(rel_tol_*b)) .or. (diff <= abs(rel_tol_*a))) &
            .or. (diff <= abs_tol_) )

    end function

!******************************************************************************

logical function allclose_rank1(a, b, rel_tol, abs_tol)
    !! Checks if two rank-1 floating point arrays of type double are close within
    !! tolerance.
    !""
    real(rp), dimension(:), intent(in) :: a
        !! `(m,)` array
    real(rp), dimension(:), intent(in) :: b
        !! `(m,)` array
    real (rp), intent(in), optional :: rel_tol
        !! Relative tolerance; default 1e-10
    real (rp), intent(in), optional :: abs_tol
        !! Absolute tolerance; default 0.0
    real (rp) :: abs_tol_, rel_tol_ 

    rel_tol_ = 1e-10_rp
    abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol

    allclose_rank1 = all(isclose(a, b, rel_tol_, abs_tol_))

    end function


logical function allclose_rank2(a, b, rel_tol, abs_tol)
    !! Checks if two rank-2 floating point arrays of type double are close within
    !! tolerance.
    !""
    real(rp), dimension(:,:), intent(in) :: a
        !! `(m,n)` array
    real(rp), dimension(:,:), intent(in) :: b
        !! `(m,n)` array
    real (rp), intent(in), optional :: rel_tol
        !! Relative tolerance; default 1e-10
    real (rp), intent(in), optional :: abs_tol
        !! Absolute tolerance; default 0.0
    real (rp) :: abs_tol_, rel_tol_ 

    rel_tol_ = 1e-10_rp
    abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol

    allclose_rank2 = all(isclose(a, b, rel_tol_, abs_tol_))

    end function


logical function allclose_rank3(a, b, rel_tol, abs_tol)
    !! Checks if two rank-3 floating point arrays of type double are close within
    !! tolerance.
    !""
    real(rp), dimension(:,:,:), intent(in) :: a
        !! `(m,n,p)` array
    real(rp), dimension(:,:,:), intent(in) :: b
        !! `(m,n,p)` array
    real (rp), intent(in), optional :: rel_tol
        !! Relative tolerance; default 1e-10
    real (rp), intent(in), optional :: abs_tol
        !! Absolute tolerance; default 0.0
    real (rp) :: abs_tol_, rel_tol_ 

    rel_tol_ = 1e-10_rp
    abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol

    allclose_rank3 = all(isclose(a, b, rel_tol_, abs_tol_))

    end function

!************************************************************************

elemental subroutine rad2deg(rad, deg)

    real(rp), intent(in) :: rad
    real(rp), intent(out) :: deg

    deg = (180.0_rp/math_pi)*rad

    end subroutine

!************************************************************************

elemental subroutine deg2rad(deg, rad)

    real(rp), intent(in) :: deg
    real(rp), intent(out) :: rad

    rad = (math_pi/180.0_rp)*deg

    end subroutine

!******************************************************************************

elemental function sinc(x) result(res)
    !! Calculates sin(x)/x elementwise.
    !
    ! At |x|=1e-4, the first omitted term from the Taylor expansion is 
    ! of the order of 1e-18, comfortably below double-precision rounding error
    ! near 1.
    real(rp), intent(in) :: x
    real(rp) :: res

    if (abs(x) < 1e-4_rp) then
        res = 1.0_rp - x*x/6.0_rp
    else
        res = sin(x)/x
    end if

    end function

!******************************************************************************

function norm(v) result(res)
    !! Calculates the 2-norm of a vector.
    real(rp), dimension(:), intent(in) :: v
    real(rp) :: res
    real(rp) :: s, snorm, tmp
    integer :: i

    ! scaling factor
    s = 0.0_rp
    do i = 1, size(v)
        s = max(s, abs(v(i)))
    end do

    if (s == 0.0_rp) then
        res = 0.0_rp ! norm of a zero vector is zero
        return
    end if

    ! calculate the square of the scaled norm
    snorm = 0.0_rp
    tmp = 0.0_rp
    do i = 1, size(v)
        tmp = v(i)/s
        snorm = snorm + tmp*tmp
    end do

    ! unscale and return
    res = s*sqrt(snorm)

    end function

!******************************************************************************

subroutine normalized(v, estat)
    !! Normalizes a vector in-place. If the magnitude of the vector is zero, 
    !! the vector is not modified and an optional error status of `1` is returned.
    real(rp), dimension(:), intent(in out) :: v
        !! A vector
    integer, intent(out), optional :: estat
        !! Error status. Equals `1` in case of a zero vector, and `0` otherwise.
    real(rp) :: s, snorm
    integer :: i

    ! scaling factor
    s = 0.0_rp
    do i = 1, size(v)
        s = max(s, abs(v(i)))
    end do

    if (s == 0.0_rp) then
        ! cannot normalize a zero vector
        if (present(estat)) estat = 1
        return
    end if

    v = v/s ! scale the vector
    snorm = sqrt(dot_product(v, v)) ! scaled norm
    v = v/snorm ! Divide by the scaled norm
    if (present(estat)) estat = 0 ! Set optional error status

    end subroutine

!************************************************************************

subroutine cross(a, b, c)
    !! Calculates the cross product between two 3-element vectors.
    real(rp), dimension(3), intent(in)  :: a
    real(rp), dimension(3), intent(in)  :: b
    real(rp), dimension(3), intent(out) :: c
        !! Cross product of `a` and `b`; **c** = **a** x **b**
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    end subroutine

!******************************************************************************

subroutine cross_mat(a, mat)
    !! Calculates the cross product matrix of a 3-element vector. The cross
    !! product matrix **A** of **a** is defined as **a** x **b** = **A** . **b**,
    !! where **b** is another 3-element vector.
    real(rp), dimension(3), intent(in)    :: a
    real(rp), dimension(3,3), intent(out) :: mat
        !! Cross product matrix of `a`
    
    mat = 0.0_rp
    mat(2,1) =  a(3)
    mat(3,1) = -a(2)
    mat(1,2) = -a(3)
    mat(3,2) =  a(1)
    mat(1,3) =  a(2)
    mat(2,3) = -a(1)

    end subroutine

!******************************************************************************

function vec_angle_3d (a, b) result(res)
    !! Returns the angle (in radian) between two non-zero vectors in 3D.

    real(rp), dimension(3), intent(in) :: a
    real(rp), dimension(3), intent(in) :: b
    real(rp) :: res
    real(rp) :: ab, a_(3), b_(3), c(3), c_mag

    a_ = a/maxval(abs(a))
    b_ = b/maxval(abs(b))
    call cross(a_, b_, c)
    c_mag = norm(c)
    ab = dot_product(a_, b_)
    res = atan2(c_mag, ab)

    end function

!******************************************************************************

function uvec_angle (u, v) result(res)
    !! Returns the angle (in radian) between two unit vectors using Kahan's half
    !! angle formula.

    real(rp), dimension(:), intent(in) :: u
        !! Unit vector
    real(rp), dimension(:), intent(in) :: v
        !! Unit vector
    real(rp) :: res
    real(rp) :: s_mag, d_mag

    s_mag = norm(u+v)
    d_mag = norm(u-v)
    res = 2.0_rp*atan2(d_mag, s_mag)

    end function

!******************************************************************************

function scalar_triple_product (a, b, c) result(res)
    !! Returns the scalar triple product **a**.(**b** x **c**)
    real(rp), dimension(3), intent(in)  :: a
    real(rp), dimension(3), intent(in)  :: b
    real(rp), dimension(3), intent(in)  :: c
    real(rp)  :: res

    res = a(1)*( b(2)*c(3)-c(2)*b(3) ) - a(2)*( b(1)*c(3)-c(1)*b(3) ) &
            + a(3)*( b(1)*c(2)-c(1)*b(2) )

    end function

!******************************************************************************

subroutine vector_triple_product (a, b, c, d)
    !! Returns the vector triple product **d** = **a** x (**b** x **c**)
    real(rp), dimension(3), intent(in)  :: a
    real(rp), dimension(3), intent(in)  :: b
    real(rp), dimension(3), intent(in)  :: c
    real(rp), dimension(3), intent(out) :: d
     !! Vector triple product

    d = b*( a(1)*c(1) + a(2)*c(2) + a(3)*c(3) )  &
        - c*( a(1)*b(1) + a(2)*b(2) + a(3)*b(3) )

    end subroutine

!******************************************************************************

subroutine outer(a, b, c)
    !! Calculates the outer product of two vectors, \(c_{ij} = a_i  b_j\).
    real(rp), dimension(:), intent(in)    :: a
        !! `(m,)` array
    real(rp), dimension(:), intent(in)    :: b
        !! `(n,)` array
    real(rp), dimension(:,:), intent(out) :: c
        !! `(m,n)` array; Outer product
    integer :: m, n
    integer :: i, j

    m = size(a)
    n = size(b)

    do j = 1, n
        do i = 1, m
            c(i,j) = a(i)*b(j)
        end do
    end do

    end subroutine

!******************************************************************************

subroutine linspace(start, finish, num, val, step)
    !! Generates evenly spaced numbers over a specified interval. Both end
    !! points are included. If `start` > `finish`, the returned step size (if
    !! `step` is present) will be negative.
    !!
    !! Values are interpolated directly from the endpoints to avoid overflow
    !! in the endpoint difference. Floating-point rounding can make adjacent
    !! spacings unequal or produce repeated values. The supplied endpoints
    !! are assigned exactly.
    !""
    real(rp), intent(in) :: start
        !! Finite starting point
    real(rp), intent(in) :: finish
        !! Finite ending point, `finish /= start`
    integer, intent(in) :: num
        !! Number of values to generate, `num >= 2`
    real(rp), dimension(:), intent(in out) :: val
        !! `(n,)` array, where `n >= num`. Generated values. Elements beyond
        !! the first `num` values are left as is.
    real(rp), intent(out), optional :: step
        !! Nominal step size, which may differ from adjacent output spacings
        !! because of rounding. For `num=2`, its calculation can overflow if
        !! the endpoint difference is not representable. No step is calculated
        !! when this argument is absent.
    real(rp) :: intervals, t, difference
    logical :: opposite_signs
    integer :: i

    if (num < 2) then
        write(*,*) '[linspace] `num` must be >= 2'
        stop
    end if

    if (num > size(val)) then
        write(*,*) '[linspace] `size(val)` must be >= `num`'
        stop
    end if

    if (start == finish) then
        write(*,*) '[linspace] `start` must not be equal to `finish`'
        stop
    end if

    intervals = real(num-1, rp)
    opposite_signs = (start < 0.0_rp .and. finish > 0.0_rp) .or. &
                     (start > 0.0_rp .and. finish < 0.0_rp)

    val(1) = start; val(num) = finish

    if (opposite_signs) then
        if (present(step)) then
            if (num == 2) then
                step = finish-start
            else
                step = finish/intervals - start/intervals
            end if
        end if
        do i = 2, num-1
            t = real(i-1, rp)/intervals
            val(i) = (1.0_rp-t)*start + t*finish
        end do
    else
        difference = finish-start
        if (present(step)) step = difference/intervals
        do i = 2, num-1
            t = real(i-1, rp)/intervals
            val(i) = start + t*difference
        end do
    end if

    end subroutine

!******************************************************************************

subroutine logspace(start, finish, num, val, base)
    !! Generates numbers spaced evenly on a log scale.
    !! In linear space, the sequence starts at `base ** start`
    !! (`base` to the power of `start`) and ends with `base ** finish`
    !!
    !! Exponents are generated by `linspace`. Floating-point rounding can
    !! make adjacent ratios unequal or produce repeated values. The endpoint
    !! powers are subject to rounding, overflow, and underflow as well.
    !""
    real(rp), intent(in) :: start
        !! Finite starting exponent; `base ** start` is the starting value
    real(rp), intent(in) :: finish
        !! Finite ending exponent, `finish /= start`; `base ** finish`
        !! is the ending value
    integer, intent(in) :: num
        !! Number of values to generate, `num >= 2`
    real(rp), dimension(:), intent(in out) :: val
        !! `(n,)` array, where `n >= num`. Generated values. Elements beyond
        !! the first `num` values are left as is.
    real(rp), intent(in), optional :: base
        !! Base of the logspace. Must be finite and positive, default is 10.
    real(rp) :: base_
    integer :: i

    if (num < 2) then
        write(*,*) 'logspace `num` must be >= 2'
        stop
    end if

    if (num > size(val)) then
        write(*,*) 'logspace `size(val)` must be >= `num`'
        stop
    end if

    if (start == finish) then
        write(*,*) 'logspace `start` must not be equal to `finish`'
        stop
    end if

    base_ = 10.0_rp
    if (present(base)) base_ = base

    call linspace(start, finish, num, val)

    do i = 1, num
        val(i) = base_**val(i)
    end do

    end subroutine

!******************************************************************************

function det(A) result(res)
    !! Returns the determinant of an (N x N) matrix, where N = 2, 3, or 4.
    !!
    !! Original routine by [David Simpson](http://www.davidgsimpson.com/software.html)
    !! @note For a general NxN matrix do an LU decomp
    !""
    real(rp), dimension(:,:), intent(in) :: A
        !! `(N,N)` array, where N = 2, 3, or 4.
    real(rp) :: res
    integer :: nrows

    nrows = size(A, 1)

    if (nrows==2) then
        res = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    else if (nrows==3) then
        res = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    else if (nrows==4) then
        res = A(1,1) * ( A(2,2) * (A(3,3)*A(4,4) - A(3,4)*A(4,3))    &
                       + A(2,3) * (A(3,4)*A(4,2) - A(3,2)*A(4,4))    &
                       + A(2,4) * (A(3,2)*A(4,3) - A(3,3)*A(4,2)) )  &
            - A(1,2) * ( A(2,1) * (A(3,3)*A(4,4) - A(3,4)*A(4,3))    &
                       + A(2,3) * (A(3,4)*A(4,1) - A(3,1)*A(4,4))    &       
                       + A(2,4) * (A(3,1)*A(4,3) - A(3,3)*A(4,1)) )  &
            + A(1,3) * ( A(2,1) * (A(3,2)*A(4,4) - A(3,4)*A(4,2))    &
                       + A(2,2) * (A(3,4)*A(4,1) - A(3,1)*A(4,4))    &
                       + A(2,4) * (A(3,1)*A(4,2) - A(3,2)*A(4,1)) )  &
            - A(1,4) * ( A(2,1) * (A(3,2)*A(4,3) - A(3,3)*A(4,2))    &       
                       + A(2,2) * (A(3,3)*A(4,1) - A(3,1)*A(4,3))    &
                       + A(2,3) * (A(3,1)*A(4,2) - A(3,2)*A(4,1)) )

     else
         error stop 'Matrix rank must be 2, 3, or 4.'
     end if

    end function

!******************************************************************************

function trace(mat) result (res)
    !! Returns the trace of a square matrix

    real(rp), intent(in) :: mat(:,:)
        !! `(N,N)` array
    real(rp) :: res
    integer :: nrows
    integer :: i

    nrows = size(mat, 1)
    res = 0.0_rp

    do i = 1, nrows
        res = res + mat(i,i)
    end do

    end function

!******************************************************************************

subroutine identity(mat_eye)
    !! Creates an identity matrix of size n x n.
    real(rp), dimension(:,:), intent(out) :: mat_eye
        !! `(n,n)` array
    integer :: nrows
    integer :: i

    nrows = size(mat_eye, 1)

    mat_eye = 0.0_rp

    do i = 1, nrows
        mat_eye(i,i) = 1.0_rp
    end do
    
    end subroutine

!******************************************************************************

subroutine get_diagonal(mat, d)
    !! Returns the diagonal elements of a square matrix.
    real(rp), dimension(:,:), intent(in) :: mat
        !! `(n,n)` array
    real(rp), dimension(:), intent(out) :: d
        !! `(n,)` array; contains the entries of the main diagonal
    integer :: n
    integer :: i
    
    n = size(mat, 1)

    do i = 1, n
        d(i) = mat(i,i)
    end do

    end subroutine

!******************************************************************************

subroutine add_transpose(mat)
    !! Adds a square matrix and its transpose in place: \(A_{ij} = A_{ij } + A_{ji}\)
    real(rp), dimension(:,:), intent(in out) :: mat
        !! `(n,n)` array
    integer :: n
    integer :: i
    integer :: j

    n = size(mat,1)

    !Dealing with the upper triangular part (including diagonal)
    do j = 1, n
        do i = 1, j
            mat(i,j) = mat(i,j) + mat(j,i)
        end do
    end do

    !Dealing with the strictly lower triangular part (i.e. excluding diagonal)
    do j = 1, (n-1)
        do i = (j+1), n
            mat(i,j) = mat(j,i)
        end do
    end do

    end subroutine

!******************************************************************************

subroutine subtract_transpose(mat)
    !! Calculates the difference of a square matrix and its transpose in place:
    !! \(A_{ij} = A_{ij } - A_{ji}\)
    real(rp), dimension(:,:), intent(in out) :: mat
    integer :: n
    integer :: i
    integer :: j

    n = size(mat,1)

    !Dealing with the strictly upper triangular part (i.e. excluding diagonal)
    do j = 2, n
        do i = 1, (j-1)
            mat(i,j) = mat(i,j) - mat(j,i)
        end do
    end do

    !Dealing with the strictly lower triangular part (i.e. excluding diagonal)
    do j = 1, (n-1)
        do i = (j+1), n
            mat(i,j) = -mat(j,i)
        end do
    end do

    !Putting the diagonal to zero
    do i = 1, n
        mat(i,i) = 0.0_rp
    end do

    end subroutine

!******************************************************************************

subroutine multiply_transpose(A, B)
    !! Multiplies a matrix with its transpose:
    !! \(\mathbf{\mathrm{B}} = \mathbf{\mathrm{A}} \cdot \mathbf{\mathrm{A}}^T\)
    !""
    real(rp), dimension(:,:), intent(in) :: A
        !! `(m,n)` array
    real(rp), dimension(:,:), intent(out) :: B
        !! `(m,m)` array
        !A is m x n, A^T is n x m, so B is m x m

    B = matmul(A, transpose(A))

    end subroutine

!******************************************************************************

function get_quad_form(A, x) result(res)
    !! Calculates the quadratic form *x^T A x*, where A is an *n x n* matrix and *x* is a
    !! vector of length *n*
    !""
    real(rp), dimension(:,:), intent(in) :: A
        !! `(n,n)` array
    real(rp), dimension(:), intent(in) :: x
        !! `(n,)` array
    real(rp) :: res
    integer :: n
    integer :: i, j

    n = size(A,2)
    res = 0.0_rp

    do j = 1, n
        do i = 1, n
            res = res + A(i,j)*x(j)*x(i)
        end do
    end do

    end function

!******************************************************************************

subroutine orth_mgs(a, estat)
    !! Orthogonalizes a set of non-zero vectors in-place using modified
    !! Gram-Schmidt orthonormalization.
    !! 
    !! *Reference:* Golub and Van Loan, Matrix Computations, 4th edition,
    !! Section 5.2.8, Algorithm 5.2.6, p. 255.

    real(rp), dimension(:,:), intent(in out) :: a
        !! `(m,n)` rank `n` array, where `m >= n`, storing the vectors
        !! column-wise. On return, the columns of the matrix are overwritten
        !! with the orthogonal vectors.
    integer, intent(out) :: estat
        !! Error code. Returns `1` if during execution a column is found to be
        !! zero, `0` otherwise. A return with `estat=1` may leave `a` partially
        !! overwritten.
    integer :: j, k, n

    n = size(a, 2)
    do k = 1, n
        call normalized(a(:,k), estat)
        if (estat == 1) return
        do j = k+1, n
            a(:,j) = a(:,j) - dot_product(a(:,k),a(:,j)) * a(:,k)
        end do
    end do

    end subroutine

!******************************************************************************

subroutine invert_mat33(a, inv_a)
    !! Inverts a 3x3 matrix with non-zero determinant.
    !!
    !! *Reference:* https://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf

    real(rp), dimension(3,3),  intent(in) :: a
        !! Matrix
    real(rp), dimension(3,3), intent(out) :: inv_a
        !! Matrix inverse
    real(rp) :: as(3,3), det_as, s

    s = maxval(abs(a)) ! scale
    as = a/s ! scaled a

    inv_a(1,1) =   as(2,2)*as(3,3) - as(2,3)*as(3,2)
    inv_a(2,1) = -(as(2,1)*as(3,3) - as(2,3)*as(3,1))
    inv_a(3,1) =   as(2,1)*as(3,2) - as(2,2)*as(3,1)

    inv_a(1,2) = -(as(1,2)*as(3,3) - as(1,3)*as(3,2))
    inv_a(2,2) =   as(1,1)*as(3,3) - as(1,3)*as(3,1)
    inv_a(3,2) = -(as(1,1)*as(3,2) - as(1,2)*as(3,1))

    inv_a(1,3) =   as(1,2)*as(2,3) - as(1,3)*as(2,2)
    inv_a(2,3) = -(as(1,1)*as(2,3) - as(1,3)*as(2,1))
    inv_a(3,3) =   as(1,1)*as(2,2) - as(1,2)*as(2,1)

    !Determinant of the scaled matrix
    det_as = as(1,1)*inv_a(1,1) + as(1,2)*inv_a(2,1) + as(1,3)*inv_a(3,1)
    inv_a = inv_a/det_as
    
    ! Unscale
    inv_a = inv_a/s

    end subroutine

!******************************************************************************

subroutine eigval_33rsym(a, ev)
    !! Calculates the eigenvalues of a finite real symmetric 3 x 3 matrix,
    !! in decreasing order. The input matrix is not modified.
    !!
    !! Uses maximum-entry scaling and Jacobi rotations to avoid overflow in
    !! intermediate products and cancellation in cubic-equation formulas near
    !! repeated eigenvalues. Diagonal matrices are handled directly.
    !! Iteration stops when the largest off-diagonal entry is at most machine
    !! epsilon times the largest absolute diagonal entry of the scaled matrix.
    !!
    !! Accuracy is relative to the matrix scale; high relative accuracy for
    !! very small eigenvalues or very small eigenvalue gaps is not guaranteed.
    !! Rescaling can overflow if an eigenvalue is outside the range of `rp`,
    !! or underflow if it is too small to represent.
    !!
    !! *Reference:* J. Kopp, Efficient numerical diagonalization of Hermitian
    !! 3 x 3 matrices, https://arxiv.org/abs/physics/0610206

    real(rp), dimension(3,3), intent(in) :: a
        !! Finite symmetric input matrix. Both triangles must be supplied.
    real(rp), dimension(3), intent(out) :: ev
        !! Eigenvalues, with `ev(1) >= ev(2) >= ev(3)`.
    real(rp) :: b(3,3), scale, off, diagonal_scale
    real(rp) :: delta, tangent, cosine, sine, x, y
    integer :: p, q, k

    if (a(1,2) == 0.0_rp .and. a(1,3) == 0.0_rp .and. &
        a(2,3) == 0.0_rp) then
        ! Preserve diagonal entries exactly, including the zero matrix.
        ev = [a(1,1), a(2,2), a(3,3)]
    else
        scale = maxval(abs(a))
        b = a/scale
        do
            ! Annihilate the largest off-diagonal entry.
            p = 1; q = 2
            if (abs(b(1,3)) > abs(b(p,q))) then
                p = 1; q = 3
            end if
            if (abs(b(2,3)) > abs(b(p,q))) then
                p = 2; q = 3
            end if
            off = b(p,q)
            diagonal_scale = max(abs(b(1,1)), abs(b(2,2)), abs(b(3,3)))
            if (abs(off) <= epsilon(1.0_rp)*diagonal_scale) exit

            ! Stable rotation tangent with magnitude at most one. Avoid
            ! dividing the diagonal difference by a possibly tiny off-diagonal.
            delta = 0.5_rp*(b(q,q)-b(p,p))
            tangent = off / (delta + sign(hypot(delta, off), delta))
            cosine = 1.0_rp/sqrt(1.0_rp+tangent*tangent)
            sine = tangent*cosine

            b(p,p) = b(p,p) - tangent*off
            b(q,q) = b(q,q) + tangent*off
            b(p,q) = 0.0_rp; b(q,p) = 0.0_rp
            k = 6-p-q ! The remaining row/column index.
            x = b(k,p); y = b(k,q)
            b(k,p) = cosine*x - sine*y
            b(p,k) = b(k,p)
            b(k,q) = sine*x + cosine*y
            b(q,k) = b(k,q)
        end do
        ev = [b(1,1), b(2,2), b(3,3)] * scale
    end if

    if (ev(1) < ev(2)) call swap(ev(1), ev(2))
    if (ev(1) < ev(3)) call swap(ev(1), ev(3))
    if (ev(2) < ev(3)) call swap(ev(2), ev(3))

    end subroutine

!******************************************************************************

end module utils_math_m
