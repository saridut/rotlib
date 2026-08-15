module m_precision

use, intrinsic ::  iso_fortran_env, only: int32, int64, real64

implicit none

!> Default integer precision
integer, parameter :: ip = int32
integer, parameter :: ip_long = int64

!> Default real precision
integer, parameter :: rp = real64

end module m_precision
