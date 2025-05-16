module precision_def

integer, parameter :: rp = selected_real_kind(11)    ! "Real precision" is double precision
integer, parameter :: quadp = selected_real_kind(20) ! Quad precision
integer, parameter :: qp = quadp                     ! Define qp so that in an emergency Bmad can be compiled with qp = rp.
integer, parameter :: sp = kind(1e0)                 ! Single precision
integer, parameter :: dp = selected_real_kind(2*precision(1e0_sp))  ! Double precision
integer, parameter :: i4_b = selected_int_kind(9)  ! Equiv to NR I4B

type global_common_struct
  logical :: mp_threading_is_safe = .true.   ! Can threading be used with MP? EG ramping is not thread safe.
  logical :: exit_on_error  = .true.         ! Exit program on error?
  integer :: debug = 0                       ! Used for debugging purpeses
end type

type (global_common_struct), save :: global_com

type named_number_struct
  character(40) :: name = ''
  real(rp) :: value = 0
end type

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function rp8(int_in) result (re_out)
!
! Routine to convert from integer to real of type rp.
! This routine is used to avoid the implicit integer to single precision that happens when
! multiplying int*real(rp).
!
! Input:
!   int_in    -- integer: Input integer.
!
! Output:
!   re_out    -- real(rp): Equiv real.
!-

pure function rp8(int_in) result (re_out)
implicit none
integer, intent(in) :: int_in
real(rp) re_out
!
re_out = real(int_in, rp)
end function rp8

end module
