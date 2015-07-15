module precision_def

integer, parameter :: rdef = selected_real_kind(11)
integer, parameter :: rdef2 = 2*rdef    ! used for nr.
integer, parameter :: rp = rdef

integer, parameter :: sp = kind(1e0)
integer, parameter :: dp = selected_real_kind(2*precision(1e0_sp))

integer, parameter :: i4_b = selected_int_kind(9)  ! Equiv to NR I4B

type global_common_struct
  logical :: be_thread_safe = .false.    ! Avoid thread unsafe practices?
  logical :: exit_on_error  = .true.     ! Exit program on error?
end type

type (global_common_struct), save :: global_com

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
