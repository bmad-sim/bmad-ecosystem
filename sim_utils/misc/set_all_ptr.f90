!+
! Subroutine set_all_ptr (a_ptr, value, delta, value_set)
!
! Routine to set real or integer value pointed to by a_ptr.
! If more than one thing is pointed to, an error is generated.
!
! Input:
!   a_ptr       -- all_pointer_struct: Pointer to a variable
!   value       -- real(rp): Value to load. An error is generated if the number of pointer components 
!                      of a_ptr that are associated is not 1.
!   delta       -- real(rp): Default False. If True, load value as a change in value.
!
! Output:
!   a_ptr       -- all_pointer_struct: Value set.
!   value_set   -- real(rp): Value set. Useful when delta = True.
!-

subroutine set_all_ptr (a_ptr, value, delta, value_set)

use sim_utils, dummy => set_all_ptr

implicit none

type (all_pointer_struct) a_ptr
real(rp) value
real(rp), optional :: value_set
logical, optional :: delta
integer n
character(*), parameter :: r_name = 'set_all_ptr'


!

n = 0

if (associated(a_ptr%r)) then
  n = n + 1
  if (logic_option(.false., delta)) then
    a_ptr%r = a_ptr%r + value
  else
    a_ptr%r = value
  endif
  if (present(value_set)) value_set = a_ptr%r
endif

if (associated(a_ptr%q)) then
  n = n + 1
  if (logic_option(.false., delta)) then
    a_ptr%q = a_ptr%q + value
  else
    a_ptr%q = value
  endif
  if (present(value_set)) value_set = a_ptr%q
endif


if (associated(a_ptr%i)) then
  n = n + 1
  if (logic_option(.false., delta)) then
    a_ptr%i = a_ptr%i + nint(value)
  else
    a_ptr%i = nint(value)
  endif
  if (present(value_set)) value_set = a_ptr%i
endif

if (n == 1) return

if (n == 0) then
  call out_io(s_error$, r_name, 'POINTER NOT ASSOCIATED! PLEASE REPORT THIS!')
  return
endif

if (n > 1) then
  call out_io(s_error$, r_name, 'MULTIPLE ASSOCIATED POINTERS! PLEASE REPORT THIS!')
  return
endif

end subroutine set_all_ptr

