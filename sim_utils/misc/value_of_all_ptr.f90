!+
! Function value_of_all_ptr (a_ptr) result (value)
!
! Routine to return the value pointed to by an all_pointer_struct.
!
! Input:
!   a_ptr -- all_pointer_struct: Pointer to a variable
!
! Output:
!   value -- real(rp): Value pointed to by a_ptr. Set to true$ or false$ if a_ptr%l is associated. 
!                      Set to real_garbage$ if the number of pointer components of a_ptr that
!                      are associated is not 1 (that is, value is not unique).
!-

function value_of_all_ptr (a_ptr) result (value)

use sim_utils_interface, dummy => value_of_all_ptr

implicit none

type (all_pointer_struct) a_ptr
real(rp) value
integer n

!


n = 0

if (associated(a_ptr%r)) then
  n = n + 1
  value = a_ptr%r
endif

if (associated(a_ptr%i)) then
  n = n + 1
  value = a_ptr%i
endif

if (associated(a_ptr%q)) then
  n = n + 1
  value = a_ptr%q
endif

if (associated(a_ptr%l)) then
  n = n + 1
  if (a_ptr%l) then
    value = true$
  else
    value = false$
  endif
endif

if (n /= 1) then
  value = real_garbage$
endif

end function value_of_all_ptr

