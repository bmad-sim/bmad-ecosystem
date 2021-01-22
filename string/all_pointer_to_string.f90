!+
! Function all_pointer_to_string (a_ptr, err) result (str)
! 
! Routine to turn the value pointed to by an all_pointer_struct into a string for printing.
!
! Input:
!   a_ptr -- all_pointer_struct:
!
! Output:
!   str   -- character(24): String representation. If zero or several components
!             of a_ptr are associated then str will be set to blank ''.
!   err   -- logical, optional: Set True if none of the pointers of a_ptr are associated.
!-

function all_pointer_to_string (a_ptr, err) result (str)

use sim_utils_struct

implicit none

type (all_pointer_struct) a_ptr
integer n
character(24) str
logical, optional :: err

!

n = 0
if (associated(a_ptr%r)) n = n + 1
if (associated(a_ptr%i)) n = n + 1
if (associated(a_ptr%l)) n = n + 1

if (n /= 1) then
  str = null_name$
  if (present(err)) err = .true.
  return
endif

if (associated(a_ptr%r)) then
  write (str, *) a_ptr%r
elseif (associated(a_ptr%i)) then
  write (str, *) a_ptr%i
elseif (associated(a_ptr%l)) then
  write (str, *) a_ptr%l
endif

if (present(err)) err = .false.

end function
