module var_length_string_mod

use sim_utils_struct

implicit none

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!+
! Function all_equal_var_str (var_str, str1) result (all_equal)
!
! Routine to return the equivalent of all(var_str%str == str1).
!
! Input:
!   var_str(:)    -- var_length_string_struct: String array.
!   str1          -- character(*): Single string to match to.
!
! Output:
!   all_equal     -- logical: True if all var_str%str == str1.
!-

function all_equal_var_str (var_str, str1) result (all_equal)

type (var_length_string_struct) var_str(:)
character(*) str1
integer i
logical all_equal

!

all_equal = .false.

do i = 1, size(var_str)
  if (var_str(i)%str /= str1) return
enddo

all_equal = .true.

end function all_equal_var_str

end module
