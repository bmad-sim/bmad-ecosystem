!+
! Function assert_equal (int_arr, err_str) result (ival)
!
! Routine to check that an array of integers all have the same value.
! If this is not true, the program will be stopped.
!
! Input:
!   int_arr(:)    -- integer: Array of integers.
!
! Output:
!   err_str       -- character(*): Error message to print before stopping.
!   ival          -- integer: Value of the integers.
!-

function assert_equal (int_arr, err_str) result (ival)

use output_mod, dummy => assert_equal

implicit none

integer, intent(in) :: int_arr(:)
integer ival, i
character(*) err_str
character(*), parameter :: r_name = 'assert_equal'

!

ival = int_arr(1)

do i = 2, size(int_arr)
  if (int_arr(i) == ival) cycle
  call out_io(s_abort$, r_name, 'ARRAY SIZE MISMATCH WITH: ' // err_str, 'PLEASE REPORT THIS.')
  stop
enddo

end function
