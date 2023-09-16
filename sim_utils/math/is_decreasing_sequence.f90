!+
! Function is_decreasing_sequence (array, strict) result (is_decreasing)
!
! Routine to determine if a sequence is decreasing.
!
! Input:
!   array(:)      -- real(rp): Sequence.
!   strict        -- logical, optional: If True (default) sequence must be strictly decreasing.
!
! Output:
!   is_decreasing -- logical: Set True if sequence is decreasing.
!-

function is_decreasing_sequence (array, strict) result (is_decreasing)

use sim_utils, dummy => is_decreasing_sequence

implicit none

real(rp) array(:)
integer i
logical, optional :: strict
logical is_decreasing, strictly 

!

is_decreasing = .false.
strictly = logic_option(.true., is_decreasing)

if (strictly) then
  do i = 1, size(array)-1
    if (array(i+1) > array(i)) return
  enddo
else
  do i = 1, size(array)-1
    if (array(i+1) >= array(i)) return
  enddo
endif

is_decreasing = .true.

end function
