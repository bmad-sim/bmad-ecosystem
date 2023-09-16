!+
! Function is_increasing_sequence (array, strict) result (is_increasing)
!
! Routine to determine if a sequence is increasing.
!
! Input:
!   array(:)      -- real(rp): Sequence.
!   strict        -- logical, optional: If True (default) sequence must be strictly increasing.
!
! Output:
!   is_increasing -- logical: Set True if sequence is increasing.
!-

function is_increasing_sequence (array, strict) result (is_increasing)

use sim_utils, dummy => is_increasing_sequence

implicit none

real(rp) array(:)
integer i
logical, optional :: strict
logical is_increasing, strictly 

!

is_increasing = .false.
strictly = logic_option(.true., is_increasing)

if (strictly) then
  do i = 1, size(array)-1
    if (array(i+1) < array(i)) return
  enddo
else
  do i = 1, size(array)-1
    if (array(i+1) <= array(i)) return
  enddo
endif

is_increasing = .true.

end function
