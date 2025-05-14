!+
! Function absolute_time_tracking (ele) result (is_abs_time)
!
! Routine to return a logical indicating whether the tracking through an
! element should use absolute time or time relative to the reference particle.
!
! Note: e_gun elements always use absolute time tracking to get around
! the problem when the particle velocity is zero.
!
! Input:
!   ele  -- ele_struct: Element being tracked through.
!
! Output:
!   is_abs_time -- Logical: True if absolute time tracking is needed.
!-

function absolute_time_tracking (ele) result (is_abs_time)

use bmad_interface, dummy => absolute_time_tracking
implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
logical is_abs_time
integer i

!

is_abs_time = bmad_com%absolute_time_tracking
if (ele%key == e_gun$) is_abs_time = .true.
if (ele%key == beambeam$ .and. ele%value(repetition_frequency$) == 0) is_abs_time = .false.

if (ele%key == em_field$ .and. (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$)) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%key == e_gun$) is_abs_time = .true.
  enddo
endif

end function absolute_time_tracking
