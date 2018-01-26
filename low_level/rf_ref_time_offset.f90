!+
! Function rf_ref_time_offset (ele) result (time)
!
! Routine to return an offset time due to the offset distance for super_slave and 
! slice_slave elemenets. This is used for calculating the RF phase
! in RF cavities. 
!
! This is only non-zero with absolute time tracking since relative time, which 
! references the time to the reference particle,  already takes this into account.
! 
! Input:
!   ele      -- ele_struct: RF Element being tracked through.
!
! Ouput:
!   time  -- Real(rp): Offset time.
!-

function rf_ref_time_offset (ele) result (time)

use bmad_routine_interface, dummy => rf_ref_time_offset
implicit none

type (ele_struct) ele
type (ele_struct), pointer :: lord
real(rp) time
character(*), parameter :: r_name = 'rf_ref_time_offset'

! 

time = 0
if (.not. absolute_time_tracking(ele)) return

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
  lord => pointer_to_lord (ele, 1)
  time = ele%value(ref_time_start$) - lord%value(ref_time_start$)
endif

end function rf_ref_time_offset
