!+
! Function rf_ref_time_offset (ele, ds) result (time)
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
!   ds       -- real(rp), optional: Distance of particle from start edge. Default is zero.
!
! Ouput:
!   time  -- Real(rp): Offset time.
!-

function rf_ref_time_offset (ele, ds) result (time)

use bmad_routine_interface, dummy => rf_ref_time_offset
implicit none

type (ele_struct) ele
type (ele_struct), pointer :: lord
real(rp), optional :: ds
real(rp) time, beta
character(*), parameter :: r_name = 'rf_ref_time_offset'

! 

time = 0
if (.not. absolute_time_tracking(ele)) return

if (present(ds)) then
  beta = ele%value(p0c$) / ele%value(E_tot$)
  time = ds / (c_light * beta)
endif

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
  lord => pointer_to_lord (ele, 1)
  time = time + ele%value(ref_time_start$) - lord%value(ref_time_start$)
endif

end function rf_ref_time_offset
