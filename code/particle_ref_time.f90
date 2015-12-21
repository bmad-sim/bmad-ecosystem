!+
! Function particle_ref_time (orbit, ele) result (time)
!
! Routine to return the reference time used to calculate the phase of
! time-dependent EM fields.
!
! Essentially: The oscillations of EM fields are synched relative to the absolute clock if 
! absolute time tracking is used or are synched relative to the reference particle
! if relative time tracking is used. See the Bmad manual for more details.
!
! Input:
!   orbit -- Coord_struct: Particle coordinates
!   ele   -- ele_struct: Element being tracked through.
!
! Ouput:
!   time  -- Real(rp): Current time.
!-

function particle_ref_time (orbit, ele) result (time)

use multipass_mod, dummy => particle_ref_time

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (ele_pointer_struct), allocatable :: chain(:)

real(rp) time
integer ix_pass, n_links
logical abs_time

character(*), parameter :: r_name = 'particle_ref_time'

! With absolute time tracking the reference time is relative to the reference time of the element.
! This way the phase does not have to be adjusted when switching between absolute and relative time tracking.
! Note: With a multipass_slave, use the reference time of the pass element.
! Note: e_gun uses absolute time tracking to get around the problem when orbit%beta = 0.

if (absolute_time_tracking(ele)) then
    call multipass_chain(ele, ix_pass, n_links, chain)
  if (ix_pass > 1) then
    time = orbit%t - chain(1)%ele%value(ref_time_start$)
  else
    time = orbit%t - ele%value(ref_time_start$)
  endif

else
  if (orbit%beta == 0) then
    call out_io (s_fatal$, r_name, 'PARTICLE IN NON E-GUN ELEMENT HAS VELOCITY = 0!')
    if (global_com%exit_on_error) call err_exit
    time = orbit%t  ! Just to keep on going
    return
  endif
  time = -orbit%vec(5) / (orbit%beta * c_light)
endif

end function particle_ref_time

