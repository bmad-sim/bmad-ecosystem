!+
! Function particle_rf_time (orbit, ele, reference_active_edge, s_rel, time_coords) result (time)
!
! Routine to return the reference time used to calculate the phase of
! time-dependent EM fields.
!
! Essentially: The oscillations of EM fields are synched relative to the absolute clock if 
! absolute time tracking is used or are synched relative to the reference particle
! if relative time tracking is used. See the Bmad manual for more details.
!
! Also see set_particle_from_rf_time which is the inverse of this routine.
!
! Input:
!   orbit       -- Coord_struct: Particle coordinates
!   ele         -- ele_struct: Element being tracked through.
!   reference_active_edge 
!               -- logical: If True, and ele is a rfcavity or lcavity, use the active edge as the reference point.
!   s_rel       -- real(rp), optional: Longitudinal position relative to the upstream edge of the element.
!                   Needed for relative time tracking when the particle is inside the element. Default is 0.
!   time_coords -- logical, optional: Default False. If True then orbit is using time based coordinates.
!
! Ouput:
!   time      -- Real(rp): Current time.
!-

function particle_rf_time (orbit, ele, reference_active_edge, s_rel, time_coords) result (time)

use equal_mod, dummy_except => particle_rf_time
use attribute_mod, only: has_attribute

implicit none

type (coord_struct) orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: ref_ele
type (ele_pointer_struct), allocatable :: chain(:)

real(rp) time, s_hard_offset, beta0
real(rp), optional :: s_rel
integer ix_pass, n_links
logical reference_active_edge, abs_time
logical, optional :: time_coords

character(*), parameter :: r_name = 'particle_rf_time'

!

ref_ele => ele
if (ref_ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) ref_ele => pointer_to_lord (ref_ele, 1)

call multipass_chain(ref_ele, ix_pass, n_links, chain)
if (ix_pass > 1) ref_ele => chain(1)%ele

! With absolute time tracking the reference time is relative to the reference time of the element.
! This way the phase does not have to be adjusted when switching between absolute and relative time tracking.
! Note: With a multipass_slave, use the reference time of the pass element.
! Note: e_gun uses absolute time tracking to get around the problem when orbit%beta = 0.

if (absolute_time_tracking(ele)) then
  if (bmad_com%absolute_time_ref_shift) then
    time = orbit%t - ref_ele%value(ref_time_start$)
  else
    time = orbit%t
  endif

else
  if (orbit%beta == 0) then
    call out_io (s_fatal$, r_name, 'PARTICLE IN NON E-GUN ELEMENT HAS VELOCITY = 0!')
    if (global_com%exit_on_error) call err_exit
    time = orbit%t  ! Just to keep on going
    return
  endif

  if (logic_option(.false., time_coords)) then
    time = orbit%dt_ref
  else
    time = -orbit%vec(5) / (orbit%beta * c_light)
  endif

  if (present(s_rel)) then
    ! The effective reference velocity is different from the velocity of the reference particle for wigglers where the reference particle
    ! is not traveling along the reference line and in elements where the reference velocity is not constant.
    if (ele%value(delta_ref_time$) == 0 .or. ele%key == patch$) then
      beta0 = ref_ele%value(p0c$) / ref_ele%value(e_tot$) ! Singular case. 
    else
      beta0 = ref_ele%value(l$) / (c_light * ref_ele%value(delta_ref_time$))
    endif
    time = time + (s_rel + ele%s_start - ref_ele%s_start) / (c_light * beta0)
  endif
endif

!

if (reference_active_edge .and. (ele%key == rfcavity$ .or. ele%key == lcavity$)) then
  s_hard_offset = (ref_ele%value(l$) - ref_ele%value(l_active$)) / 2  
  beta0 = ele%value(p0c_start$) / ele%value(E_tot_start$)
  time = time - s_hard_offset / (c_light * beta0)
endif

end function particle_rf_time

