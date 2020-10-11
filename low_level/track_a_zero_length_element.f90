!+
! Subroutine track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)
!
! Subroutine to do tracking through a zero length elemnt.
! This routine is called by track1_runge_kutta and track1_time_runge_kutta since 
! these routines cannot handle zero length elements with finite magnetic multipoles.
!
! Input:
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element tracked through.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   end_orb    -- Coord_struct: Ending coords.
!   err_flag   -- Logical: Set True if there is an error. False otherwise.
!   track      -- Track_struct, optional: Structure holding the track information.
!- 

subroutine track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)

use fringe_mod, except_dummy => track_a_zero_length_element

implicit none

type (coord_struct) :: start_orb, end_orb
type (lat_param_struct), target, intent(inout) :: param
type (ele_struct), target, intent(inout) :: ele
type (track_struct), optional :: track
type (fringe_edge_info_struct) fringe_info

real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), amp
integer ix_pole_max
logical err_flag

! Since Bmad always uses electric field parameters and never integrated electric field parameters,
! the kick due to any electric field over zero length must be zero.

if (present(track)) call save_a_step (track, ele, param, .false., start_orb, 0.0_rp)

end_orb = start_orb
call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false.)

nullify(fringe_info%hard_ele)
fringe_info%particle_at = first_track_edge$
call apply_element_edge_kick(end_orb, fringe_info, ele, param, .false.)

!

call multipole_ele_to_ab(ele, .false., ix_pole_max, an, bn, magnetic$, include_kicks$)
if (ix_pole_max > -1) then
  amp = ac_kicker_amp(ele, end_orb)
  call ab_multipole_kicks (an, bn, ix_pole_max, param%particle, ele, end_orb, magnetic$, amp)
endif

!

fringe_info%particle_at = second_track_edge$
call apply_element_edge_kick(end_orb, fringe_info, ele, param, .false.)

call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false.)

if (present(track)) call save_a_step (track, ele, param, .false., start_orb, 0.0_rp)

end subroutine
