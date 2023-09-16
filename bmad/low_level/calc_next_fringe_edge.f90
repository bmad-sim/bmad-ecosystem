!+
! Subroutine calc_next_fringe_edge (track_ele, s_edge_body, fringe_info, orbit, init_needed, time_tracking)
!
! Routine to locate the next "hard edge" in an element when a hard edge model is being used. 
! This routine is used by integration tracking routines like Runge-Kutta.
! This routine is called repeatedly as the integration routine tracks through the element.
! If the element is a super_slave, there are potentially many hard edges.
!
! Note: An element that has a longitudinal field component or an electric field is considered to 
! have a fringe at both ends even though the slice ends may not actually abut the ends of the full element.
! This is done to make computed transfer matrices symplectic. Non-symplectic matrices
! will mess up Twiss calculations. EG: Think of the case of a solenoid with a superimposed marker.
!
! Note: For other elements, especially quadrupoles, always applying edge effects is potentially problematic 
! due to the soft edge kick not being being exactly the reverse going from inside to outside and vice versa.
! So applying an edge kick could be confusing since a superimposed marker would shift the tracking.
!
! Input:
!   track_ele       -- ele_struct: Element being tracked through.
!   orbit           -- coord_struct: Particle position
!   init_needed     -- logical, optional: If present and True then initialize.
!   time_tracking   -- logical, optional: If present and True then this routine is being called by the
!                        time Runge-Kutta tracker. Default is False.
!
! Output:
!   s_edge_body -- Real(rp): S position of next hard edge in track_ele body frame.
!                     If there are no more hard edges then s_edge_body will be set to ele%value(l$) if 
!                     orbit%direction*orbit%time_dir*ele%orientation = 1, and set to 0 otherwise.
!   fringe_info    -- fringe_field_info_struct: Information on the next fringe to track through.
!-

subroutine calc_next_fringe_edge (track_ele, s_edge_body, fringe_info, orbit, init_needed, time_tracking)

use equal_mod, dummy => calc_next_fringe_edge

implicit none

type (ele_struct), target :: track_ele
type (fringe_field_info_struct) fringe_info
type (ele_struct), pointer :: lord
type (coord_struct) :: orbit

real(rp) s_edge_body
integer i, leng_sign
logical, optional :: init_needed, time_tracking

character(*), parameter :: r_name = 'calc_next_finge_edge'

! Situations with different fringe plancements to consider:
!   - Slice_slaves.
!   - Jumbo super_slaves.
!   - Super_slaves with multiple lords can have multiple fringes at the same s-position.
!   - An x_offset in an sbend will shift the s-position of the ends of the bend.
!   - An element with edge drifts (RF cavity, eg) will have fringes not at the element ends.
!   - With time Runge-Kutta a particle can turn around.
!   - An element may have negative length.

! Init if needed.
! Keep track of where particle is with respect to element edge using %location = entrance_end$, exit_end$, inside$.
! upsteam_end$ means particle is outside of the element on the upsteam side, and 
! downsteam_end$ means outside the element on the downstream side.
! The routine apply_element_edge_kick will modify %location as appropriate when the particle is tracked through an edge.

leng_sign = 1
if (track_ele%value(l$) < 0 .and. .not. logic_option(.false., time_tracking)) leng_sign = -leng_sign

if (logic_option(.false., init_needed)) call init_fringe_info (fringe_info, track_ele, orbit, leng_sign)

! patch element does not have a fringe and the edge is, by convention, at s = 0 independent of the 
! direction of travel.

if (track_ele%key == patch$) then
  s_edge_body = 0
  return
endif

! Initially set s_edge_body to be the "last" edge consistant with with particle's direction and ele orientaiton.

if (orbit%direction*orbit%time_dir * track_ele%orientation == 1) then
  s_edge_body = track_ele%value(l$)
else
  s_edge_body = 0
endif

if (.not. fringe_info%has_fringe) return

!-----------------------------------------

! Now test all other edges to see if they are nearer the particle.

nullify (fringe_info%hard_ele)
fringe_info%particle_at = none$

if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
  do i = 1, track_ele%n_lord
    lord => pointer_to_lord(track_ele, i)
    if (lord%key == overlay$ .or. lord%key == group$) cycle
    call does_this_ele_contain_the_next_edge (lord, i, track_ele, leng_sign, orbit, s_edge_body, fringe_info)
  enddo

else
  call does_this_ele_contain_the_next_edge (track_ele, 1, track_ele, leng_sign, orbit, s_edge_body, fringe_info)
endif

!---------------------------------------------------------------------------
contains

subroutine does_this_ele_contain_the_next_edge (this_ele, ix_this_ele, track_ele, leng_sign, orbit, s_edge_body, fringe_info)

type (ele_struct), target :: this_ele, track_ele
type (fringe_field_info_struct), target :: fringe_info
type (coord_struct) orbit

real(rp) s_this_edge, s1, s2, s_hard_entrance, s_hard_exit, s_off, s_edge_body, ds_small, s_orb, this_leng, track_leng
integer this_end, ix_this_ele, leng_sign

!

if (fringe_info%location(ix_this_ele) == nowhere$) return

! This_ele is either a lord or the track_ele itself.
! Remamber: element length can be less than zero.

if (track_ele%orientation == 1) then
  s_off = this_ele%s_start - track_ele%s_start
  s_orb = orbit%s - track_ele%s_start
else
  s_off = track_ele%s - this_ele%s
  s_orb = track_ele%s - orbit%s
endif

this_leng = this_ele%value(l$)
track_leng = track_ele%value(l$)

select case (this_ele%key)
case (rfcavity$, lcavity$)
  s1 = s_off + (this_leng - this_ele%value(l_active$)) / 2  ! Distance from entrance end to active edge
  s2 = s_off + (this_leng + this_ele%value(l_active$)) / 2  ! Distance from entrance end to the other active edge
case default
  s1 = s_off         ! Distance from track_ele entrance end to entrance hard edge
  s2 = s_off + this_leng  ! Distance from track_ele entrance end to the exit hard edge
end select

! Make sure that fringes are not outside of track_ele. This is done so calculated maps are symplectic.

s1 = max(s1, 0.0_rp)
s2 = min(s2, track_leng)

!

if (leng_sign > 0) then
  s_hard_entrance = min(s1, s2)
  s_hard_exit     = max(s1, s2)
else
  s_hard_entrance = max(s1, s2)
  s_hard_exit     = min(s1, s2)
endif

ds_small = bmad_com%significant_length / 100

if (orbit%direction*orbit%time_dir * track_ele%orientation == 1) then
  select case (fringe_info%location(ix_this_ele))
  case (entrance_end$)
    ! e_gun does not have an entrance edge
    if (this_ele%key == e_gun$) then
      s_this_edge = s_hard_exit
      this_end = second_track_edge$
    else
      s_this_edge = s_hard_entrance
      this_end = first_track_edge$
      if (s_this_edge < min(0.0_rp, track_leng) .or. s_this_edge > max(0.0_rp, track_leng)) then
        s_this_edge = s_hard_exit
        this_end = second_track_edge$
      endif
    endif
  case (inside$)
    s_this_edge = s_hard_exit
    this_end = second_track_edge$
  case (exit_end$)
    return
  case default
    call err_exit
  end select

  if (leng_sign * s_this_edge > leng_sign * s_edge_body + ds_small) return

!

else
  select case (fringe_info%location(ix_this_ele))
  case (entrance_end$)
    return
  case (inside$)
    s_this_edge = s_hard_entrance
    this_end = second_track_edge$
  case (exit_end$)
    s_this_edge = s_hard_exit
    this_end = first_track_edge$
    if (s_this_edge < min(0.0_rp, track_leng) .or. s_this_edge > max(0.0_rp, track_leng)) then
      s_this_edge = s_hard_entrance
      this_end = second_track_edge$
    endif
  case default
    call err_exit
  end select

  if (leng_sign * s_this_edge < leng_sign * s_edge_body - ds_small) return
endif

! 

if (s_this_edge < min(0.0_rp, track_leng) .or. s_this_edge > max(0.0_rp, track_leng)) return

fringe_info%hard_ele => this_ele
fringe_info%particle_at = this_end
fringe_info%hard_location => fringe_info%location(ix_this_ele)

s_edge_body = s_this_edge
fringe_info%s_edge_hard = s_edge_body - s_off
fringe_info%ds_edge = s_edge_body + s_off - s_orb

end subroutine does_this_ele_contain_the_next_edge

end subroutine calc_next_fringe_edge

