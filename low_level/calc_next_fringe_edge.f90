!+
! Subroutine calc_next_fringe_edge (track_ele, s_edge_body, fringe_info, orbit, init_needed, time_tracking)
!
! Routine to locate the next "hard edge" in an element when a hard edge model is being used. 
! This routine is used by integration tracking routines like Runge-Kutta.
! This routine is called repeatedly as the integration routine tracks through the element.
! If the element is a super_slave, there are potentially many hard edges.
!
! Rule: When track_ele is a super_slave, the edges of track_ele may be inside the field of a super_lord.
! In this case, the hard edge is applied at the edges to make the particle look like it is
! outside the field. This is done to make the tracking symplectic from entrance edge to exit
! edge. [Remember that Bmad's coordinates are not canonical inside a field region.]
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
!                     orbit%direction*ele%orientation = 1, and set to 0 otherwise.
!   fringe_info    -- fringe_edge_info_struct: Information on the next fringe to track through.
!-

subroutine calc_next_fringe_edge (track_ele, s_edge_body, fringe_info, orbit, init_needed, time_tracking)

use equal_mod, dummy => calc_next_fringe_edge

implicit none

type (ele_struct), target :: track_ele
type (fringe_edge_info_struct) fringe_info
type (ele_struct), pointer :: lord
type (coord_struct) :: orbit

real(rp) s_edge_body, s_orb
integer i, dir
logical, optional :: init_needed, time_tracking

character(*), parameter :: r_name = 'calc_next_finge_edge'

! Situations with different fringe plancements to consider:
!   - Slice_slaves.
!   - jumbo super_slaves.
!   - super_slaves with multiple lords.
!   - an x_offset in an sbend will shift the s-position of the ends of the bend.
!   - An element with edge drifts (RF cavity, eg) will have fringes not at the element ends.
!   - With time Runge-Kutta a particle can turn around.
!   - An element may have negative length.

! Init if needed.
! Keep track of where particle is with respect to element edge using %location = upsteam_end$, downstream_end$, inside$.
! upsteam_end$ means particle is outside of the element on the upsteam side, and 
! downsteam_end$ means outside the element on the downstream side.
! The routine apply_element_edge_kick will modify %location as appropriate when the particle is tracked through an edge.

s_orb = orbit%s - track_ele%s_start
dir = 1
if (track_ele%value(l$) < 0 .and. .not. logic_option(.false., time_tracking)) dir = -dir

if (logic_option(.false., init_needed)) then
  nullify(fringe_info%hard_ele)

  if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
    call re_allocate(fringe_info%location, track_ele%n_lord)
    do i = 1, track_ele%n_lord
      lord => pointer_to_lord(track_ele, i)
      if (lord%key == overlay$ .or. lord%key == group$) cycle
      call init_this_ele (lord, i, dir)
    enddo

  else
    call re_allocate(fringe_info%location, 1)
    call init_this_ele (track_ele, 1, dir)
  endif
endif

! patch element does not have a fringe and the edge is, by convention, at s = 0 independent of the 
! direction of travel.

if (track_ele%key == patch$) then
  s_edge_body = 0
  return
endif

! Find next hard edge. 
! Initially set s_edge_body to be the "last" edge consistant with with particle's direction and ele orientaiton.

if (orbit%direction == 1) then
  s_edge_body = track_ele%value(l$)
else
  s_edge_body = 0
endif

! Now test all other edges to see if they are nearer the particle.

nullify (fringe_info%hard_ele)
fringe_info%particle_at = none$

if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
  do i = 1, track_ele%n_lord
    lord => pointer_to_lord(track_ele, i)
    if (lord%key == overlay$ .or. lord%key == group$) cycle
    call does_this_ele_contain_the_next_edge (lord, i, track_ele, dir, orbit, s_edge_body, s_orb, fringe_info)
  enddo

else
  call does_this_ele_contain_the_next_edge (track_ele, 1, track_ele, dir, orbit, s_edge_body, s_orb, fringe_info)
endif

!-------------------------------------------------------------------------
contains 

subroutine init_this_ele (this_ele, ix_loc, dir)

type (ele_struct) this_ele
real(rp) s_off, s1, s2, s_hard_upstream, s_hard_downstream, ds_small
integer ix_loc, dir

!


s_off = this_ele%s_start - track_ele%s_start
s1 = s_off + (this_ele%value(l$) - hard_edge_model_length(this_ele)) / 2 
s2 = s_off + (this_ele%value(l$) + hard_edge_model_length(this_ele)) / 2 

if (dir > 0) then
  s_hard_upstream   = min(s1, s2)
  s_hard_downstream = max(s1, s2)
else
  s_hard_upstream   = max(s1, s2)
  s_hard_downstream = min(s1, s2)
endif

ds_small = bmad_com%significant_length

if (orbit%direction == 1) then
  if (orbit%location == upstream_end$) then
    fringe_info%location(ix_loc) = upstream_end$

  elseif (orbit%location == inside$) then
    if (dir * s_hard_upstream > dir * (s_orb + ds_small)) then
      fringe_info%location(ix_loc) = upstream_end$
    elseif (dir * s_hard_downstream > dir * s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = downstream_end$
    endif

  else
    call out_io (s_fatal$, r_name, 'CONFUSED FORWARD DIRECTION INITIALIZATION!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

else  ! orbit%direction = -1
  if (orbit%location == downstream_end$) then
      fringe_info%location(ix_loc) = downstream_end$

  elseif (orbit%location == inside$) then
    if (dir * s_hard_downstream < dir * (s_orb - ds_small)) then
      fringe_info%location(ix_loc) = downstream_end$
    elseif (dir * s_hard_upstream < dir * s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = upstream_end$
    endif

  else
    call out_io (s_fatal$, r_name, 'CONFUSED REVERSE DIRECTION INITIALIZATION!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

end subroutine init_this_ele

!---------------------------------------------------------------------------
! contains

subroutine does_this_ele_contain_the_next_edge (this_ele, ix_loc, track_ele, dir, orbit, s_edge_body, s_orb, fringe_info)

type (ele_struct), target :: this_ele, track_ele
type (fringe_edge_info_struct), target :: fringe_info
type (coord_struct) orbit

real(rp) s_this_edge, s1, s2, s_hard_upstream, s_hard_downstream, s_off, s_edge_body, ds_small, s_orb, leng
integer this_end, ix_loc, dir, rel_dir

! Remamber: element length can be less than zero.

leng = this_ele%value(l$)
rel_dir = dir * orbit%direction
s_off = this_ele%s_start - track_ele%s_start

s1 = s_off + (leng - hard_edge_model_length(this_ele)) / 2 
s2 = s_off + (leng + hard_edge_model_length(this_ele)) / 2 

! With a solenoid must always apply the fringe kick due to the longitudinal field. 
! If not done the matrix calc will not be symplectic.
! For other elements, especially quadrupoles, this is problematic due to the soft edge kick not being being exactly the reverse
! going from inside to outside and vice versa (it is confusing if a superimposed marker shifts the tracking).

if (track_ele%key == solenoid$ .or. track_ele%key == sol_quad$) then
  s1 = max(s1, 0.0_rp)
  s2 = min(s2, track_ele%value(l$))
endif

!

if (dir > 0) then
  s_hard_upstream   = min(s1, s2)
  s_hard_downstream = max(s1, s2)
else
  s_hard_upstream   = max(s1, s2)
  s_hard_downstream = min(s1, s2)
endif

ds_small = bmad_com%significant_length / 100

if (orbit%direction == 1) then
  select case (fringe_info%location(ix_loc))
  case (upstream_end$)
    ! e_gun does not have an entrance edge
    if (this_ele%key == e_gun$) then
      s_this_edge = s_hard_downstream
      this_end = second_track_edge$
    else
      s_this_edge = s_hard_upstream
      this_end = first_track_edge$
    endif
  case (inside$)
    s_this_edge = s_hard_downstream
    this_end = second_track_edge$
  case (downstream_end$)
    if (track_ele%orientation == -1) s_edge_body = leng - s_edge_body
    return
  case default
    call err_exit
  end select

  if (dir * s_this_edge > dir * s_edge_body + ds_small) return

!

else
  select case (fringe_info%location(ix_loc))
  case (upstream_end$)
    if (track_ele%orientation == -1) s_edge_body = leng - s_edge_body
    return
  case (inside$)
    s_this_edge = s_hard_upstream
    this_end = second_track_edge$
  case (downstream_end$)
    s_this_edge = s_hard_downstream
    this_end = first_track_edge$
  case default
    call err_exit
  end select

  if (dir * s_this_edge < dir * s_edge_body - ds_small) return
endif

! This looks like the next hard edge

fringe_info%hard_ele => this_ele
fringe_info%particle_at = this_end
fringe_info%hard_location => fringe_info%location(ix_loc)

s_edge_body = s_this_edge
fringe_info%s_edge_hard = s_edge_body - s_off
fringe_info%ds_edge = s_edge_body + s_off - s_orb

if (track_ele%orientation == -1) then
  s_edge_body = leng - s_edge_body
  fringe_info%s_edge_hard = leng - fringe_info%s_edge_hard
  fringe_info%ds_edge = -fringe_info%ds_edge
endif

end subroutine does_this_ele_contain_the_next_edge

end subroutine calc_next_fringe_edge

