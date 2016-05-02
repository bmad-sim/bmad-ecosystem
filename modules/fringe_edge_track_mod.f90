module fringe_edge_track_mod

use bmad_utils_mod

implicit none

private does_this_ele_contain_the_next_edge 

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine calc_next_fringe_edge (track_ele, track_direction, s_edge_track, fringe_info, init_needed, orbit)
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
!   track_direction -- integer: +1 -> +s direction, -1 -> -s direction.
!   init_needed     -- logical, optional: If present and True then initialize.
!   orbit           -- coord_struct, optional: Used for initialization. Must be present if init_needed = True.
!
! Output:
!   s_edge_track -- Real(rp): S position of next hard edge in track_ele frame.
!                     If there are no more hard edges then s_edge_track will be set to ele%value(l$) if 
!                     track_direction = 1, and set to 0 if track_direction = -1.
!   fringe_info    -- fringe_edge_info_struct: Information on the next fringe to track through.
!-

subroutine calc_next_fringe_edge (track_ele, track_direction, s_edge_track, fringe_info, init_needed, orbit)

type (ele_struct), target :: track_ele
type (fringe_edge_info_struct) fringe_info
type (ele_struct), pointer :: lord
type (coord_struct), optional :: orbit

real(rp) s_edge_track
integer i, dir, track_direction
logical, optional :: init_needed

character(*), parameter :: r_name = 'calc_next_finge_edge'

! Situations with different fringe plancements to consider:
!   - Slice_slaves.
!   - jumbo super_slaves.
!   - super_slaves with multiple lords.
!   - an x_offset in an sbend will shift the s-position of the ends of the bend.
!   - An element with edge drifts (RF cavity, eg) will have fringes not at the element ends.

! Init if needed.
! Keep track of where particle is with respect to element edge using ele%ixx = upsteam_end$, downstream_end$, inside$
! upsteam_end$ means outside element on the upsteam side, downsteam_end$ manes outside element on the downstream side.
! The routine apply_element_edge_kick will modify ele%ixx as appropriate when the particle is tracked through an edge.

dir = track_direction

if (logic_option(.false., init_needed)) then
  nullify(fringe_info%hard_ele)
  if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
    call re_allocate(fringe_info%location, track_ele%n_lord)
    do i = 1, track_ele%n_lord
      lord => pointer_to_lord(track_ele, i)
      call init_this_ele (lord, i)
    enddo
  else
    call re_allocate(fringe_info%location, 1)
    call init_this_ele (track_ele, 1)
  endif
endif

! Find next hard edge. 

if (dir == 1) then
  s_edge_track = track_ele%value(l$)
else
  s_edge_track = 0
endif

nullify (fringe_info%hard_ele)
fringe_info%particle_at = none$

if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
  do i = 1, track_ele%n_lord
    lord => pointer_to_lord(track_ele, i)
    if (lord%key == overlay$ .or. lord%key == group$) cycle
    call does_this_ele_contain_the_next_edge (lord, i, track_ele, dir, s_edge_track, fringe_info)
  enddo
else
  call does_this_ele_contain_the_next_edge (track_ele, 1, track_ele, dir, s_edge_track, fringe_info)
endif

!-------------------------------------------------------------------------
contains 

subroutine init_this_ele (this_ele, ix_loc)

type (ele_struct) this_ele
real(rp) s_off, s1, s2, s_hard_upstream, s_hard_downstream, s_orb, ds_small
integer ix_loc

!

s_off = (this_ele%s - this_ele%value(l$)) - (track_ele%s - track_ele%value(l$))
s1 = s_off + (this_ele%value(l$) - hard_edge_model_length(this_ele)) / 2 
s2 = s_off + (this_ele%value(l$) + hard_edge_model_length(this_ele)) / 2 

s_hard_upstream   = min(s1, s2)
s_hard_downstream = max(s1, s2)

s_orb = orbit%s - (track_ele%s - track_ele%value(l$)) - track_ele%value(z_offset_tot$)
ds_small = bmad_com%significant_length

if (dir == 1) then
  if (orbit%location == upstream_end$) then
    if (s_hard_upstream > s_orb - ds_small) then
      fringe_info%location(ix_loc) = upstream_end$
    elseif (s_hard_downstream > s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = downstream_end$
    endif

  elseif (orbit%location == inside$) then
    if (s_hard_upstream > s_orb + ds_small) then
      fringe_info%location(ix_loc) = upstream_end$
    elseif (s_hard_downstream > s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = downstream_end$
    endif

  else
    call out_io (s_fatal$, r_name, 'CONFUSED FORWARD DIRECTION INITIALIZATION!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

else  ! dir = -1
  if (orbit%location == downstream_end$) then
    if (s_hard_downstream < s_orb + ds_small) then
      fringe_info%location(ix_loc) = downstream_end$
    elseif (s_hard_upstream < s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = upstream_end$
    endif

  elseif (orbit%location == inside$) then
    if (s_hard_downstream < s_orb - ds_small) then
      fringe_info%location(ix_loc) = downstream_end$
    elseif (s_hard_downstream < s_orb) then
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

end subroutine calc_next_fringe_edge

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine does_this_ele_contain_the_next_edge (this_ele, ix_loc, track_ele, dir, s_edge_track, fringe_info)

type (ele_struct), target :: this_ele, track_ele
type (fringe_edge_info_struct), target :: fringe_info

real(rp) s_this_edge, s1, s2, s_hard_upstream, s_hard_downstream, s_off, s_edge_track, ds_small
integer this_end, dir, ix_loc

! Remamber: element length can be less than zero.

s_off = (this_ele%s - this_ele%value(l$)) - (track_ele%s - track_ele%value(l$))
s1 = s_off + (this_ele%value(l$) - hard_edge_model_length(this_ele)) / 2 
s2 = s_off + (this_ele%value(l$) + hard_edge_model_length(this_ele)) / 2 

s_hard_upstream   = min(s1, s2)
s_hard_downstream = max(s1, s2)

ds_small = bmad_com%significant_length / 100

if (dir == 1) then
  select case (fringe_info%location(ix_loc))
  case (upstream_end$)
    s_this_edge = s_hard_upstream
    this_end = first_track_edge$
  case (inside$)
    s_this_edge = s_hard_downstream
    this_end = second_track_edge$
  case (downstream_end$)
    return
  case default
    call err_exit
  end select
  if (s_this_edge > s_edge_track + ds_small) return

else
  select case (fringe_info%location(ix_loc))
  case (upstream_end$)
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
  if (s_this_edge < s_edge_track - ds_small) return
endif

! This looks like the next hard edge

s_edge_track = s_this_edge

fringe_info%hard_ele => this_ele
fringe_info%particle_at = this_end
fringe_info%s_edge_hard = s_edge_track - s_off
fringe_info%hard_location => fringe_info%location(ix_loc)

end subroutine does_this_ele_contain_the_next_edge

end module
