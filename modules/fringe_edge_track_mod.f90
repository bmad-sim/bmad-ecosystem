module fringe_edge_track_mod

use bmad_utils_mod

implicit none

private does_this_ele_contain_the_next_edge 

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine calc_next_fringe_edge (track_ele, track_direction, s_edge_track, fringe_info, orbit, init_needed)
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
!   orbit           -- coord_struct: Particle position
!
! Output:
!   s_edge_track -- Real(rp): S position of next hard edge in track_ele frame.
!                     If there are no more hard edges then s_edge_track will be set to ele%value(l$) if 
!                     track_direction = 1, and set to 0 if track_direction = -1.
!   fringe_info    -- fringe_edge_info_struct: Information on the next fringe to track through.
!-

subroutine calc_next_fringe_edge (track_ele, track_direction, s_edge_track, fringe_info, orbit, init_needed)

type (ele_struct), target :: track_ele
type (fringe_edge_info_struct) fringe_info
type (ele_struct), pointer :: lord
type (coord_struct) :: orbit

real(rp) s_edge_track, s_orb
integer i, dir, track_direction, num_lords
logical, optional :: init_needed
logical has_z_offset

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
! Keep track of where particle is with respect to element edge using ele%ixx = upsteam_end$, downstream_end$, inside$
! upsteam_end$ means outside element on the upsteam side, downsteam_end$ manes outside element on the downstream side.
! The routine apply_element_edge_kick will modify ele%ixx as appropriate when the particle is tracked through an edge.

dir = track_direction
s_orb = orbit%s - (track_ele%s_start + track_ele%value(z_offset_tot$))

if (logic_option(.false., init_needed)) then
  nullify(fringe_info%hard_ele)

  if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
    call re_allocate(fringe_info%location, track_ele%n_lord)
    num_lords = 0; has_z_offset = .false.
    do i = 1, track_ele%n_lord
      lord => pointer_to_lord(track_ele, i)
      if (lord%key == overlay$ .or. lord%key == group$) cycle
      num_lords = num_lords + 1
      ! Can handle jumbo super_lord with z_offset but not regular super_lord.
      if (lord%value(z_offset_tot$) /= 0 .and. lord%value(lord_pad1$) == 0 .and. &
                                               lord%value(lord_pad2$) == 0) has_z_offset = .true.
      call init_this_ele (lord, i)
    enddo

    if (num_lords > 1 .and. has_z_offset) then
      call out_io (s_error$, r_name, 'TRACKING INTEGRATION THROUGH A SLAVE ELEMENT WITH MULTIPLE LORDS AND', &
                                     'ANY LORD HAVING A FINITE Z_OFFSET WILL NOT BE ACCURATE!', &
                                     'FOR ELEMENT: ' // track_ele%name)
    endif

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
    call does_this_ele_contain_the_next_edge (lord, i, track_ele, dir, s_edge_track, s_orb, fringe_info)
  enddo
else
  call does_this_ele_contain_the_next_edge (track_ele, 1, track_ele, dir, s_edge_track, s_orb, fringe_info)
endif

!-------------------------------------------------------------------------
contains 

subroutine init_this_ele (this_ele, ix_loc)

type (ele_struct) this_ele
real(rp) s_off, s1, s2, s_hard_upstream, s_hard_downstream, ds_small
integer ix_loc

!

s_off = this_ele%s_start - track_ele%s_start
s1 = s_off + (this_ele%value(l$) - hard_edge_model_length(this_ele)) / 2 
s2 = s_off + (this_ele%value(l$) + hard_edge_model_length(this_ele)) / 2 

s_hard_upstream   = min(s1, s2)
s_hard_downstream = max(s1, s2)

ds_small = bmad_com%significant_length

if (dir == 1) then
  if (orbit%location == upstream_end$) then
    fringe_info%location(ix_loc) = upstream_end$

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
      fringe_info%location(ix_loc) = downstream_end$

  elseif (orbit%location == inside$) then
    if (s_hard_downstream < s_orb - ds_small) then
      fringe_info%location(ix_loc) = downstream_end$
    elseif (s_hard_upstream < s_orb) then
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

subroutine does_this_ele_contain_the_next_edge (this_ele, ix_loc, track_ele, dir, s_edge_track, s_orb, fringe_info)

type (ele_struct), target :: this_ele, track_ele
type (fringe_edge_info_struct), target :: fringe_info

real(rp) s_this_edge, s1, s2, s_hard_upstream, s_hard_downstream, s_off, s_edge_track, ds_small, s_orb
integer this_end, dir, ix_loc

! Remamber: element length can be less than zero.

s_off = this_ele%s_start - track_ele%s_start
s1 = s_off + (this_ele%value(l$) - hard_edge_model_length(this_ele)) / 2 
s2 = s_off + (this_ele%value(l$) + hard_edge_model_length(this_ele)) / 2 

s_hard_upstream   = min(s1, s2)
s_hard_downstream = max(s1, s2)

ds_small = bmad_com%significant_length / 100

if (dir == 1) then
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
    return
  case default
    call err_exit
  end select

  if (s_this_edge > s_edge_track + ds_small) return

!

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
fringe_info%ds_edge = s_edge_track + s_off - s_orb

end subroutine does_this_ele_contain_the_next_edge

end module
