!+
! Module bmad_utils_mod
!
! Module for some core utility routines.
!-

module bmad_utils_mod

use basic_attribute_mod
use equality_mod

implicit none

private pointer_to_ele1, pointer_to_ele2
private pointer_to_branch_given_name, pointer_to_branch_given_ele

!+
! Function pointer_to_branch
!
! Routine to return a pointer to the lattice branch associated with a given name
! or a given element.
!
! This routine is an overloaded name for:
!   pointer_to_branch_given_ele (ele) result (branch_ptr, parameter_is_branch0)
!   pointer_to_branch_given_name (branch_name, lat) result (branch_ptr)
!
! The lattice branch *associated* with a given element is not necessarily the
! branch where the element is *located*. For example, all lords live in branch #0.
! But the branch associated with a super_lord element is the branch of its slaves.
!
! To get the branch where the element is located, simply use ele%ix_branch.
! 
! Note: Result is ambiguous if ele argument is associated with multiple branches 
! which can happen, for example, with overlay elements.
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   ele         -- Ele_struct: Element contained in the branch.
!   branch_name -- Character(*): May be a branch name or a branch index.
!   lat         -- Lat_struct: Lattice to search.
!   parameter_is_branch0 -- logical, optional: If True, 'PARAMETER' is taken to be
!                     an alternative name for branch(0). Default is False.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the branch.
!                   Nullified if there is no associated branch.
!-

interface pointer_to_branch
  module procedure pointer_to_branch_given_ele
  module procedure pointer_to_branch_given_name
end interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_ele (...)
!
! Routine to return a pointer to an element.
! pointer_to_ele is an overloaded name for:
!     Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
!     Function pointer_to_ele2 (lat, ele_loc_id) result (ele_ptr)
!
! Also see:
!   pointer_to_slave
!   pointer_to_lord
!
! Module needed:
!   use bmad_utils_mod
!
! Input:
!   lat       -- lat_struct: Lattice.
!   ix_ele    -- Integer: Index of element in lat%branch(ix_branch)
!   ix_branch -- Integer: Index of the lat%branch(:) containing the element.
!   ele_loc   -- Lat_ele_loc_struct: Location identification.
!
! Output:
!   ele_ptr  -- Ele_struct, pointer: Pointer to the element. 
!-

interface pointer_to_ele
  module procedure pointer_to_ele1
  module procedure pointer_to_ele2
end interface

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function relative_tracking_charge (orbit, param) result (rel_charge)
!
! Routine to determine the relative charge/mass of the particle being 
! tracked relative to the charge of the reference particle.
!
! Input:
!   orbit -- coord_struct: Particle position structure.
!   param -- lat_param_struct: Structure holding the reference particle id.
!
! Output:
!   rel_charge -- real(rp): Relative charge/mass
!-

function relative_tracking_charge (orbit, param) result (rel_charge)

type (coord_struct) orbit
type (lat_param_struct) param
real(rp) rel_charge

!

if (orbit%species == photon$ .or. param%particle == photon$) then
  rel_charge = 0
else
  rel_charge = (charge_of(orbit%species) / mass_of(orbit%species)) / &
                 (charge_of(param%particle) / mass_of(param%particle))
endif

end function relative_tracking_charge

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine set_orbit_to_zero (orbit, n1, n2, ix_noset)
!
! Routine to set the orbit of a subset orbit(n1:n2) of a coord_struct array to zero.
!
! Input:
!   n1        -- integer: Lower bound of orbit(:) array subset.
!   n2        -- integer: Upper bound of orbit(:) array subset.
!   ix_noset  -- integer, optional: If present then orbit(ix_noset) will not be zeroed.
!
! Output:
!   orbit(:)  -- coord_struct: Array with particle positions in the range orbit(n1:n2)
!                   set to zero except for orbit(ix_noset).
!-

subroutine set_orbit_to_zero (orbit, n1, n2, ix_noset)

type (coord_struct) orbit(0:)
integer n, n1, n2
integer, optional :: ix_noset

! Note: n will never equal -1 so the if statement does what it should do.

do n = n1, n2
  if (n == integer_option(-1, ix_noset)) cycle 
  orbit(n) = coord_struct()
enddo

end subroutine set_orbit_to_zero

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function lord_edge_aligned (slave, slave_edge, lord) result (is_aligned)
!
! Routine to determine if the edge of a super_lord is aligned with a given edge 
! of a super_slave or slice_slave.
!
! Input:
!   slave       -- ele_struct: Slave element.
!   slave_edge  -- integer: End under consideration: entrance_end$, exit_end$, in_between$, etc.
!   lord        -- ele_struct: Lord element.
!
! Output:
!   is_aligned  -- integer: True if a lord edge is aligned with the slave edge.
!                   If slave_edge is not entrance_end$ nor exit_end$ then is_aligned is False.
!- 

function lord_edge_aligned (slave, slave_edge, lord) result (is_aligned)

type (ele_struct), target :: slave, lord
type (branch_struct), pointer :: branch
integer slave_edge, ix_slave
real(rp) s_lord
logical is_aligned
character(*), parameter :: r_name = 'lord_edge_aligned'

! 

select case (stream_ele_end(slave_edge, slave%orientation))
case (upstream_end$)
  s_lord = lord%s_start
  branch => slave%branch
  if (associated(branch)) then
    if (s_lord < branch%ele(0)%s) s_lord = s_lord + branch%param%total_length
  endif
  is_aligned = (abs(slave%s_start - s_lord) < bmad_com%significant_length)

case (downstream_end$)
  is_aligned = (abs(slave%s - lord%s) < bmad_com%significant_length)

case default
  is_aligned = .false.
end select

end function lord_edge_aligned

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function at_this_ele_end (now_at, where_at) result (is_at_this_end)
!
! Routine to determine if an aperture or fringe field is present.
!
! Input:
!   now_at      -- Integer: Which end is under consideration: entrance_end$, exit_end$, surface$, or in_between$.
!   where_at    -- Integer: Which ends have the aperture or fringe field: entrance_end$, exit_end$, 
!                     continuous$, both_ends$, no_aperture$, surface$.
!
! Output:
!   is_at_this_end   -- Logical: True if at this end. False otherwise.
!- 

function at_this_ele_end (now_at, where_at) result (is_at_this_end)

integer now_at, where_at
logical is_at_this_end

!

if (where_at == no_aperture$) then
  is_at_this_end = .false.
  return
endif

if (now_at == surface$ .or. where_at == surface$) then
  is_at_this_end = (now_at == where_at)
  return
endif

if (where_at == continuous$) then
  is_at_this_end = .true.
  return
endif

!

select case (now_at)
case (entrance_end$)
  select case (where_at)
  case (entrance_end$, both_ends$); is_at_this_end = .true.
  case default;                     is_at_this_end = .false.
  end select

case (exit_end$)
  select case (where_at)
  case (exit_end$, both_ends$); is_at_this_end = .true.
  case default;                 is_at_this_end = .false.
  end select

case (in_between$)
  is_at_this_end = .false. ! where_at will not be continuous$
end select

end function at_this_ele_end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine set_fringe_on_off (fringe_at, ele_end, on_or_off) 
!
! Routine to modify a ele%value(fringe_at$) setting to either turn on or turn off a fringe
! at either the entrance or exit end of the element. 
!
! Input:
!   fringe_at       -- real(rp): Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
!   ele_end         -- integer: Element edge: entrance_end$ or exit_end$
!   on_or_off       -- integer: Turn on$ or off$
!
! Output:
!   fringe_at       -- real(rp): Modified fringe setting.
!- 

subroutine set_fringe_on_off (fringe_at, ele_end, on_or_off) 

integer ele_end, on_or_off
real(rp) fringe_at

!

select case (on_or_off)

! On

case (on$)

  select case (ele_end)
  case (entrance_end$)
    select case (nint(fringe_at))
    case (no_end$);   fringe_at = entrance_end$
    case (exit_end$); fringe_at = both_ends$
    end select

  case (exit_end$)
    select case (nint(fringe_at))
    case (no_end$);       fringe_at = exit_end$
    case (entrance_end$); fringe_at = both_ends$
    end select

  case default
    call err_exit  ! Should not be here
  end select

! Off

case (off$)

  select case (ele_end)
  case (entrance_end$)
    select case (nint(fringe_at))
    case (both_ends$);    fringe_at = exit_end$
    case (entrance_end$); fringe_at = no_end$
    end select

  case (exit_end$)
    select case (nint(fringe_at))
    case (both_ends$); fringe_at = entrance_end$
    case (exit_end$);  fringe_at = no_end$
    end select

  case default
    call err_exit  ! Should not be here
  end select

case default
  call err_exit  ! Should not be here
end select


end subroutine set_fringe_on_off

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function physical_ele_end (track_end, track_direction, ele_orientation, return_stream_end) result (physical_end)
!
! Rotine to determine which physical end of an element a particle is at given 
! the position in terms of upstream/downstream and the element's orientation
!
! Input:
!   track_end         -- Integer: first_track_edge$, second_track_edge$, surface$, or in_between$
!   track_direction   -- Integer: +1 or -1
!   ele_orientation   -- Integer: Either 1 = Normal or -1 = element reversed.
!   return_stream_end -- logical, optional: If True return the stream end instead of the physical end.
!                          Default is False.
!
! Output:
!   physical_end     -- Integer: Return_stream_end ->  Possibilities
!                                False             ->  entrance_end$, exit_end$, surface$, or in_between$
!                                True              ->  upstream_end$, downstream_end$ 
!-

function physical_ele_end (track_end, track_direction, ele_orientation, return_stream_end) result (physical_end)

integer track_end, track_direction, ele_orientation, physical_end
logical, optional :: return_stream_end
character(*), parameter :: r_name  = 'physical_ele_end'

!

if (track_end == surface$ .or. track_end == in_between$) then
  physical_end = track_end
  return
endif

!

if (logic_option(.false., return_stream_end)) then
  select case (ele_orientation * track_direction)
  case (1) 
    select case (track_end)
    case (first_track_edge$);   physical_end = upstream_end$
    case (second_track_edge$);  physical_end = downstream_end$
    end select

  case (-1)
    select case (track_end)
    case (first_track_edge$);   physical_end = downstream_end$
    case (second_track_edge$);  physical_end = upstream_end$
    end select
  end select

  return
endif

!

select case (ele_orientation * track_direction)

case (1) 
  select case (track_end)
  case (first_track_edge$);   physical_end = entrance_end$
  case (second_track_edge$); physical_end = exit_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD TRACK_END: \i0\ ', i_array = [track_end])
    if (global_com%exit_on_error) call err_exit
  end select

case (-1) 
  select case (track_end)
  case (first_track_edge$);   physical_end = exit_end$
  case (second_track_edge$); physical_end = entrance_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD TRACK_END: \i0\ ', i_array = [track_end])
    if (global_com%exit_on_error) call err_exit
  end select

case default
  call out_io (s_fatal$, r_name, 'BAD ELEMENT ORIENTATION: \2i4\ ', i_array = [ele_orientation, track_direction])
  if (global_com%exit_on_error) call err_exit

end select

end function physical_ele_end 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function stream_ele_end (physical_end, ele_orientation) result (stream_end)
!
! Rotine to determine which stream end of an element a particle is at given 
! the position in terms of the physical end and the element's orientation
!
! Input:
!   physical_end     -- Integer: entrance_end$, exit_end$, surface$, etc.
!   ele_orientation  -- Integer: Either 1 = Normal or -1 = element reversed.
!
! Output:
!   stream_end       -- Integer: upstream_end$, downstream_end$, or set equal
!                         to physical_end if physical_end is neither entrance_end$
!                         nor exit_end$
!-

function stream_ele_end (physical_end, ele_orientation) result (stream_end)

integer stream_end, ele_orientation, physical_end
character(*), parameter :: r_name  = 'stream_ele_end'

!

if (physical_end /= entrance_end$ .and. physical_end /= exit_end$) then
  stream_end = physical_end
  return
endif

!

select case (ele_orientation)

case (1) 
  select case (physical_end)
  case (entrance_end$);   stream_end = upstream_end$
  case (exit_end$);       stream_end = downstream_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD PHYSICAL_END: \i0\ ', i_array = [physical_end])
    if (global_com%exit_on_error) call err_exit
  end select

case (-1) 
  select case (physical_end)
  case (entrance_end$);   stream_end = downstream_end$
  case (exit_end$);       stream_end = upstream_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD PHYSICAL_END: \i0\ ', i_array = [physical_end])
    if (global_com%exit_on_error) call err_exit
  end select

case default
  call out_io (s_fatal$, r_name, 'BAD ELEMENT ORIENTATION: \i0\ ', i_array = [ele_orientation])
  if (global_com%exit_on_error) call err_exit

end select

end function stream_ele_end 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_if_s_in_bounds (branch, s, err_flag, translated_s)
!
! Routine to check if a given longitudinal position s is within the bounds of a given branch of a lattice.
! For linear branches the bounds are normally [0, branch_length].
! For circular branches negative s values do make sense so the bounds 
!   are normally [-branch_length, branch_length].
!
! "Normally" means that starting s-position in the branch is zero. This routine does
! adjust for non-zero starting s-positions.
!
! This routine will bomb the program if global_com%exit_on_error is True.
!
! Optionally: translated_s is a translated longitudinal position which is normally
! in the range [0, branch_length].
!
! Moduels needed:
!   use bmad
!
! Input:
!   branch        -- branch_struct: Branch
!   s             -- Real(rp): longitudinal position in the given branch.
!   
! Output:
!   err_flag      -- Logical: Set True if s position is out-of-bounds. False otherwise.
!   translated_s  -- Real(rp), optional: position translated to the range [0, branch_length]
!-

subroutine check_if_s_in_bounds (branch, s, err_flag, translated_s)

type (branch_struct) branch

real(rp) s, ss, s_min, s_max, ds_fudge, s_bound
real(rp), optional :: translated_s

logical err_flag

character(24), parameter :: r_name = 'check_if_s_in_bounds'

! Setup

s_min = branch%ele(0)%s
s_max = branch%ele(branch%n_ele_track)%s 
ds_fudge = bmad_com%significant_length
err_flag = .false.
ss = s

! Check

if (s > s_max + ds_fudge) then
  err_flag = .true.
  s_bound = s_max
elseif (branch%param%geometry == closed$) then
  if (s < s_min - (s_max - s_min) - ds_fudge) then
    err_flag = .true.
    s_bound = s_min - (s_max - s_min)
  endif
  if (s < s_min) ss = s + (s_max - s_min)
elseif (s < s_min - ds_fudge) then
  err_flag = .true.
  s_bound = s_min
endif

! Finish

if (err_flag) then
  call out_io (s_fatal$, r_name, &
        'S-POSITION \es20.12\ PAST EDGE OF LATTICE. ' , &
        'PAST LATTICE EDGE AT: \es20.12\ ', r_array = [s, s_bound])
  if (global_com%exit_on_error) call err_exit
endif

if (present(translated_s)) translated_s = ss

end subroutine check_if_s_in_bounds 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_loc_to_string (ele, show_branch0) result (str)
!
! Routine to encode an element's location into a string.
! Example output:
!   "34"     ! Input: lat%ele(34) which is equivalent to lat%branch(0)%ele(34)
!   "0>>34"  ! Same as above if show_branch0 is set to True.
!   "1>>56"  ! Input: lat%branch(1)%ele(56).
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Element in a lattice
!   show_branch0 -- Logical, optional: Explicitly show branch for main 
!                     lattice elements? Default is False.
!
! Output:
!   str(10)     -- Character: Output string. Left justified.
!-

function ele_loc_to_string (ele, show_branch0) result (str)

type (ele_struct) ele
logical, optional :: show_branch0

character(10) str

!

if (ele%ix_branch == 0 .and. .not. logic_option(.false., show_branch0)) then
  write (str, '(i0)') ele%ix_ele
else
  write (str, '(i0, a, i0)') ele%ix_branch, '>>', ele%ix_ele
endif

end function ele_loc_to_string 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function on_a_girder(ele) result (is_on_girder)
!
! Routine to determine if an element is being supported by a girder element.
!
! Input:
!   ele     -- ele_struct: Element to check.
!
! Output:
!   is_on_girder -- Logical: True if supported. False otherwise.
!- 

function on_a_girder(ele) result (is_on_girder)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
logical is_on_girder
integer i

!

is_on_girder = .false.

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  if (lord%key /= girder$) cycle
  is_on_girder = .true.
  return
enddo

end function on_a_girder

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_controller_controls (contrl, name, err)
!
! Routine to check for problems when setting up group or overlay controllers.
!
! Modules needed:
!   use bmad
!
! Input:
!   contrl(:)   -- Control_struct: control info. 1 element for each slave.
!   name        -- Character(*): Lord name. Used for error reporting.
!
! Output:
!   err         -- Logical: Set true if there is a problem. False otherwise.
!-

subroutine check_controller_controls (contrl, name, err)

type (control_struct), target :: contrl(:)
type (control_struct), pointer :: c1, c2
integer i, j
logical err
character(*) name
character(40) :: r_name = 'check_controller_controls'

!

err = .true.

do i = 1, size(contrl)
  c1 => contrl(i)
  do j = i+1, size(contrl)
    c2 => contrl(j)
    if (c1%slave == c2%slave .and. c1%ix_attrib == c2%ix_attrib) then
      call out_io (s_error$, r_name, 'DUPLICATE SLAVE CONTROL FOR LORD: ' // name)
      return
    endif
  enddo
enddo

err = .false.

end subroutine check_controller_controls

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function particle_is_moving_forward (orbit) result (is_moving_forward)
!
! Routine to determine if a particle is moving in the forward +s direction.
! If not moving forward it is dead or is moving backward.
!
! Remember: +s and +z directions are counteraligned if element being tracked 
! through is reversed.
!
! Input:
!   orbit     -- coord_struct: Particle coordinates
!
! Output:
!   is_moving_forward -- Logical: True if moving forward. False otherwise.
!-

function particle_is_moving_forward (orbit) result (is_moving_forward)

type (coord_struct) orbit
integer particle
logical is_moving_forward

!

is_moving_forward = (orbit%state == alive$) .and. (orbit%direction == 1)

end function particle_is_moving_forward

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function particle_is_moving_backwards (orbit) result (is_moving_backwards)
!
! Routine to determine if a particle is moving in the backward -s direction.
! If not moving backward it is dead or is moving backward.
!
! Remember: +s and +z directions are counteraligned if element being tracked 
! through is reversed.
!
! Input:
!   orbit  -- coord_struct: Particle coordinates
!
! Output:
!   is_moving_backwards -- Logical: True if moving backward. False otherwise.
!-

function particle_is_moving_backwards (orbit) result (is_moving_backwards)

type (coord_struct) orbit
logical is_moving_backwards

!

is_moving_backwards = (orbit%state == alive$) .and. (orbit%direction == -1)

end function particle_is_moving_backwards

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_has_nonzero_kick (ele) result (has_kick)
!
! Function to determine if an element has nonzero kick values.
! Kicks are something hkick$, bl_vkick$, etc.
! See also: zero_ele_kicks, ele_has_offset, zero_ele_offsets.
!
! Modules needed:
!   use bmad
!
! Input
!   ele -- Ele_struct: Element with possible nonzero kicks.
!
! Output:
!   ele -- Ele_struct: Element with no kicks.
!-

function ele_has_nonzero_kick (ele) result (has_kick)

type (ele_struct) ele
logical has_kick

!

has_kick = .false.

if (has_hkick_attributes(ele%key)) then
  if (ele%value(bl_hkick$) /= 0) has_kick = .true.
  if (ele%value(bl_vkick$) /= 0) has_kick = .true.

elseif (has_kick_attributes(ele%key)) then
  if (ele%value(bl_kick$) /= 0) has_kick = .true.
endif

if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
  if (ele%value(coupler_strength$) /= 0) has_kick = .true.
endif

end function ele_has_nonzero_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine zero_ele_kicks (ele)
!
! Subroutine to zero any kick attributes like hkick$, bl_vkick$, etc.
! See also: ele_has_nonzero_kick, ele_has_offset, zero_ele_offsets.
!
! Modules needed:
!   use bmad
!
! Input
!   ele -- Ele_struct: Element with possible nonzero kicks.
!
! Output:
!   ele -- Ele_struct: Element with no kicks.
!-

subroutine zero_ele_kicks (ele)

type (ele_struct) ele

!

if (has_hkick_attributes(ele%key)) then
  ele%value(hkick$) = 0
  ele%value(vkick$) = 0
  ele%value(bl_hkick$) = 0
  ele%value(bl_vkick$) = 0

elseif (has_kick_attributes(ele%key)) then
  ele%value(kick$) = 0
  ele%value(bl_kick$) = 0
endif

if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
  ele%value(coupler_strength$) = 0
endif

end subroutine zero_ele_kicks

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_has_offset (ele) result (has_offset)
!
! Function to tell if an element has a non-zero offset, pitch or tilt.
! Also see: zero_ele_offsets, zero_ele_kicks, ele_has_nonzero_kick
!
! Modules needed:
!   use bmad
!
! Input
!   ele -- Ele_struct: Element with possible nonzero offsets.
!
! Output:
!   has_offset -- Logical: Set true is element has a non-zero offset.
!-

function ele_has_offset (ele) result (has_offset)

type (ele_struct) ele
logical has_offset

!

has_offset = .false.
if (.not. has_orientation_attributes(ele)) return

select case (ele%key)
case (sbend$)
  if (ele%value(roll_tot$) /= 0) has_offset = .true.
  if (ele%value(ref_tilt_tot$) /= 0) has_offset = .true.
case (mirror$, multilayer_mirror$, crystal$)
  if (ele%value(tilt_tot$) /= 0) has_offset = .true.
  if (ele%value(ref_tilt_tot$) /= 0) has_offset = .true.
case default
  if (ele%value(tilt_tot$) /= 0) has_offset = .true.
end select

if (ele%value(x_pitch_tot$) /= 0) has_offset = .true.
if (ele%value(y_pitch_tot$) /= 0) has_offset = .true.
if (ele%value(x_offset_tot$) /= 0) has_offset = .true.
if (ele%value(y_offset_tot$) /= 0) has_offset = .true.
if (ele%value(z_offset_tot$) /= 0) has_offset = .true.

end function ele_has_offset

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine zero_ele_offsets (ele)
!
! Subroutine to zero the offsets, pitches, tilt and ref_tilt of an element.
! Also see: ele_has_offset, zero_ele_kicks, ele_has_nonzero_kick
!
! Modules needed:
!   use bmad
!
! Input
!   ele -- Ele_struct: Element with possible nonzero offsets, etc.
!
! Output:
!   ele -- Ele_struct: Element with no (mis)orientation.
!-

subroutine zero_ele_offsets (ele)

type (ele_struct) ele

!

if (.not. has_orientation_attributes(ele)) return

select case (ele%key)
case (sbend$)
  ele%value(roll$) = 0
  ele%value(roll_tot$) = 0
  ele%value(ref_tilt$) = 0
  ele%value(ref_tilt_tot$) = 0
case (mirror$, multilayer_mirror$, crystal$)
  ele%value(tilt$) = 0
  ele%value(tilt_tot$) = 0
  ele%value(ref_tilt$) = 0
  ele%value(ref_tilt_tot$) = 0
case default
  ele%value(tilt$) = 0
  ele%value(tilt_tot$) = 0
end select

ele%value(x_pitch$) = 0
ele%value(y_pitch$) = 0
ele%value(x_offset$) = 0
ele%value(y_offset$) = 0
ele%value(z_offset$) = 0

ele%value(x_pitch_tot$) = 0
ele%value(y_pitch_tot$) = 0
ele%value(x_offset_tot$) = 0
ele%value(y_offset_tot$) = 0
ele%value(z_offset_tot$) = 0

end subroutine zero_ele_offsets

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+ 
! Function equivalent_taylor_attributes (ele_taylor, ele2) result (equiv)
!
! Subroutine to see if two elements are equivalent in terms of attributes so
! that their Taylor Maps would be the same. 
!
! This routine is used to see if a taylor map from one element may be 
! used for another and thus save some computation time. elements of type taylor
! are considered *never* to be equivalent since their maps are never computed.
!
! Modules needed:
!   use bmad
!
! Input: 
!   ele_taylor -- Ele_struct: Element with a Taylor map
!   ele2       -- Ele_struct: Element that might receive the Taylor map from ele_taylor.
!
! Output:
!   equiv -- logical: True if elements are equivalent.
!-

function equivalent_taylor_attributes (ele_taylor, ele2) result (equiv)

type (ele_struct) :: ele_taylor, ele2

integer it

logical equiv
logical vmask(num_ele_attrib$), vnot(num_ele_attrib$)

!

equiv = .false.

if (ele_taylor%key == taylor$) return  
if (ele_taylor%key /= ele2%key) return
if (ele_taylor%sub_key /= ele2%sub_key) return
if (ele_taylor%taylor_map_includes_offsets .neqv. ele2%taylor_map_includes_offsets) return
if (ele_taylor%value(integrator_order$) /= ele2%value(integrator_order$)) return

vmask = .true.
vmask(delta_ref_time$) = .false.
vmask(ref_time_start$) = .false.
if ((ele_taylor%key == wiggler$ .or. ele_taylor%key == undulator$) .and. ele_taylor%sub_key == map_type$) then
  vmask( [k1$, rho$, b_max$] ) = .false.  ! These are dependent attributes.
endif
if (.not. ele_taylor%taylor_map_includes_offsets) then
  vmask( [x_offset$, y_offset$, z_offset$, tilt$, x_pitch$, &
            y_pitch$, x_offset_tot$, y_offset_tot$, z_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$] ) = .false.
endif

vnot = (ele_taylor%value /= ele2%value)
vnot = vnot .and. vmask
if (any(vnot)) return

if (associated(ele_taylor%cartesian_map) .neqv. associated(ele2%cartesian_map)) return
if (associated(ele_taylor%cartesian_map)) then
  if (.not. all(ele_taylor%cartesian_map == ele2%cartesian_map)) return
endif

equiv = .true.

end function equivalent_taylor_attributes 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine clear_lat_1turn_mats (lat)
!
! Subroutine to clear the 1-turn matrices in the lat structure:
!   lat%param%t1_no_RF
!   lat%param%t1_with_RF
! This will force any routine dependent upon these to do a remake.
!
! Modules needed:
!   use bmad
!
! Output:
!   lat -- lat_struct: Lat with 1-turn matrices cleared.
!-

subroutine clear_lat_1turn_mats (lat)

type (lat_struct) lat

lat%param%t1_no_RF = 0
lat%param%t1_with_RF = 0

end subroutine clear_lat_1turn_mats

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine match_ele_to_mat6 (ele, err_flag, twiss_ele)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of a match element. 
!
! Note: ele%taylor%term will be deallocated if the xfer map has changed.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele       -- ele_struct: Match element.
!   twiss_ele -- ele_struct, optional: If present and ele%value(match_end$) is True, 
!                  then use the Twiss parameters in this element as the upstream end Twiss.
!
! Output:
!   ele       -- ele_struct: Match element.
!     %vec0(6)   -- Real(rp): 0th order part of the transfer map.
!     %mat6(6,6) -- Real(rp): Transfer matrix (1st order part of xfer map).
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!-

subroutine match_ele_to_mat6 (ele, err_flag, twiss_ele)

type (ele_struct), target :: ele, ele0, ele1
type (ele_struct), optional, target :: twiss_ele
type (ele_struct), pointer :: t_ele

real(rp), pointer :: v(:)
real(rp) orb0(6), orb1(6)
real(rp) old_mat6(6,6), old_vec0(6)
logical err_flag

! Error Check

old_mat6 = ele%mat6
old_vec0 = ele%vec0

if (ele%value(beta_a1$) == 0 .or. ele%value(beta_b1$) == 0) then
  ele%mat6 = 0
  ele%vec0 = 0
  err_flag = .true.
  call kill_taylor (ele%taylor)
  return
endif

! Match_end

if (is_true(ele%value(match_end$))) then
  if (present(twiss_ele)) then
    t_ele => twiss_ele
  else
    t_ele => pointer_to_next_ele (ele, -1)
  endif
  ele%value(beta_a0$)    = t_ele%a%beta
  ele%value(beta_b0$)    = t_ele%b%beta
  ele%value(alpha_a0$)   = t_ele%a%alpha
  ele%value(alpha_b0$)   = t_ele%b%alpha
  ele%value(eta_x0$)     = t_ele%x%eta
  ele%value(eta_y0$)     = t_ele%y%eta
  ele%value(etap_x0$)    = t_ele%x%etap
  ele%value(etap_y0$)    = t_ele%y%etap
  ele%value(gamma_c$)    = t_ele%gamma_c
  ele%value(c_11$:c_22$) = [t_ele%c_mat(1,1), t_ele%c_mat(1,2), t_ele%c_mat(2,1), t_ele%c_mat(2,2)]
endif

! Special case where match_end is set but there is no beginning beta value yet.
! In this case, just return the unit matrix. 
! This is not an error since it is important for lat_make_mat6 to keep on computing matrices.

if (is_true(ele%value(match_end$)) .and. (ele%value(beta_a0$) == 0 .or. ele%value(beta_b0$) == 0)) then
  call mat_make_unit (ele%mat6)
  ele%vec0 = 0
  err_flag = .false.
  call kill_taylor (ele%taylor)
  return
endif

!

if (ele%value(beta_a0$) == 0 .or. ele%value(beta_b0$) == 0) then
  ele%mat6 = 0
  err_flag = .true.
  call kill_taylor (ele%taylor)
  return
endif

!

err_flag = .false.

v => ele%value

ele0%a%beta   = v(beta_a0$)
ele0%a%alpha  = v(alpha_a0$)
ele0%a%phi    = 0
ele0%x%eta    = v(eta_x0$)
ele0%x%etap   = v(etap_x0$)

ele0%b%beta   = v(beta_b0$)
ele0%b%alpha  = v(alpha_b0$)
ele0%b%phi    = 0
ele0%y%eta    = v(eta_y0$)
ele0%y%etap   = v(etap_y0$)

ele1%a%beta   = v(beta_a1$)
ele1%a%alpha  = v(alpha_a1$)
ele1%a%phi    = v(dphi_a$)
ele1%x%eta    = v(eta_x1$)
ele1%x%etap   = v(etap_x1$)

ele1%b%beta   = v(beta_b1$)
ele1%b%alpha  = v(alpha_b1$)
ele1%b%phi    = v(dphi_b$)
ele1%y%eta    = v(eta_y1$)
ele1%y%etap   = v(etap_y1$)

ele0%c_mat(1,:) = [v(c_11$), v(c_12$)]
ele0%c_mat(2,:) = [v(c_21$), v(c_22$)]
ele0%gamma_c    = v(gamma_c$)

ele1%c_mat = 0 
ele1%gamma_c = 1

ele0%name = ele%name
ele1%name = ele%name

orb0 = [v(x0$), v(px0$), v(y0$), v(py0$), v(z0$), v(pz0$)]
orb1 = [v(x1$), v(px1$), v(y1$), v(py1$), v(z1$), v(pz1$)]

call transfer_mat_from_twiss (ele0, ele1, orb0, orb1, ele%mat6)

! Kick part

ele%vec0 = orb1 - matmul (ele%mat6, orb0)

if (any(ele%mat6 /= old_mat6) .or. any(ele%vec0 /= old_vec0)) call kill_taylor (ele%taylor)

end subroutine match_ele_to_mat6

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine calc_super_slave_key (lord1, lord2, slave, create_jumbo_slave)
!
! Function to decide what the type of element a super_slave should be based
! upon the types of its lords.
!
! Modules needed:
!   use bmad
!
! Input:
!   lord1     -- Ele_struct: First slave.
!     %key
!     %sub_key
!   lord2     -- Ele_struct: Second slave.
!     %key
!     %sub_key
!   create_jumbo_slave
!              -- Logical, optional: If True then slave%key will be set to em_field.
!                   Default is False. 
!
! Output:
!   slave -- Ele_struct: Super_slave element.
!     %key        -- Set to -1 if there is an error.
!     %sub_key
!-

subroutine calc_super_slave_key (lord1, lord2, slave, create_jumbo_slave)

type (ele_struct), target :: lord1, lord2, slave
integer key1, key2
logical, optional :: create_jumbo_slave

!

key1 = lord1%key
key2 = lord2%key

slave%key = -1  ! Default if no superimpose possible
slave%sub_key = 0

! control elements, etc. cannot be superimposed.

select case (key1)
case (overlay$, group$, girder$, taylor$, match$, patch$, fiducial$, floor_shift$, multipole$, ab_multipole$)
  return
end select

select case (key2)
case (overlay$, group$, girder$, taylor$, match$, patch$, fiducial$, floor_shift$, multipole$, ab_multipole$)
  return
end select

! Superimposing two of like kind...

if (key1 == key2) then
  select case (key1)
  case (sbend$)
    ! Bad
  case (wiggler$, undulator$, rfcavity$)
    slave%key = em_field$
    slave%sub_key = const_ref_energy$
  case (lcavity$)
    slave%key = em_field$
    slave%sub_key = nonconst_ref_energy$
  case (em_field$)
    slave%key = em_field$
    if (lord1%sub_key == nonconst_ref_energy$ .or. lord2%sub_key == nonconst_ref_energy$) then
      slave%sub_key = nonconst_ref_energy$
    else
      slave%sub_key = const_ref_energy$
    endif
  case default
    slave%key = key1
  end select
  return
endif

! If one element is a drift then slave%key = key of other element.

if (key1 == drift$) then
  slave%key = key2
  slave%sub_key = lord2%sub_key
  return
endif

if (key2 == drift$) then
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
endif

! sbend elements are problematical due to the different reference orbit so cannot superimpose them.

if (key1 == sbend$ .or. key2 == sbend$) return

! If one element is a pipe then slave%key = key of other element.

if (any(key1 == [pipe$])) then
  slave%key = key2
  slave%sub_key = lord2%sub_key
  return
endif

if (any(key2 == [pipe$])) then
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
endif

! If one element is a rcollimator, monitor, or instrument then slave%key = key of other element.

if (lord1%aperture_type == elliptical$ .or. lord2%aperture_type == elliptical$ ) slave%aperture_type = elliptical$

if (any(key1 == [ecollimator$, rcollimator$, monitor$, instrument$])) then
  slave%key = key2
  slave%sub_key = lord2%sub_key
  return
endif

if (any(key2 == [ecollimator$, rcollimator$, monitor$, instrument$])) then
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
endif

! If one element is a kicker then slave%key = key of other element.

if (any(key1 == [kicker$, hkicker$, vkicker$])) then
  if (any(key2 == [kicker$, hkicker$, vkicker$])) then
    slave%key = kicker$
  else
    slave%key = key2
  endif
  return
endif

if (any(key2 == [kicker$, hkicker$, vkicker$])) then
  slave%key = key1
  slave%sub_key = lord1%sub_key
  return
endif

! General case...

! em_field wanted

if (logic_option(.false., create_jumbo_slave)) then
  slave%key = em_field$
  if (key1 == lcavity$ .or. key2 == lcavity$) then
    slave%sub_key = nonconst_ref_energy$
  elseif (key1 == em_field$) then
    slave%sub_key = lord1%sub_key
  elseif (key2 == em_field$) then
    slave%sub_key = lord2%sub_key
  else
    slave%sub_key = const_ref_energy$
  endif
  return
endif

!

select case (key1)

case (quadrupole$,  solenoid$, sol_quad$) 
  select case (key2)
  case (quadrupole$);    slave%key = sol_quad$
  case (solenoid$);      slave%key = sol_quad$
  case (sol_quad$);      slave%key = sol_quad$
  case (bend_sol_quad$); slave%key = bend_sol_quad$
  case (sbend$);         slave%key = bend_sol_quad$
  end select

case (bend_sol_quad$)
  select case (key2)
  case (quadrupole$);    slave%key = bend_sol_quad$
  case (solenoid$);      slave%key = bend_sol_quad$
  case (sol_quad$);      slave%key = bend_sol_quad$
  case (sbend$);         slave%key = bend_sol_quad$
  end select
end select

if (slave%key /= -1) return  ! Have found something

! Only thing left is to use em_field type element.

slave%key = em_field$
if (key1 == lcavity$ .or. key2 == lcavity$) then
  slave%sub_key = nonconst_ref_energy$
elseif (key1 == em_field$) then
  slave%sub_key = lord1%sub_key
elseif (key2 == em_field$) then
  slave%sub_key = lord2%sub_key
else
  slave%sub_key = const_ref_energy$
endif

end subroutine calc_super_slave_key

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function gradient_shift_sr_wake (ele, param) result (grad_shift)
! 
! Function to return the shift in the accelerating gradient due to:
!   1) Short range longitudinal wake forces.
!   2) Adjustment in the external applied RF power attempting to compensate for the SR wake loss.
!
! Note: This routine will only return a non-zero value when bmad_com%sr_wakes_on = True
!
! Module needed:
!   use bmad
!
! Input:
!   ele           -- ele_struct: Lcavity element.
!   param         -- lat_param_struct: Lattice parameters
!     %n_part        -- Number of particles in a bunch
!     %particle      -- Type of particle
!
! Output:
!   grad_shift -- Real(rp): Shift in gradient
!-

function gradient_shift_sr_wake (ele, param) result (grad_shift)

type (ele_struct) ele
type (lat_param_struct) param
real(rp) grad_shift

! 

if (bmad_com%sr_wakes_on .and. ele%value(l$) /= 0) then
  grad_shift = ele%value(e_loss$) * param%n_part * abs(charge_of(param%particle)) * &
                                                                  e_charge / ele%value(l$) 
else
  grad_shift = 0
endif

end function gradient_shift_sr_wake

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_ele_status_stale (ele, status_group, set_slaves)
!
! Routine to set a status flags to stale in an element and the corresponding 
! ones for any slaves the element has.
!
! Also the branch%param structure of the branch the element is in is set.
!
! For example: status_group = ref_energy_group$ sets stale:
!   ele%bookkeeping_state%ref_energy 
!   ele%bookkeeping_state%floor_position
!   ele%bookkeeping_state%mat6
! See the code for more details.
! 
! Output:
!   ele           -- ele_struct: Element.
!     %bookkeeping_state   -- Status block to set.
!   status_group  -- Integer: Which flag groups to set. Possibilities are:
!                      attribute_group$, control_group$, floor_position_group$, s_position_group$, 
!                      s_and_floor_position_group$, ref_energy_group$, or mat6_group$, all_groups$
!   set_slaves    -- Logical, optional: If present and False then do not set
!                      the status for any slaves. Default is True.
!-

recursive subroutine set_ele_status_stale (ele, status_group, set_slaves)

type (bookkeeping_state_struct), pointer :: state
type (ele_struct), target :: ele
type (ele_struct), pointer :: slave
integer status_group, i
logical, optional :: set_slaves

! Only set overall lattice status flags if the element is part of a lattice.


if (ele%ix_ele > -1 .and. associated(ele%branch)) then
  ! If a lord
  if (ele%ix_branch == 0 .and. ele%ix_ele > ele%branch%n_ele_track) then
    state => ele%branch%lat%lord_state
  else
    state => ele%branch%param%bookkeeping_state
  endif
else
  nullify(state)
endif

!

select case (status_group)

case (attribute_group$)
  call set_attributes_status
  call set_mat6_status
  call set_ptc_status

case (control_group$)
  call set_control_status

case (floor_position_group$)
  call set_floor_position_status
  call set_mat6_status
  call set_ptc_status

case (s_and_floor_position_group$)
  call set_s_position_status
  call set_floor_position_status
  call set_mat6_status
  call set_ptc_status

case (s_position_group$)
  call set_s_position_status
  call set_mat6_status
  call set_ptc_status

case (ref_energy_group$)
  call set_ref_energy_status
  call set_mat6_status
  call set_attributes_status ! EG: k1 <--> b1_gradient calc needed 
  call set_ptc_status

case (mat6_group$)
  call set_mat6_status

case (rad_int_group$)
  call set_rad_int_status

case (all_groups$)
  call set_attributes_status
  call set_control_status
  call set_ref_energy_status
  call set_floor_position_status
  call set_s_position_status
  call set_mat6_status
  call set_rad_int_status
  call set_ptc_status

case default
   if (global_com%exit_on_error) call err_exit   ! Should not be here

end select

! Set slave

if (logic_option(.true., set_slaves)) then
  do i = 1, ele%n_slave
    slave => pointer_to_slave (ele, i)
    call set_ele_status_stale (slave, status_group)
  enddo
endif

!----------------------------------------------------------------------------
contains

subroutine set_attributes_status
  if (ele%key == overlay$) return
  if (ele%key == group$) return
  ele%bookkeeping_state%attributes = stale$
  if (associated(state)) state%attributes = stale$
end subroutine set_attributes_status

!----------------------------------------------------------------------------
! contains

subroutine set_control_status
  if (ele%lord_status == not_a_lord$ .and. ele%n_lord == 0) return
  ele%bookkeeping_state%control = stale$
  if (associated(state)) state%control = stale$
end subroutine set_control_status

!----------------------------------------------------------------------------
! contains

subroutine set_floor_position_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%floor_position = stale$
  if (associated(state)) state%floor_position = stale$
end subroutine set_floor_position_status

!----------------------------------------------------------------------------
! contains

subroutine set_s_position_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%s_position = stale$
  if (associated(state)) state%s_position = stale$
end subroutine set_s_position_status

!----------------------------------------------------------------------------
! contains

subroutine set_ref_energy_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%ref_energy = stale$
  if (associated(state)) state%ref_energy = stale$
end subroutine set_ref_energy_status

!----------------------------------------------------------------------------
! contains

subroutine set_rad_int_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%rad_int = stale$
  if (associated(state)) state%rad_int = stale$
end subroutine set_rad_int_status

!----------------------------------------------------------------------------
! contains

subroutine set_ptc_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%ptc = stale$
  if (associated(state)) state%ptc = stale$
end subroutine set_ptc_status

!----------------------------------------------------------------------------
! contains

! Ignore if the element does not have an associated linear transfer map
! Also the branch status does not get set since the transfer map calc
! must always check a branch due to possible reference orbit shifts.

subroutine set_mat6_status
  if (ele%key == overlay$) return
  if (ele%key == group$) return
  if (ele%lord_status == multipass_lord$) return
  ele%bookkeeping_state%mat6 = stale$
end subroutine set_mat6_status

end subroutine set_ele_status_stale 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_lords_status_stale (ele, stat_group, flag)
!
! Routine to recursively set the status flag of all slaves of an element.
!
! Input:
!   ele        -- ele_struct: Element
!   stat_group -- Integer: which status group to set. floor_position_group$, etc.
!                   See set_ele_status_stale for more details.
!   flag       -- Logical, optional: Do not use. For determining recursion depth.
!
! Output:
!   ele%lat    -- Lat_struct: Lattice with status flags of lords of ele set.
!-

recursive subroutine set_lords_status_stale (ele, stat_group, flag)

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer stat_group, i
logical, optional :: flag

! First time through the flag argument will not be present.
! Do not set status first time through since this is the original element.
! That is, only want to set the flags of the lords.

if (present(flag)) call set_ele_status_stale (ele, stat_group, .false.)

do i = 1, ele%n_lord
  lord => pointer_to_lord (ele, i)
  call set_lords_status_stale (lord, stat_group, .true.)
enddo

end subroutine set_lords_status_stale

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0) result (branch_ptr)
!
! Function to point to the named lattice branch.
! This routine is overloaded by the routine: pointer_to_branch.
! See pointer_to_branch for more details.
! 
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   branch_name -- Character(*): May be a branch name or a branch index.
!   lat         -- Lat_struct: Lattice to search.
!   parameter_is_branch0 -- logical, optional: If True, 'PARAMETER' is taken to be
!                     an alternative name for branch(0). Default is False.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the nameed branch.
!                   Nullified if there is no such branch.
!-

function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0) result (branch_ptr)

type (branch_struct), pointer :: branch_ptr
type (lat_struct), target :: lat

integer ib, ios
logical, optional :: parameter_is_branch0

character(*) branch_name
character(40) b_name
character(*), parameter :: r_name = 'pointer_to_branch_given_name'

! parameter

call str_upcase (b_name, branch_name)

if (logic_option(.false., parameter_is_branch0) .and. b_name == 'PARAMETER') then
  branch_ptr => lat%branch(0)
  return
endif

! Init in case of error

nullify(branch_ptr)

! Is index.

if (is_integer(trim(b_name))) then
  read (b_name, *, iostat = ios) ib
  if (ib < 0 .or. ib > ubound(lat%branch, 1)) then
    !! call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: ' // b_name)
    return
  endif
  branch_ptr => lat%branch(ib)

! Is name.

else
  do ib = lbound(lat%branch, 1), ubound(lat%branch, 1)
    if (lat%branch(ib)%name == b_name) then
      branch_ptr => lat%branch(ib)
      return
    endif
  enddo
  !! call out_io (s_error$, r_name, 'BRANCH NAME NOT FOUND: ' // branch_name)
  return
endif

end function pointer_to_branch_given_name

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch_given_ele (ele) result (branch_ptr)
!
! Function to point to the lattice branch associated with an element.
! This routine is overloaded by the routine: pointer_to_branch.
! See pointer_to_branch for more details.
!
! Note: Result is ambiguous if the element is associated with multiple branches which
! can happen, for example, with overlay lord elements.
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   ele      -- Ele_struct: Element.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the branch associated with the element.
!                   Nullified if the element is not associated with a lattice.
!-

recursive function pointer_to_branch_given_ele (ele) result (branch_ptr)

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch_ptr

! Now associated with a lattice case

if (.not. associated(ele%branch)) then
  nullify(branch_ptr)
  return
endif

! Not a lord case

if (ele%n_slave == 0) then
  branch_ptr => ele%branch
  return
endif

! If a lord then look to the first slave for the associated branch

branch_ptr => pointer_to_branch_given_ele(pointer_to_slave(ele, 1))

end function pointer_to_branch_given_ele

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine remove_lord_slave_link (lord, slave)
!
! Routine to remove the pointers between a lord and a slave.
! Note: This routine will not modify lord%lord_status and slave%slave_status.
!
! Input:
!   lord  -- ele_struct: Lord element
!   slave -- ele_struct: Slave element
!
! Output:
!   lord  -- ele_struct: Lord element with link info removed
!   slave -- ele_struct: Slave element with link info removed
!   lat   -- Lattice_struct: Lattice obtaind from lord%lat and slave%lat pointers.
!     %control(:)  -- Array modified to remove lord/slave link.
!     %ic(:)       -- Array modified to remove lord/slave link.
!-

subroutine remove_lord_slave_link (lord, slave)

type (ele_struct), target :: lord, slave
type (ele_struct), pointer :: ele
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch

integer ic_out, icon_out, ib, i, n, ic2_lord

!

if (.not. associated(lord%branch, slave%branch)) call err_exit  ! Should not be

lat => lord%branch%lat

! Find lat%control(:) and lat%ic(:) elements associated with this link

ic2_lord = slave%ic1_lord + slave%n_lord - 1
do ic_out = slave%ic1_lord, ic2_lord
  icon_out = lat%ic(ic_out)
  if (lat%control(icon_out)%lord%ix_ele == lord%ix_ele) exit
enddo

if (icon_out == ic2_lord+1) call err_exit ! Should not be

! Compress lat%control and lat%ic arrays.

n = lat%n_control_max
lat%control(icon_out:n-1) = lat%control(icon_out+1:n)
lat%n_control_max = n - 1

n = lat%n_ic_max
lat%ic(ic_out:n-1) = lat%ic(ic_out+1:n)
lat%n_ic_max = n - 1

do i = 1, n - 1
  if (lat%ic(i) > icon_out) lat%ic(i) = lat%ic(i) - 1
enddo

! Correct info in all elements

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)

    if (ele%ix1_slave >  icon_out) ele%ix1_slave = ele%ix1_slave - 1
    if (ele%ic1_lord >  ic_out) ele%ic1_lord = ele%ic1_lord - 1
  enddo
enddo

slave%n_lord = slave%n_lord - 1
if (slave%n_lord == 0) then
  slave%ic1_lord = 0
  slave%slave_status = not_a_child$
endif

lord%n_slave = lord%n_slave - 1
if (lord%n_slave == 0) then
  lord%ix1_slave = 0
  lord%lord_status = not_a_lord$
endif

end subroutine remove_lord_slave_link 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_next_ele (this_ele, offset, skip_beginning, follow_fork) result (next_ele)
!
! Function to return a pointer to the N^th element relative to this_ele
! in the array of elements in a lattice branch.
! 
! If the this_ele is a super_lord element, the appropriate element in the tracking 
! part of the lattice is returned.
!
! If this_ele is a lord element but not a super_lord then this is an error since
! it is not clear what to do in this case.
!
! This routine will always wrap around between branch end and branch beginning
! independent of the branch geometry. That is, it is the responsibility of the
! calling routine to check this if important.
!
! Notice that the first element in a lattice is the beginning element with index 0.
!
! Input:
!   this_ele       -- ele_struct: Starting element.
!   offset         -- integer, optional: +1 -> return next element, +2 -> element 
!                       after that, etc. Can be negative. Default = +1.
!   skip_beginning -- logical, optional: If True then skip beginning element #0
!                       when wrapping around. Default is False.
!   follow_fork    -- logical, optional: If True then fork at any fork element.
!
!   next_ele -- ele_struct, pointer: Element after this_ele (if offset = 1).
!                Nullified if there is an error. EG bad this_ele.
!-

function pointer_to_next_ele (this_ele, offset, skip_beginning, follow_fork) result (next_ele)

type (ele_struct), target :: this_ele
type (ele_struct), pointer :: next_ele
type (ele_struct), pointer :: an_ele
type (branch_struct), pointer :: branch

integer, optional :: offset
integer i, ix_ele, n_off
logical, optional :: skip_beginning, follow_fork

! If element is not associated with a lat_struct then do nothing.

next_ele => null()

if (.not. associated(this_ele%branch)) return   ! Not associated with a lat_struct
if (this_ele%ix_ele < 0 .or. this_ele%ix_ele > this_ele%branch%n_ele_max) return  ! Slice_slave?

! Compute offset

n_off = integer_option(+1, offset)
if (n_off == 0) then
  next_ele => this_ele
  return
endif

! If a lord then point to a slave. 
! Initially point to the first or last slave element so that this routine
! will return a pointer to an element that is not a slave of this_ele.

if (this_ele%ix_ele > this_ele%branch%n_ele_track) then  ! Is a lord
  if (this_ele%lord_status /= super_lord$) return    ! Error
  if (n_off > 0) then
    an_ele => pointer_to_slave(this_ele, this_ele%n_slave)
  else
    an_ele => pointer_to_slave(this_ele, 1)
  endif
else
  an_ele => this_ele
endif

! Apply offset.
! If follow_fork = True then must check all elements in between for a possible fork element.

ix_ele = an_ele%ix_ele
branch => an_ele%branch

if (n_off > 0) then

  do i = 1, abs(n_off)

    if (logic_option(.false., follow_fork) .and. (an_ele%key == fork$ .or. an_ele%key == photon_fork$)) then
      ix_ele = nint(an_ele%value(ix_to_element$))
      branch => an_ele%branch%lat%branch(nint(an_ele%value(ix_to_branch$)))
    else
      ix_ele = ix_ele + 1
    endif

    if (ix_ele == branch%n_ele_track + 1) then
      if (logic_option(.false., skip_beginning)) then
        ix_ele = 1
      else
        ix_ele = 0
      endif
    endif

    an_ele => branch%ele(ix_ele)

  enddo

! Going backwards

else
  do i = 1, abs(n_off)

    if (logic_option(.false., follow_fork) .and. nint(an_ele%value(direction$)) == -1 .and. & 
                            (an_ele%key == fork$ .or. an_ele%key == photon_fork$)) then
      ix_ele = nint(an_ele%value(ix_to_element$))
      branch => an_ele%branch%lat%branch(nint(an_ele%value(ix_to_branch$)))
    else
      ix_ele = ix_ele - 1
    endif

    if (logic_option(.false., skip_beginning)) then
      if (ix_ele == 0) ix_ele = branch%n_ele_track
    else
      if (ix_ele == -1) ix_ele = branch%n_ele_track
    endif

    an_ele => branch%ele(ix_ele)

  enddo
endif

next_ele => an_ele

end function pointer_to_next_ele

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
!
! Function to return a pointer to an element in a lattice.
! This routine is overloaded by pointer_to_ele.
! See pointer_to_ele for more details.
!-

function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_ptr

integer ix_branch, ix_ele

!

ele_ptr => null()

if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) return
if (ix_ele < 0 .or. ix_ele > lat%branch(ix_branch)%n_ele_max) return

ele_ptr => lat%branch(ix_branch)%ele(ix_ele)

end function pointer_to_ele1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
!
! Function to return a pointer to an element in a lattice.
! This routine is overloaded by pointer_to_ele.
! See pointer_to_ele for more details.
!-

function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_ptr
type (lat_ele_loc_struct) ele_loc

!

ele_ptr => null()

if (ele_loc%ix_branch < 0 .or. ele_loc%ix_branch > ubound(lat%branch, 1)) return
if (ele_loc%ix_ele < 0 .or. ele_loc%ix_ele > lat%branch(ele_loc%ix_branch)%n_ele_max) return

ele_ptr => lat%branch(ele_loc%ix_branch)%ele(ele_loc%ix_ele)

end function pointer_to_ele2

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_has_constant_ds_dt_ref (ele) result (is_const)
!
! Function to determine if an element has a constant longitudinal reference velocity.
! When in doubt, the assumption is that the longitudinal velocity is not constant.
!
! Module needed:
!   use bmad
!
! Input:
!   ele -- ele_struct: Element.
!
! Output:
!   is_const -- Logical: True if reference velocity must be a constant.
!-

function ele_has_constant_ds_dt_ref (ele) result (is_const)

type (ele_struct) ele
logical is_const

! Anything with longitudinal electric fields or anything
! where the "zero-orbit" is not a straight line down the middle
! has a varying ds/dt(ref).

select case (ele%key)
case (lcavity$, custom$, hybrid$, wiggler$, undulator$, rfcavity$, em_field$)
  is_const = .false.
case default
  is_const = .true.
end select

end function ele_has_constant_ds_dt_ref

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function tracking_uses_end_drifts (ele, use_hard_edge_model) result (has_drifts)
!
! Function to determine if the tracking for an element uses a "hard edge model"
! where the tracking looks like (drift, model, drift). For example,
! RF cavity fields with ele%field_calc = bmad_standard$ use a hard edge model where
! the length of the cavity is c_light / (2 * freq).
!
! If the use_hard_edge model argument is set to False, this function will always return False.
!
! Module needed:
!   use bmad
!
! Input:
!   ele                   -- ele_struct: Element.
!   use_hard_edge_drifts  -- logical, optional: Use bmad hard edge model for cavities, etc?
!                              Default is set by bmad_com%use_hard_edge_drifts.
!                              Default bmad_com%use_hard_edge_drifts is True.
!
! Output:
!   has_drifts -- Logical: True if tracking uses end drifts.
!-

function tracking_uses_end_drifts (ele, use_hard_edge_model) result (has_drifts)

type (ele_struct) ele
logical has_drifts
logical, optional :: use_hard_edge_model

!

has_drifts = .false.
if (.not. logic_option(bmad_com%use_hard_edge_drifts, use_hard_edge_model)) return

select case (ele%key)
case (lcavity$, rfcavity$, solenoid$)
    if (ele%field_calc == bmad_standard$) has_drifts = .true.
    if (ele%value(l_hard_edge$) == ele%value(l$)) has_drifts = .false.
end select

end function tracking_uses_end_drifts

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine create_hard_edge_drift (ele_in, which_end, drift_ele)
!
! Routine to create the drift element for the end drifts of an element.
! The end drifts are present, for example, when doing runge_kutta tracking
! through an rf_cavity with field_calc = bmad_standard. In this case, the
! field model is a pi-wave hard-edge resonator whose length may not match
! the length of the element and so particles must be drifted from the edge
! of the element to the edge of the field model.
!
! Input:
!   ele_in     -- ele_struct: Input element.
!   which_end  -- Integer: Which end is being created. upstream_end$ or downstream_end$.
!                   For an Lcavity one can have differences in reference energy.
!
! Output:
!   drift_ele  -- ele_struct: drift elment.
!-

subroutine create_hard_edge_drift (ele_in, which_end, drift_ele)

type (ele_struct) ele_in, drift_ele
real(rp) E_tot, p0c
integer which_end

!

select case (which_end)
case (upstream_end$)
  e_tot = ele_in%value(e_tot_start$)
  p0c   = ele_in%value(p0c_start$)
  drift_ele%name                   = 'drift1_' // ele_in%name(1:33)
case (downstream_end$)
  e_tot = ele_in%value(e_tot$)
  p0c   = ele_in%value(p0c$)
  drift_ele%name                   = 'drift2_' // ele_in%name(1:33)
case default
  if (global_com%exit_on_error) call err_exit
end select

drift_ele%key                    = drift$

drift_ele%value                  = 0
drift_ele%value(p0c$)            = p0c
drift_ele%value(e_tot$)          = e_tot
drift_ele%value(p0c_start$)      = p0c
drift_ele%value(e_tot_start$)    = e_tot
drift_ele%value(l$)              = (ele_in%value(l$) - hard_edge_model_length(ele_in)) / 2 
drift_ele%value(ds_step$)        = drift_ele%value(l$)
drift_ele%value(num_steps$)      = 1
drift_ele%value(delta_ref_time$) = drift_ele%value(l$) * e_tot / (c_light * p0c)
drift_ele%orientation            = ele_in%orientation

end subroutine create_hard_edge_drift 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function hard_edge_model_length (ele) result (l_hard)
!
! Input:
!   ele -- ele_struct: Element
!
! Output:
!   l_hard -- real(rp): Length of the hard edge model.
!-

function hard_edge_model_length (ele) result (l_hard)

type (ele_struct) ele
real(rp) l_hard

!

if (bmad_com%use_hard_edge_drifts) then
  l_hard = ele%value(l_hard_edge$)
else
  l_hard = ele%value(l$)
endif

end function hard_edge_model_length

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function absolute_time_tracking (ele) result (is_abs_time)
!
! Routine to return a logical indicating whether the tracking through an
! element should use absolute time or time relative to the reference particle.
!
! Note: e_gun elements always use aboslute time tracking to get around
! the problem when the particle velocity is zero.
!
! Input:
!   ele  -- ele_struct: Element being tracked through.
!
! Output:
!   is_abs_time -- Logical: True if absolute time tracking is needed.
!-

function absolute_time_tracking (ele) result (is_abs_time)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
logical is_abs_time
integer i

!

is_abs_time = bmad_com%absolute_time_tracking_default
if (associated(ele%branch)) is_abs_time = ele%branch%lat%absolute_time_tracking
if (ele%key == e_gun$) is_abs_time = .true.

if (ele%key == em_field$ .and. (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$)) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%key == e_gun$) is_abs_time = .true.
  enddo
endif

end function absolute_time_tracking

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
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

type (coord_struct) orbit
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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function rf_is_on (branch) result (is_on)
!
! Routine to check if any rfcavity is powered in a branch.
!
! Input:
!   branch -- branch_struct: Lattice branch to check.
!
! Output:
!   is_on  -- Logical: True if any rfcavity is powered. False otherwise.
!-

function rf_is_on (branch) result (is_on)

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele

integer i
logical is_on

!

is_on = .false.

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%key == rfcavity$ .and. ele%is_on .and. ele%value(voltage$) /= 0) then
    is_on = .true.
    return
  endif
enddo

end function rf_is_on

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine canonical_to_angle_coords (orbit)
!
! Routine to convert from canonical (x, px, y, py, z, pz) coordinates to
! angle (x, x', y, y', z, z') coordinates.
!
! Input:
!   orbit -- coord_struct: Orbit in canonical coordinates.
!
! Output:
!   orbit -- coord_struct: Orbit in angular coordinates.
!-

subroutine canonical_to_angle_coords (orbit)

type (coord_struct) orbit

!

orbit%vec(2) = orbit%vec(2) / (1 + orbit%vec(6))
orbit%vec(4) = orbit%vec(4) / (1 + orbit%vec(6))

end subroutine canonical_to_angle_coords

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine angle_to_canonical_coords (orbit)
!
! Routine to convert from angle (x, x', y, y', z, z') coordinates to
! canonical (x, px, y, py, z, pz) coordinates.
!
! Input:
!   orbit -- coord_struct: Orbit in angular coordinates.
!
! Output:
!   orbit -- coord_struct: Orbit in canonical coordinates.
!-

subroutine angle_to_canonical_coords (orbit)

type (coord_struct) orbit

!

orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))
orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))

end subroutine angle_to_canonical_coords

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_value_has_changed (ele, list, abs_tol, set_old) result (has_changed)
!
! Routine to see if a parameter value in a lattice element has changed significantly.
! A relative tolerance of small_rel_change$ = 1d-14 is also added to abs_tol.
!
! Input:
!   ele         -- ele_struct: Element under consideration.
!   list(:)     -- integer: List of indexes of ele%value(:) array to check.
!   abs_tol(:)  -- real(rp): List of values such that if the change in parameter value is
!                    less than this it is not considered to have changed significantly.
!   set_old     -- logical: If True then set ele%old_value(j) = ele%value(j) for j in list
!
! Output:
!   ele         -- ele_struct: ele%old_value may be set depending upon setting of set_old
!   has_changed -- logical: Set True if a value has changed significantly.
!-

function ele_value_has_changed (ele, list, abs_tol, set_old) result (has_changed)

type (ele_struct) ele
integer list(:)
integer i, j
real(rp) abs_tol(:)
logical set_old, has_changed

!

has_changed = .false.
do i = 1, size(list)
  j = list(i)
  if (abs(ele%value(j) - ele%old_value(j)) > abs_tol(i) + &
                    small_rel_change$ * (abs(ele%value(j) + abs(ele%old_value(j))))) then
    has_changed = .true.
    exit
  endif
enddo

if (set_old) ele%old_value(list) = ele%value(list)

end function ele_value_has_changed

end module
