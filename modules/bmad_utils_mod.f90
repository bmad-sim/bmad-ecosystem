!+
! Module bmad_utils_mod
!
! Module for subroutines that use bmad_struct structures but do not
! call other routines in bmad_interface.
!
! ALSO: THESE ROUTINES DO NOT HAVE ACCESS TO THE OVERLOADED
! ELE1 = ELE2 AND LAT1 = LAT2.
!-

module bmad_utils_mod

use make_mat6_mod
use basic_attribute_mod
use basic_bmad_interface

private pointer_to_ele1, pointer_to_ele2
private pointer_to_branch_given_name, pointer_to_branch_given_ele

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord (...)
!
! Routine to initialize a coord_struct. 
!
! This routine is an overloaded name for:
!   Subroutine init_coord1 (orb, vec, ele, at_downstream_end, particle, E_photon, t_ref_offset, shift_vec6)
!   Subroutine init_coord2 (orb, orb_in, ele, at_downstream_end, particle, t_ref_offset, shift_vec6)
!
! Exception: If ele is an init_ele (branch%ele(0)), orb%p0c is shifted to ele%value(p0c$).
! Additionally, if ele is an init_ele, and vec is zero or not present, orb%vec(6) is shifted
! so that the particle's energy is maintained at ele%value(p0c_start$).
!
! Modules needed:
!   use bmad
!
! Input:
!   orb_in       -- Coord_struct: Input orbit.
!   vec(6)       -- real(rp), optional: Coordinate vector. If not present then taken to be zero.
!   ele          -- ele_struct, optional: Particle is initialized to start from the entrance end of ele
!   at_downstream_end  -- Logical, optional: Particle is at entrance or exit end of the element?
!                     Must be present if ele argument is present.
!                     Default is False.
!   particle     -- Integer, optional: Particle type (electron$, etc.). 
!   E_photon     -- real(rp), optional: Photon energy if particle is a photon. Ignored otherwise.
!   t_ref_offset -- real(rp), optional: Offset of the reference time. This is non-zero when
!                     there are multiple bunches and the reference time for a particular particle
!                     is pegged to the time of the center of the bunch.
!   shift_vec6   -- Logical, optional: If present and False, prevent the shift of orb%vec(6).
!
! Output:
!   orb -- Coord_struct: Initialized coordinate.
!                 Note: For photons, orb%vec(6) is computed as sqrt(1 - vec(2)^2 - vec(4)^2) if needed.
!-

interface init_coord
  module procedure init_coord1
  module procedure init_coord2
end interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine reallocate_coord (...)
!
! Routine to allocate or reallocate at allocatable coord_struct array.
! reallocate_coord is an overloaded name for:
!   reallocate_coord_n (coord, n_coord)
!   reallocate_coord_lat (coord, lat, ix_branch)
!
! Subroutine to allocate an allocatable coord_struct array to at least:
!     coord(0:n_coord)                            if n_coord arg is used.
!     coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.
!
! The old coordinates are saved
! If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
! In any case, coord(n)%vec for n > 0 is set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord(:)  -- Coord_struct, allocatable: Allocatable array.
!   n_coord   -- Integer: Minimum array upper bound wanted.
!   lat       -- lat_struct: Lattice 
!   ix_branch -- Integer, optional: Branch to use. Default is 0 (main branch).
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

interface reallocate_coord
  module procedure reallocate_coord_n
  module procedure reallocate_coord_lat
end interface

!+
! Function pointer_to_branch
!
! Routine to return a pointer to the lattice branch associated with a given name
! or a given element.
!
! This routine is an overloaded name for:
!   pointer_to_branch_given_ele (ele) result (branch_ptr)
!   pointer_to_branch_given_name (branch_name, lat) result (branch_ptr)
!
! The lattice branch *associated* with a given element is not necessarily the
! branch where the element is *located*. For example, all lords live in branch #0.
! But the branch associated with a super_lord element is the branch of its slaves.
!
! To get the branch where the element is located, simply use ele%ix_branch.
! 
! Note: Result is ambiguous if ele argument is associated with multiple branches 
! which can happen, for example, with overlay_lord elements.
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   ele         -- Ele_struct: Element contained in the branch.
!   branch_name -- Character(*): May be a branch name or a branch index.
!   lat         -- Lat_struct: Lattice to search.
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
!   pointer_to_field_ele
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
! Function physical_ele_end (stream_end, ele_orientation) result (physical_end)
!
! Rotine to determine which physical end of an element a particle is at given 
! the position in terms of upstream/downstream and the element's orientation
!
! Input:
!   stream_end       -- Integer: Either upstream_end$ or downstream_end$
!   ele_orientation  -- Integer: Either 1 = Normal or -1 = element reversed.
!
! Output:
!   physical_end     -- Integer: Either entrance_end$ or exit_end$
!-

function physical_ele_end (stream_end, ele_orientation) result (physical_end)

implicit none

integer stream_end, ele_orientation, physical_end

!

select case (ele_orientation)

case (1) 
  select case (stream_end)
  case (upstream_end$);   physical_end = entrance_end$
  case (downstream_end$); physical_end = exit_end$
  end select

case (-1) 
  select case (stream_end)
  case (upstream_end$);   physical_end = exit_end$
  case (downstream_end$); physical_end = entrance_end$
  end select

end select

end function physical_ele_end 

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

implicit none

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
  if (global_com%type_out) call out_io (s_fatal$, r_name, &
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
! Subroutine transfer_twiss (ele_in, ele_out)
!
! Routine to transfer the twiss parameters from one element to another.
!
! Moduels needed:
!   use bmad
!
! Input:
!   ele_in   -- Ele_struct: Element with existing Twiss parameters.
!
! Output:
!   ele_out  -- Ele_struct: Element receiving the Twiss parameters.
!-

subroutine transfer_twiss (ele_in, ele_out)

implicit none

type (ele_struct) ele_in, ele_out

!

ele_out%x       = ele_in%x
ele_out%y       = ele_in%y
ele_out%a       = ele_in%a
ele_out%b       = ele_in%b
ele_out%z       = ele_in%z
ele_out%c_mat   = ele_in%c_mat
ele_out%gamma_c = ele_in%gamma_c

end subroutine transfer_twiss

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

implicit none

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

implicit none

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
    if (c1%ix_slave == c2%ix_slave .and. c1%ix_branch == c2%ix_branch .and. &
                                         c1%ix_attrib == c2%ix_attrib) then
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
! Subroutine init_coord1 (orb, vec, ele, at_downstream_end, particle, E_photon, t_ref_offset, shift_vec6)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord1 (orb, vec, ele, at_downstream_end, particle, E_photon, t_ref_offset, shift_vec6)

implicit none

type (coord_struct) orb, orb2
type (coord_struct), save :: init_orb
type (ele_struct), optional, target :: ele

real(rp), optional :: vec(:), E_photon, t_ref_offset
real(rp) p0c, e_tot, ref_time

integer species
integer, optional :: particle
logical, optional :: at_downstream_end, shift_vec6

character(16), parameter :: r_name = 'init_coord1'

! Use temporary orb2 so if actual arg for vec, particle, or E_photon
! is part of the orb actual arg things do not get overwriten.

orb2 = init_orb                   ! See definition of coord_struct for default values.

orb2%state = alive$
orb2%p0c = 0

! Set %vec

if (present(vec)) then
  orb2%vec = vec
else
  orb2%vec = 0
endif

! Set %location

orb2%location = upstream_end$
if (present(at_downstream_end)) then
  if (at_downstream_end) orb2%location = downstream_end$
endif

! set species

if (present(particle)) then
  species = particle
elseif (present(ele)) then
  if (associated (ele%branch)) species = ele%branch%param%particle
elseif (init_orb%state == not_set$) then
  species = positron$
endif

! Energy values

if (present(ele)) then
  if (.not. present(at_downstream_end)) then
    call out_io (s_fatal$, r_name, 'Rule: "at_downstream_end" argument must be present if "ele" argument is.')
    call err_exit
  endif
  if (at_downstream_end .or. ele%key == init_ele$) then
    p0c = ele%value(p0c$)
    e_tot = ele%value(e_tot$)
    ref_time = ele%ref_time
    orb2%s = ele%s
  else
    p0c = ele%value(p0c_start$)
    e_tot = ele%value(e_tot_start$)
    ref_time = ele%value(ref_time_start$)
    orb2%s = ele%s - ele%value(l$)
  endif
endif

! Photon

if (species == photon$) then

  if (present(ele)) orb2%p0c = p0c

  if (present(E_photon)) then
    if (E_photon /= 0) orb2%p0c = E_photon
  endif

  if (orb2%vec(6) >= 0) orb2%vec(6) = sqrt(1 - orb2%vec(2)**2 - orb2%vec(4)**2)
  orb2%beta = 1

endif

! If ele is present...

orb2%ix_ele = -1

if (present(ele)) then

  orb2%ix_ele = ele%ix_ele
  if (ele%slave_status == slice_slave$) orb2%ix_ele = ele%lord%ix_ele

  if (ele%key == init_ele$) orb2%location = downstream_end$

  if (species /= photon$) then

    orb2%p0c = p0c

    ! Only time p0c_start /= p0c for an init_ele is when there is an e_gun present in the branch.
    if (ele%key == init_ele$ .and. logic_option(.true., shift_vec6)) then
      orb2%vec(6) = orb2%vec(6) + (ele%value(p0c_start$) - ele%value(p0c$)) / ele%value(p0c$)
    endif

    if (orb2%vec(6) == 0) then
      orb2%beta = p0c / e_tot
    else
      call convert_pc_to (p0c * (1 + orb2%vec(6)), species, beta = orb2%beta)
    endif

    ! Do not set %t if %beta = 0 since %t may be a good value.

    if (orb2%beta == 0) then
      if (orb2%vec(5) /= 0) then
        call out_io (s_error$, r_name, 'Z-POSITION IS NONZERO WITH BETA = 0.', &
                                       'THIS IS NONSENSE SO SETTING Z TO ZERO.')
        orb2%vec(5) = 0
      endif
    else
      orb2%t = ref_time - orb2%vec(5) / (orb2%beta * c_light)
      if (present(t_ref_offset)) orb2%t = orb2%t + t_ref_offset
    endif
  endif

endif

orb = orb2

end subroutine init_coord1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord2 (orb, orb_in, ele, at_downstream_end, t_ref_offset, shift_vec6)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord2 (orb, orb_in, ele, at_downstream_end, particle, t_ref_offset, shift_vec6)

implicit none

type (coord_struct) orb, orb_in, orb_save
type (ele_struct), optional :: ele
real(rp), optional :: t_ref_offset
integer, optional :: particle
logical, optional :: at_downstream_end, shift_vec6

!

orb_save = orb_in  ! Needed if actual args orb and orb_in are the same.

call init_coord1 (orb, orb_in%vec, ele, at_downstream_end, particle, orb_in%p0c, t_ref_offset, shift_vec6)

orb%spin      = orb_save%spin
orb%e_field_x = orb_save%e_field_x
orb%e_field_y = orb_save%e_field_y
orb%phase_x   = orb_save%phase_x
orb%phase_y   = orb_save%phase_y
orb%charge    = orb_save%charge
if (orb%beta == 0) orb%t = orb_save%t

end subroutine init_coord2

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

implicit none

type (coord_struct) orbit
integer particle
logical is_moving_forward

!

is_moving_forward = (orbit%state == alive$) .and. (orbit%p0c > 0)

end function particle_is_moving_forward

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function particle_is_moving_backward (orbit) result (is_moving_backward)
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
!   is_moving_backward -- Logical: True if moving backward. False otherwise.
!-

function particle_is_moving_backward (orbit) result (is_moving_backward)

implicit none

type (coord_struct) orbit
logical is_moving_backward

!

is_moving_backward = (orbit%state == alive$) .and. (orbit%p0c < 0)

end function particle_is_moving_backward

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)
!
! Function to convert a character string  (eg: "drift") to an index (eg: drift$).
!
! Modules needed:
!   use bmad
!
! Input:
!   key_str        -- Character(*): Name of the key. Result is case insensitive.
!   abbrev_allowed -- Logical, optional: Abbreviations (eg: "quad") allowed?
!                       Default is False. At least 3 characters are needed 
!                       (except for rfcavity elements) if True.
!
! Output:
!   key_index -- Integer: Index of the key. Set to -1 if key_name not recognized.
!-

function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)

implicit none

character(*) key_str
character(16) name

logical, optional :: abbrev_allowed
logical abbrev

integer key_index
integer i, n_name, n_match

!

n_match = 0
key_index = -1
if (key_str == '') return

call str_upcase (name, key_str)
call string_trim (name, name, n_name)

abbrev = logic_option(.false., abbrev_allowed)

do i = 1, n_key$
  if (abbrev .and. (n_name > 2 .or. name(1:2) == "RF")) then
    if (name(:n_name) == key_name(i)(1:n_name)) then
      key_index = i
      n_match = n_match + 1
    endif
  else
    if (name == key_name(i)) then
      key_index = i
      return
    endif
  endif
enddo

if (abbrev .and. n_match > 1) key_index = -1  ! Multiple matches are not valid

end function key_name_to_key_index 

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

implicit none

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

implicit none

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

end subroutine zero_ele_kicks

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_has_offset (ele) result (haz_offset)
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
!   haz_offset -- Logical: Set true is element has a non-zero offset.
!-

function ele_has_offset (ele) result (haz_offset)

implicit none

type (ele_struct) ele
logical haz_offset

!

haz_offset = .false.
if (.not. has_orientation_attributes(ele)) return

if (ele%value(tilt_tot$) /= 0) haz_offset = .true.
if (ele%value(x_pitch_tot$) /= 0) haz_offset = .true.
if (ele%value(y_pitch_tot$) /= 0) haz_offset = .true.
if (ele%value(x_offset_tot$) /= 0) haz_offset = .true.
if (ele%value(y_offset_tot$) /= 0) haz_offset = .true.
if (ele%value(z_offset_tot$) /= 0) haz_offset = .true.

end function ele_has_offset

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine zero_ele_offsets (ele)
!
! Subroutine to zero the offsets, pitches and tilt of an element.
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

implicit none

type (ele_struct) ele

!

if (.not. has_orientation_attributes(ele)) return

ele%value(tilt$) = 0
ele%value(x_pitch$) = 0
ele%value(y_pitch$) = 0
ele%value(x_offset$) = 0
ele%value(y_offset$) = 0
ele%value(z_offset$) = 0

ele%value(tilt_tot$) = 0
ele%value(x_pitch_tot$) = 0
ele%value(y_pitch_tot$) = 0
ele%value(x_offset_tot$) = 0
ele%value(y_offset_tot$) = 0
ele%value(z_offset_tot$) = 0

end subroutine zero_ele_offsets

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, dbeta, err_flag)
!
! Routine to calculate the momentum, etc. from a particle's total energy.
!
! Modules needed:
!   use bmad
!
! Input:
!   E_tot    -- Real(rp): Total energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma    -- Real(rp), optional: Gamma factor. Set to -1 for photons.
!   kinetic  -- Real(rp), optional: Kinetic energy
!   beta     -- Real(rp), optional: velocity / c_light
!   pc       -- Real(rp), optional: Particle momentum
!   brho     -- Real(rp), optional: Nominal B_field*rho_bend
!   dbeta    -- Real(rp), optional: 1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, dbeta, err_flag)

implicit none

real(rp), intent(in) :: E_tot
real(rp), intent(out), optional :: kinetic, beta, pc, brho, gamma, dbeta
real(rp) pc_new, mc2, g2

integer, intent(in) :: particle
logical, optional :: err_flag

character(24) :: r_name = 'convert_total_energy_to'

!

if (present(err_flag)) err_flag = .true.

mc2 = mass_of(particle)
if (E_tot < mc2) then
  call out_io (s_abort$, r_name, 'ERROR: TOTAL ENERGY IS LESS THAN REST MASS:\f10.0\ ', E_tot)
  if (global_com%exit_on_error) call err_exit
  return
endif

pc_new = E_tot * sqrt(1.0 - (mc2/E_tot)**2)
if (present(pc))     pc     = pc_new
if (present(beta))    beta    = pc_new / E_tot  
if (present(kinetic)) kinetic = E_tot - mc2
if (present(brho))    brho    = pc_new / c_light

if (present(gamma)) then
  if (mc2 == 0) then
    gamma = -1
  else
    gamma   = E_tot / mc2
  endif
endif

if (present(dbeta)) then
  if (E_tot/mc2 > 100) then
    g2 = (E_tot / mc2)**2
    dbeta = 1/(2*g2) + 1/(8*g2**2)
  else
    dbeta = 1 - pc_new / E_tot
  endif
endif

if (present(err_flag)) err_flag = .false.

end subroutine convert_total_energy_to

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, dbeta, err_flag)
!
! Routine to calculate the energy, etc. from a particle's momentum.
!
! Modules needed:
!   use bmad
!
! Input:
!   pc       -- Real(rp): Particle momentum
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   E_tot    -- Real(rp), optional: Total energy of the particle.
!   gamma    -- Real(rp), optional: Gamma factor.
!   kinetic  -- Real(rp), optional: Kinetic energy
!   beta     -- Real(rp), optional: velocity / c_light
!   brho     -- Real(rp), optional: Nominal B_field*rho_bend
!   dbeta    -- Real(rp), optional: 1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, dbeta, err_flag)

implicit none

real(rp), intent(in) :: pc
real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma, dbeta
real(rp) g2, mc2, E_tot_this 

integer, intent(in) :: particle
logical, optional :: err_flag

character(20) :: r_name = 'convert_pc_to'

!

if (present(err_flag)) err_flag = .false.

mc2 = mass_of(particle)
E_tot_this = sqrt(pc**2 + mc2**2)

if (present(E_tot))   E_tot   = E_tot_this
if (present(beta))    beta    = pc / E_tot_this
if (present(kinetic)) kinetic = E_tot_this - mc2
if (present(brho))    brho    = pc / c_light
if (present(gamma))   gamma   = E_tot_this / mc2

if (present(dbeta)) then
  if (E_tot/mc2 > 100) then
    g2 = (E_tot_this / mc2)**2
    dbeta = 1/(2*g2) + 1/(8*g2**2)
  else
    dbeta = 1 - pc / E_tot_this
  endif
endif

if (present(err_flag)) err_flag = .false.

end subroutine convert_pc_to

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine transfer_lat_parameters (lat_in, lat_out)
!
! Subroutine to transfer the lat parameters (such as lat%name, lat%param, etc.)
! from one lat to another. The only stuff that is not transfered are things
! that are (or have) pointers or arrays
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat with parameters set.
!-

subroutine transfer_lat_parameters (lat_in, lat_out)

implicit none

type (lat_struct), intent(in) :: lat_in
type (lat_struct) :: lat_out

!

lat_out%use_name                  = lat_in%use_name
lat_out%lattice                   = lat_in%lattice
lat_out%input_file_name           = lat_in%input_file_name
lat_out%title                     = lat_in%title
lat_out%a                         = lat_in%a
lat_out%b                         = lat_in%b
lat_out%z                         = lat_in%z
lat_out%param                     = lat_in%param
lat_out%lord_state                = lat_in%lord_state
lat_out%beam_start                = lat_in%beam_start
lat_out%pre_tracker               = lat_in%pre_tracker
lat_out%version                   = lat_in%version
lat_out%n_ele_track               = lat_in%n_ele_track
lat_out%n_ele_max                 = lat_in%n_ele_max
lat_out%n_control_max             = lat_in%n_control_max
lat_out%n_ic_max                  = lat_in%n_ic_max
lat_out%input_taylor_order        = lat_in%input_taylor_order
lat_out%absolute_time_tracking    = lat_in%absolute_time_tracking
lat_out%rf_auto_scale_phase       = lat_in%rf_auto_scale_phase
lat_out%rf_auto_scale_amp         = lat_in%rf_auto_scale_amp
lat_out%use_ptc_layout            = lat_in%use_ptc_layout

end subroutine transfer_lat_parameters

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)
!
! Subroutine to transfer a Taylor map from one element to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele_in       -- Ele_struct: Element with the Taylor map.
!   taylor_order -- Integer, optional: Order to truncate the Taylor map at.
!
! Output:
!   ele_out      -- Ele_struct: Element receiving the Taylor map truncated to
!                     order taylor_order.
!-

subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)

implicit none

type (ele_struct) ele_in, ele_out
integer, optional :: taylor_order
integer it, ix, k 

!

do it = 1, 6

  if (present(taylor_order)) then
    ix = 0
    do k = 1, size(ele_in%taylor(it)%term)
      if (sum(ele_in%taylor(it)%term(k)%expn(:)) > taylor_order) cycle
      ix = ix + 1
    enddo
  else
    ix = size(ele_in%taylor(it)%term)
  endif

  if (.not. associated(ele_out%taylor(it)%term)) allocate (ele_out%taylor(it)%term(ix))
  if (size(ele_out%taylor(it)%term) /= ix) allocate (ele_out%taylor(it)%term(ix))

  ix = 0
  do k = 1, size(ele_in%taylor(it)%term)
    if (present(taylor_order)) then
      if (sum(ele_in%taylor(it)%term(k)%expn(:)) > taylor_order) cycle
    endif
    ix = ix + 1
    ele_out%taylor(it)%term(ix) = ele_in%taylor(it)%term(k)
  enddo

enddo

ele_out%taylor(:)%ref = ele_in%taylor(:)%ref

end subroutine transfer_ele_taylor

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_lat (lat, n)
!
! Subroutine to initialize a BMAD lat.
! 
! Modules needed:
!   use bmad
!
! Input:
!   n    -- Integer, optional: Upper bound lat%ele(0:) array is initialized to.
!
! Output:
!   lat -- lat_struct: Initialized lat.
!-

subroutine init_lat (lat, n)

implicit none

type (lat_struct)  lat
type (lat_param_struct), save :: param0

integer, optional :: n

!

call init_attribute_name_array
call deallocate_lat_pointers (lat)
if (present(n)) call allocate_lat_ele_array(lat, n)
call init_ele (lat%ele_init)

call reallocate_control (lat, 100)

lat%title = ' '
lat%use_name = ' '
lat%lattice = ' '
lat%input_file_name = ' '

lat%param = param0
call set_status_flags (lat%param%bookkeeping_state, stale$)

call init_mode_info (lat%a)
call init_mode_info (lat%b)
call init_mode_info (lat%z)

lat%n_ele_track = 0
lat%n_ele_max = 0
lat%n_control_max = 0
lat%n_ic_max = 0
lat%input_taylor_order = 0
lat%version = -1
lat%absolute_time_tracking   = bmad_com%absolute_time_tracking_default
lat%rf_auto_scale_phase      = bmad_com%rf_auto_scale_phase_default
lat%rf_auto_scale_amp        = bmad_com%rf_auto_scale_amp_default
lat%use_ptc_layout           = bmad_com%use_ptc_layout_default

call allocate_branch_array (lat, 0)

!----------------------------------------
contains

subroutine init_mode_info (t)
type (mode_info_struct) t
t%tune = 0
t%emit = 0
t%chrom = 0
end subroutine init_mode_info

end subroutine init_lat

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

implicit none

type (ele_struct) :: ele_taylor, ele2

integer it

logical equiv
logical vmask(num_ele_attrib$), vnot(num_ele_attrib$)

!

equiv = .false.

if (ele_taylor%key == taylor$) return  
if (ele_taylor%key /= ele2%key) return
if (ele_taylor%sub_key /= ele2%sub_key) return
if (ele_taylor%map_with_offsets .neqv. ele2%map_with_offsets) return
if (ele_taylor%value(integrator_order$) /= ele2%value(integrator_order$)) return

vmask = .true.
vmask(delta_ref_time$) = .false.
vmask(ref_time_start$) = .false.
if (ele_taylor%key == wiggler$ .and. ele_taylor%sub_key == map_type$) then
  vmask( [k1$, rho$, b_max$] ) = .false.  ! These are dependent attributes.
endif
if (.not. ele_taylor%map_with_offsets) then
  vmask( [x_offset$, y_offset$, z_offset$, tilt$, x_pitch$, &
            y_pitch$, x_offset_tot$, y_offset_tot$, z_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$] ) = .false.
endif

vnot = (ele_taylor%value /= ele2%value)
vnot = vnot .and. vmask
if (any(vnot)) return

if (associated(ele_taylor%wig) .neqv. associated(ele2%wig)) return
if (associated(ele_taylor%wig)) then
  if (size(ele_taylor%wig%term) /= size(ele2%wig%term)) return
  do it = 1, size(ele_taylor%wig%term)
    if (ele_taylor%wig%term(it)%coef  /= ele2%wig%term(it)%coef)  cycle
    if (ele_taylor%wig%term(it)%kx    /= ele2%wig%term(it)%kx)    cycle
    if (ele_taylor%wig%term(it)%ky    /= ele2%wig%term(it)%ky)    cycle
    if (ele_taylor%wig%term(it)%kz    /= ele2%wig%term(it)%kz)    cycle
    if (ele_taylor%wig%term(it)%phi_z /= ele2%wig%term(it)%phi_z) cycle
  enddo
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

implicit none

type (lat_struct) lat

lat%param%t1_no_RF = 0
lat%param%t1_with_RF = 0

end subroutine clear_lat_1turn_mats


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_ele (ele1, ele2, nullify_pointers)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal operator.
! The result is that ele2's pointers will point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1             -- Ele_struct:
!   nullify_pointers -- Logical, optional: If present and True then nullify the 
!                         pointers in ele2 except for the ele2%lat and ele2%lord pointers. 
!                         This gives a "bare bones" copy where one does not have to 
!                         worry about deallocating allocated structure components later.
!
! Output:
!   ele2 -- Ele_struct:
!-

subroutine transfer_ele (ele1, ele2, nullify_pointers)

type (ele_struct), target :: ele1
type (ele_struct) :: ele2
logical, optional :: nullify_pointers

!

ele2 = ele1

if (logic_option (.false., nullify_pointers)) then
  call deallocate_ele_pointers (ele2, .true.)
  ele2%branch => ele1%branch  ! Reinstate
  ele2%lord   => ele1%lord    ! Reinstate
endif

end subroutine transfer_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_eles (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1(:) -- Ele_struct:
!
! Output:
!   ele2(:) -- Ele_struct:
!-

subroutine transfer_eles (ele1, ele2)

type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(inout) :: ele2(:)

ele2 = ele1

end subroutine transfer_eles

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branch (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   branch1 -- Branch_struct:
!
! Output:
!   branch2 -- Branch_struct:
!-

subroutine transfer_branch (branch1, branch2)

type (branch_struct) :: branch1
type (branch_struct) :: branch2

!

branch2 = branch1

end subroutine transfer_branch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branches (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   branch1(:) -- Branch_struct:
!
! Output:
!   branch2(:) -- Branch_struct:
!-

subroutine transfer_branches (branch1, branch2)

type (branch_struct) :: branch1(:)
type (branch_struct) :: branch2(:)

branch2 = branch1

end subroutine transfer_branches

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_lat (lat1, lat2)
!
! Subroutine to set lat2 = lat1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end lat2's pointers point to the same memory as lat1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   lat1 -- lat_struct:
!
! Output:
!   lat2 -- lat_struct:
!-

subroutine transfer_lat (lat1, lat2)

type (lat_struct), intent(in) :: lat1
type (lat_struct), intent(out) :: lat2

lat2 = lat1

end subroutine transfer_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_n (coord, n_coord)
!
! Subroutine to allocate an allocatable  coord_struct array.
! This is an overloaded subroutine. See reallocate_coord.
!-

subroutine reallocate_coord_n (coord, n_coord)

type (coord_struct), allocatable :: coord(:)
type (coord_struct), allocatable :: old(:)

integer, intent(in) :: n_coord
integer i, n_old

!

if (allocated (coord)) then

  n_old = ubound(coord, 1)
  if (n_old >= n_coord) return
  allocate(old(0:n_old))

  do i = 0, n_old
    old(i) = coord(i)
  enddo

  deallocate (coord)
  allocate (coord(0:n_coord))

  do i = 0, n_old
    coord(i) = old(i)
  enddo

  deallocate(old)

else
  allocate (coord(0:n_coord))
  call init_coord (coord(0))
endif

end subroutine reallocate_coord_n

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_lat (coord, lat, ix_branch)
!
! Subroutine to allocate an allocatable  coord_struct array.
! This is an overloaded subroutine. See reallocate_coord.
!-

subroutine reallocate_coord_lat (coord, lat, ix_branch)

type (coord_struct), allocatable :: coord(:)
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer, optional :: ix_branch

!

branch => lat%branch(integer_option(0, ix_branch))

if (allocated(coord)) then
  call reallocate_coord_n (coord, branch%n_ele_max)
else
  allocate (coord(0:branch%n_ele_max))
  call init_coord (coord(0), ele = branch%ele(0), at_downstream_end = .true.)
endif

end subroutine reallocate_coord_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord_array (coord_array, lat)
!
! Subroutine to allocate an allocatable coord_array_struct array to
! the proper size for a lattice.
!
! Note: Any old coordinates are not saved except for coord_array(:)%orb(0).
! If, at input, coord_array is not allocated, coord_array(:)%orb(0)%vec is set to zero.
! In any case, all other %vec components are set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord(:) -- Coord_struct, allocatable: Allocatable array.
!   lat      -- lat_struct: 
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord_array (coord_array, lat)

implicit none

type (coord_array_struct), allocatable :: coord_array(:)
type (lat_struct) lat
type (coord_struct), allocatable :: start(:)

integer i, j, nb

!

if (.not. allocated(lat%branch)) return
nb = ubound(lat%branch, 1)

if (allocated (coord_array)) then
  if (size(coord_array) /= nb + 1) then
    call reallocate_coord(start, nb)
    do i = 0, nb
      start(i) = coord_array(i)%orb(0)
    enddo
    deallocate (coord_array)
    allocate (coord_array(0:nb))
    do i = 0, nb
      call reallocate_coord (coord_array(i)%orb, lat%branch(i)%n_ele_max)
      coord_array(i)%orb(0) = start(i)
    enddo
  endif
else
  allocate (coord_array(0:nb))
  do i = 0, nb
    call reallocate_coord (coord_array(i)%orb, lat%branch(i)%n_ele_max)
  enddo
endif

end subroutine reallocate_coord_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_control(lat, n) 
!
! Function to reallocate the lat%control(:) and lat%ic(:) arrays.
! The old data in the arrays will be saved.
! 
! Modules needed:
!   use bmad
!
! Input:
!   lat  -- Lat_struct: Lattice.
!   n    -- Integer: Array size for lat%control(:) and lat%ic(:).
!
! Output:
!   lat  -- Lat_struct: Lattice.
!     %control(:) -- Control Array with size at least n.
!     %ic(:)      -- Control Array.
!-

subroutine reallocate_control (lat, n)

implicit none

type (lat_struct) lat
type (control_struct), allocatable :: control(:)
integer, intent(in) :: n
integer n_old

!

if (.not. allocated(lat%control)) then
  allocate (lat%control(n), lat%ic(n))
  return
endif

n_old = size(lat%control)
if (n_old >= n) return

allocate (control(n_old))
control = lat%control

deallocate (lat%control)
allocate (lat%control(n))
lat%control(1:n_old) = control
deallocate (control)

call re_allocate(lat%ic, max(n, size(lat%ic) + n - n_old))
lat%ic(n_old+1:) = 0

end subroutine reallocate_control

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine deallocate_ele_pointers (ele, nullify_only, nullify_branch, dealloc_poles)
!
! Subroutine to deallocate the pointers in an element.
! Note: ele%branch is always nullified. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele            -- ele_struct: Element with pointers.
!   nullify_only   -- Logical, optional: If present and True: Nullify & do not deallocate.
!   nullify_branch -- Logical, optional: Nullify ele%branch? Default is True.
!   dealloc_poles  -- Logical, optional: Dealloc ele%a_pole, ele%b_pole? Default is True.
!
! Output:
!   ele -- Ele_struct: Element with deallocated pointers.
!-

subroutine deallocate_ele_pointers (ele, nullify_only, nullify_branch, dealloc_poles)

implicit none

type (ele_struct), target :: ele
type (em_field_mode_struct), pointer :: mode
logical, optional, intent(in) :: nullify_only, nullify_branch, dealloc_poles
integer i

! %lord and %lat never point to something that has been allocated for the element
! so just nullify these pointers.

if (logic_option(.true., nullify_branch)) nullify (ele%branch)
nullify (ele%lord)

! nullify

if (logic_option (.false., nullify_only)) then
  nullify (ele%wig)
  nullify (ele%rad_int_cache)
  nullify (ele%r)
  nullify (ele%descrip)
  nullify (ele%a_pole, ele%b_pole)
  nullify (ele%rf_wake)
  nullify (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
            ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)
  nullify (ele%ptc_genfield)
  nullify (ele%ptc_fibre)
  nullify (ele%mode3)
  nullify (ele%wall3d)
  nullify (ele%em_field)
  return
endif

! Normal deallocate.

if (associated (ele%a_pole) .and. logic_option(.true., dealloc_poles)) &
                                     deallocate (ele%a_pole, ele%b_pole)
if (associated (ele%rad_int_cache))  deallocate (ele%rad_int_cache)
if (associated (ele%r))              deallocate (ele%r)
if (associated (ele%descrip))        deallocate (ele%descrip)
if (associated (ele%mode3))          deallocate (ele%mode3)

if (associated (ele%rf_wake)) then
  if (associated (ele%rf_wake%sr_table))      deallocate (ele%rf_wake%sr_table)
  if (associated (ele%rf_wake%sr_mode_long))  deallocate (ele%rf_wake%sr_mode_long)
  if (associated (ele%rf_wake%sr_mode_trans)) deallocate (ele%rf_wake%sr_mode_trans)
  if (associated (ele%rf_wake%lr))            deallocate (ele%rf_wake%lr)
  deallocate (ele%rf_wake)
endif

call deallocate_wall3d_pointer (ele%wall3d)

if (associated (ele%em_field)) then
  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    if (associated (mode%map)) then
      mode%map%n_link = mode%map%n_link - 1
      if (mode%map%n_link == 0) deallocate (ele%em_field%mode(i)%map)
    endif
    if (associated (mode%grid)) then
      mode%grid%n_link = mode%grid%n_link - 1
      if (mode%grid%n_link == 0) deallocate (ele%em_field%mode(i)%grid)
    endif
  enddo
  deallocate (ele%em_field)
endif

if (associated(ele%wig)) then
  ele%wig%n_link = ele%wig%n_link - 1
  if (ele%wig%n_link == 0) then
    deallocate (ele%wig)
  else
    nullify (ele%wig)
  endif
endif

if (associated (ele%taylor(1)%term)) deallocate &
         (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
         ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)

call kill_ptc_genfield (ele%ptc_genfield)

end subroutine deallocate_ele_pointers

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine kill_ptc_genfield (ptc_genfield)
!
! Subroutine to kill a ptc_genfield.
!
! Modules needed:
!   use bmad
!
! Input:
!   ptc_genfield -- Genfield, pointer: ptc_genfield to kill.
!
! Output:
!   ptc_genfield -- Genfield, pointer: Killed ptc_genfield.
!-

subroutine kill_ptc_genfield (ptc_genfield)

use tpsalie_analysis, only: kill 

implicit none

type (genfield), pointer :: ptc_genfield

!

if (associated(ptc_genfield)) then
  call kill (ptc_genfield)
  deallocate (ptc_genfield)
endif

end subroutine kill_ptc_genfield

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine init_ele (ele, key, sub_key, ix_ele, ix_branch, branch)
!
! Subroutine to initialize a Bmad element. Element is initialized to be free
! (not a lord or slave) and all %values set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   key       -- Integer, optional: Key to initialize to. EG: quadrupole$, etc.
!   sub_key   -- Integer, optional: Sub-key to initialize to.
!   ix_ele    -- Integer, optional: ix_ele index to initalize to. Default = -1.
!   ix_branch -- Integer, optional: Branch index to initalize to. Default = 0.
!   branch    -- branch_struct: Branch to point ele%branch to. Otherwise ele%branch is nullified.
!
! Output:
!   ele -- Ele_struct: Initialized element.
!-

subroutine init_ele (ele, key, sub_key, ix_ele, ix_branch, branch)

implicit none

type (ele_struct)  ele
type (branch_struct), optional, target :: branch
integer, optional :: key, sub_key
integer, optional :: ix_branch, ix_ele

!

call deallocate_ele_pointers (ele)
if (present(branch)) ele%branch => branch

ele%type = ' '
ele%alias = ' '
ele%name = '<Initialized>'
ele%component_name = ' '

ele%key = integer_option (0, key)
ele%sub_key = integer_option (0, sub_key)
if (present(key)) call set_ele_defaults(ele)

ele%value(:) = 0
ele%old_value(:) = 0
ele%map_ref_orb_in = 0
ele%map_ref_orb_out = 0
ele%time_ref_orb_in = 0
ele%time_ref_orb_out = 0

ele%lord_status = not_a_lord$
ele%slave_status = free$
ele%ix_value = 0
ele%ic1_lord = 0
ele%ic2_lord = -1
ele%n_lord = 0
ele%ix1_slave = 0
ele%ix2_slave = -1
ele%n_slave = 0
ele%ix_pointer = 0
ele%s = 0
ele%ref_time = 0
ele%ix_branch = 0
ele%ix_ele = -1
ele%orientation       = 1

ele%ixx = 0
ele%iyy = 0

call set_status_flags (ele%bookkeeping_state, stale$)

if (present(ix_branch)) ele%ix_branch = ix_branch
if (present(ix_ele)) ele%ix_ele = ix_ele

call init_floor (ele%floor)

ele%mat6_calc_method     = bmad_standard$
ele%tracking_method      = bmad_standard$
ele%spin_tracking_method = bmad_standard$
ele%field_calc           = bmad_standard$
ele%ptc_integration_type = matrix_kick$

ele%is_on             = .true.
ele%multipoles_on     = .true.
ele%scale_multipoles  = .true.
ele%symplectify       = .false.
ele%map_with_offsets  = .true.
ele%on_a_girder       = .false.
ele%csr_calc_on       = .true.
ele%logic             = .false.
ele%mode_flip         = .false.
ele%field_master      = .false.
ele%offset_moves_aperture = .false.

ele%aperture_type = rectangular$
ele%aperture_at   = downstream_end$

! init Twiss

ele%c_mat = 0
ele%gamma_c = 1.0

ele%x%eta  = 0
ele%x%etap = 0

ele%y%eta  = 0
ele%y%etap = 0

ele%a%beta     = 0
ele%a%alpha    = 0
ele%a%gamma    = 0
ele%a%eta      = 0
ele%a%etap     = 0
ele%a%phi      = 0
ele%a%sigma    = 0
ele%a%emit     = 0

ele%b%beta     = 0
ele%b%alpha    = 0
ele%b%gamma    = 0
ele%b%eta      = 0
ele%b%etap     = 0
ele%b%phi      = 0
ele%b%sigma    = 0
ele%b%emit     = 0

ele%z%beta     = 0
ele%z%alpha    = 0
ele%z%gamma    = 0
ele%z%eta      = 0
ele%z%etap     = 1
ele%z%phi      = 0
ele%z%sigma    = 0
ele%z%emit     = 0

end subroutine init_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+ 
! Subroutine init_floor (floor)
!
! Routine to initialize a floor_position_struct to zero.
!
! Output:
!   floor -- Floor_position_struct: Floor coordinates to init.
!-

subroutine init_floor (floor)

implicit none

type (floor_position_struct) floor

!

floor%r     = 0
floor%theta = 0
floor%phi   = 0
floor%psi   = 0

end subroutine init_floor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:)  -- Element array to reallocate.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 100 if ele is not allocated.
!   ix_branch   -- Integer, optional: Branch index. Default is 0.
!
! Output:
!   lat         -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:) -- Ele_struct, pointer: Resized element array.
!-

subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch)

implicit none

type (lat_struct), target :: lat
integer, optional :: upper_bound
integer, optional :: ix_branch
integer ix_br, i

!

ix_br = integer_option (0, ix_branch)

if (ix_br == 0) then
  call allocate_element_array (lat%ele, upper_bound, .true.)
  if (allocated(lat%branch)) then
    do i = 0, ubound(lat%ele, 1)
      lat%ele(i)%branch => lat%branch(0)
    enddo
    lat%branch(0)%ele => lat%ele
  endif

else
  call allocate_element_array (lat%branch(ix_br)%ele, upper_bound, .true.)
  do i = 0, ubound(lat%branch(ix_br)%ele, 1)
    lat%branch(ix_br)%ele(i)%branch => lat%branch(ix_br)
  enddo
  lat%branch(ix_br)%ele%ix_branch = ix_br
endif


end subroutine allocate_lat_ele_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_element_array (ele, upper_bound, init_ele0)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Note: Use allocate_lat_ele_array instead for all ele(:) arrays that
!       are part of a lattice.
!   
!
! Modules needed:
!   use bmad
!
! Input:
!   ele(:)      -- Ele_struct, pointer: Element array.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 100 if ele is not allocated.
!   init_ele0   -- Logical, optional: If present and True and ele(:) array has not been allocated then set:
!                     ele(0)%name = 'BEGINNING'
!                     ele(0)%key = init_ele$
!                     ele(0)%mat6 = unit matrix
!
! Output:
!   ele(:)      -- Ele_struct, pointer: Allocated element array.
!-

subroutine allocate_element_array (ele, upper_bound, init_ele0)

implicit none

type (ele_struct), pointer :: ele(:)
type (ele_struct), pointer :: temp_ele(:)

integer, optional :: upper_bound
integer curr_ub, ub, i

logical, optional :: init_ele0

! get new size

ub = 10
if (associated (ele)) ub = max (int(1.3*size(ele)), ub)
if (present(upper_bound))  ub = upper_bound

!  save ele if present

if (associated (ele)) then
  if (ub == ubound(ele, 1)) return
  curr_ub = min(ub, ubound(ele, 1))
  allocate (temp_ele(0:curr_ub))
  call transfer_eles (ele(0:curr_ub), temp_ele)
  do i = curr_ub+1, ubound(ele, 1)
    call deallocate_ele_pointers(ele(i))
  enddo
  deallocate (ele)
  allocate(ele(0:ub))
  call transfer_eles (temp_ele(0:curr_ub), ele(0:curr_ub))
  deallocate (temp_ele)
else
  curr_ub = -1
  allocate(ele(0:ub))
endif

! 

do i = curr_ub+1, ub
  call init_ele (ele(i))
  ele(i)%ix_ele = i
end do

if (logic_option(.false., init_ele0) .and. curr_ub == -1) then
  ele(0)%name = 'BEGINNING'
  ele(0)%key = init_ele$
  call mat_make_unit (ele(0)%mat6)
  call set_ele_defaults(ele(0))
endif

end subroutine allocate_element_array

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_branch_array (lat, upper_bound, lat)
!
! Subroutine to allocate or re-allocate an branch array.
! The old information is saved.
! The lower bound is always 0.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- Lat_struct: 
!     %branch(:)  -- Branch array to be allocated.
!   upper_bound -- Integer: Desired upper bound.
! 
! Output:
!   lat         -- Lat_struct: 
!     %branch(:)  -- Allocated branch array.
!-

subroutine allocate_branch_array (lat, upper_bound)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (branch_struct), pointer :: temp_branch(:)

integer :: upper_bound
integer curr_ub, ub, i, j

character(20) :: r_name = 'allocate_branch_array'

!  save branch if present

ub = upper_bound
if (allocated (lat%branch)) then
  if (ub == ubound(lat%branch, 1)) return
  curr_ub = min(ub, ubound(lat%branch, 1))
  allocate (temp_branch(0:curr_ub))
  call transfer_branches (lat%branch(0:curr_ub), temp_branch)
  do i = curr_ub+1, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers(lat%branch(i)%ele)
    deallocate(lat%branch(i)%n_ele_track)
    deallocate(lat%branch(i)%n_ele_max)
  enddo
  deallocate (lat%branch)
  allocate(lat%branch(0:ub))
  call transfer_branches (temp_branch(0:curr_ub), lat%branch(0:curr_ub))
  deallocate (temp_branch)
else
  curr_ub = -1
  allocate(lat%branch(0:ub))
  lat%branch(0)%ele            => lat%ele
  lat%branch(0)%param          => lat%param
  lat%branch(0)%a              => lat%a
  lat%branch(0)%b              => lat%b
  lat%branch(0)%z              => lat%z
  lat%branch(0)%n_ele_track    => lat%n_ele_track
  lat%branch(0)%n_ele_max      => lat%n_ele_max
  if (associated(lat%ele)) then
    do i = 0, ubound(lat%ele, 1)
      lat%ele(i)%branch => lat%branch(0)
    enddo
  endif
endif

! 

do i = curr_ub+1, ub
  branch => lat%branch(i)
  branch%lat => lat
  branch%name = ''
  branch%ix_branch = i
  branch%ix_from_branch = -1
  branch%ix_from_ele = -1
  if (i == 0) cycle
  allocate(branch%n_ele_track)
  allocate(branch%n_ele_max)
  allocate(branch%param)
  allocate(branch%a, branch%b, branch%z)
  branch%param = lat%param
  call set_status_flags (branch%param%bookkeeping_state, stale$)
end do

do i = 0, ub
  branch => lat%branch(i)
  if (.not. associated (branch%ele)) cycle
  do j = 0, ubound(branch%ele, 1)
    branch%ele(j)%branch => branch
  enddo
enddo

end subroutine allocate_branch_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine deallocate_lat_pointers (lat)
!
! Subroutine to deallocate the pointers in a lat.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Lat with pointers.
!
! Output:
!   lat -- lat_struct: Lat with deallocated pointers.
!-

subroutine deallocate_lat_pointers (lat)

implicit none

type (lat_struct) lat
integer i

!

if (associated (lat%ele)) then
  call deallocate_ele_array_pointers (lat%ele)
  call deallocate_ele_pointers (lat%ele_init)
endif

if (allocated(lat%control))  deallocate (lat%control)
if (allocated(lat%ic))       deallocate (lat%ic)

! Do not need to deallocate stuff in lat%branch(0) since
! these pointers have been deallocated above.

if (allocated (lat%branch)) then
  do i = 1, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers (lat%branch(i)%ele)
    deallocate (lat%branch(i)%param, lat%branch(i)%a, lat%branch(i)%b, lat%branch(i)%z)
    call deallocate_wall3d_pointer (lat%branch(i)%wall3d)
  enddo
  deallocate (lat%branch)
endif

!

lat%n_ele_track  = -1
lat%n_ele_max  = -1

end subroutine deallocate_lat_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_ele_array_pointers (eles)
!
! Routine to deallocate the pointers of all the elements in an 
! element array and the array itself.
!
! Modules needed:
!   use bmad
!
! Input:
!   eles(:) -- Ele_struct, pointer: Array of elements.
!
! Output:
!   eles(:) -- Ele_struct, pointer: Deallocated array.
!-

subroutine deallocate_ele_array_pointers (eles)

implicit none

type (ele_struct), pointer :: eles(:)
integer i

!

do i = lbound(eles, 1), ubound(eles, 1)
  call deallocate_ele_pointers (eles(i))
enddo

deallocate (eles)

end subroutine deallocate_ele_array_pointers

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_wall3d_pointer (wall3d)
!
! Routine to deallocate a wall3d pointer.
!
! Input:
!   wall3d -- wall3d_struct, pointer: Pointer to wall3d structure.
!
! Output:
!   wall3d -- wall3d_struct, pointer: deallocated
!-

subroutine deallocate_wall3d_pointer (wall3d)

implicit none

type (wall3d_struct), pointer :: wall3d

!

if (associated (wall3d)) then
  wall3d%n_link = wall3d%n_link - 1
  if (wall3d%n_link == 0) then
    deallocate (wall3d%section)
    deallocate (wall3d)
  else
    nullify(wall3d)
  endif
endif

end subroutine deallocate_wall3d_pointer

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine transfer_mat_from_twiss (ele1, ele2, m)
!
! Subroutine to make a 6 x 6 transfer matrix from the twiss parameters
! at two points.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele1 -- Ele_struct: Element with twiss parameters for the starting point.
!     %a, %b -- a-mode and b-mode Twiss paramters
!       %beta   -- Beta parameter.
!       %alpha  -- Alpha parameter.
!       %phi    -- Phase at initial point.
!     %x  %y -- dispersion values
!       %eta    -- Dispersion at initial point.
!       %etap   -- Dispersion derivative at initial point.
!     %c_mat(2,2) -- Coupling matrix
!   ele2 -- Ele_struct: Element with twiss parameters for the ending point.
!
! Output:
!   m(6,6) -- Real(rp): Transfer matrix between the two points.
!-

subroutine transfer_mat_from_twiss (ele1, ele2, m)

implicit none

type (ele_struct) ele1, ele2

real(rp) m(6,6), v_mat(4,4), v_inv_mat(4,4), det
character(20) :: r_name = 'transfer_mat_from_twiss'

! Error check

if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then
  call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT: ' // ele1%name)
  if (global_com%exit_on_error) call err_exit
endif

if (ele2%a%beta == 0 .or. ele2%b%beta == 0) then
  call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT: ' // ele2%name)
  if (global_com%exit_on_error) call err_exit
endif

! Transfer matrices without coupling or dispersion

call mat_make_unit (m)
call transfer_mat2_from_twiss (ele1%a, ele2%a, m(1:2,1:2))
call transfer_mat2_from_twiss (ele1%b, ele2%b, m(3:4,3:4))

! Add in coupling

if (any(ele1%c_mat /= 0)) then
  det = determinant (ele1%c_mat)
  ele1%gamma_c = sqrt(1-det)
  call make_v_mats (ele1, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (m(1:4,1:4), v_inv_mat)
endif

if (any(ele2%c_mat /= 0)) then
  det = determinant (ele2%c_mat)
  ele2%gamma_c = sqrt(1-det)
  call make_v_mats (ele2, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (v_mat, m(1:4,1:4))
endif

! Add in dispersion.

m(1:4,6) = [ele2%x%eta, ele2%x%etap, ele2%y%eta, ele2%y%etap] - &
        matmul (m(1:4,1:4), [ele1%x%eta, ele1%x%etap, ele1%y%eta, ele1%y%etap]) 

! The m(5,x) terms follow from the symplectic condition.

m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)


end subroutine transfer_mat_from_twiss

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine match_ele_to_mat6 (ele, vec0, mat6, err_flag)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele -- Ele_struct: Match element.
!     %value(beta_a0$) -- Beta_a at the start
!
! Output:
!   vec0(6)   -- Real(rp): Currently just set to zero.
!   mat6(6,6) -- Real(rp): Transfer matrix.
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!-

subroutine match_ele_to_mat6 (ele, vec0, mat6, err_flag)

implicit none

type (ele_struct), target :: ele, ele0, ele1

real(rp) mat6(6,6), vec0(6)
real(rp), pointer :: v(:)

logical err_flag

! Special case where match_end is set but there is no beginning beta value yet.
! In this case, just return the unit matrix and set the err_flag.

if (ele%value(match_end$) /= 0 .and. (ele%value(beta_a0$) == 0 .or. ele%value(beta_b0$) == 0)) then
  call mat_make_unit (mat6)
  vec0 = 0
  err_flag = .true.
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

call transfer_mat_from_twiss (ele0, ele1, mat6)

! Kick part

vec0 = [v(x1$), v(px1$), v(y1$), v(py1$), v(z1$), v(pz1$)] - &
       matmul (mat6, [v(x0$), v(px0$), v(y0$), v(py0$), v(z0$), v(pz0$)])

end subroutine match_ele_to_mat6

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wig (wig_in, wig_out)
!
! Subroutine to point wig_out => wig_in
!
! Modules needed:
!   use bmad
!
! Input:
!   wig_in  -- Wig_struct, pointer: Input wiggler field.
!
! Output:
!   wig_out -- Wig_struct, pointer: Output wiggler field.
!-

subroutine transfer_wig (wig_in, wig_out)

implicit none

type (wig_struct), pointer :: wig_in, wig_out

!

if (.not. associated(wig_in) .and. .not. associated(wig_out)) return
if (associated(wig_in, wig_out)) return

! If both associated must be pointing to different memory locations

if (associated(wig_in) .and. associated(wig_out)) then
  wig_out%n_link = wig_out%n_link - 1
  if (wig_out%n_link == 0) then
    deallocate (wig_out%term)
    deallocate (wig_out)
  endif
  wig_out => wig_in
  wig_out%n_link = wig_out%n_link + 1

elseif (associated(wig_out)) then 
  wig_out%n_link = wig_out%n_link - 1
  if (wig_out%n_link == 0) then
    deallocate (wig_out%term)
    deallocate (wig_out)
  else
    nullify (wig_out)
  endif

elseif (associated(wig_in)) then 
  wig_out => wig_in
  wig_out%n_link = wig_out%n_link + 1
endif

end subroutine transfer_wig

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wall3d (wall3d_in, wall3d_out)
!
! Subroutine to point wall3d_out => wall3d_in
!
! Modules needed:
!   use bmad
!
! Input:
!   wall3d_in  -- Wall3d_struct, pointer: Input wall3dgler field.
!
! Output:
!   wall3d_out -- Wall3d_struct, pointer: Output wall3dgler field.
!-

subroutine transfer_wall3d (wall3d_in, wall3d_out)

implicit none

type (wall3d_struct), pointer :: wall3d_in, wall3d_out

!

if (.not. associated(wall3d_in) .and. .not. associated(wall3d_out)) return
if (associated(wall3d_in, wall3d_out)) return

! If both associated must be pointing to different memory locations

if (associated(wall3d_out)) call deallocate_wall3d_pointer(wall3d_out)

if (associated(wall3d_in)) then 
  wall3d_out => wall3d_in
  wall3d_out%n_link = wall3d_out%n_link + 1
endif

end subroutine transfer_wall3d

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_em_field (field_in, field_out)
!
! Subroutine to transfer the field info from one struct to another.
! In the end will have:
!     field_out%map  => field_in%map
!     field_out%grid => field_in%grid
!
! Modules needed:
!   use bmad
!
! Input:
!   field_in -- Field_struct, pointer: Input RF field.
!
! Output:
!   field_out -- Field_struct, pointer: Output RF field.
!-

subroutine transfer_em_field (field_in, field_out)

implicit none

type (em_fields_struct), pointer :: field_in, field_out
type (em_field_mode_struct), pointer :: mode, mode_in, mode_out

integer i

! Rule: If field_in or field_out is associated then %mode must be allocated

if (.not. associated(field_in) .and. .not. associated(field_out)) return

! field_in exists and field_out does not exist: Create field_out.

if (.not. associated(field_out)) then
  call init_em_field (field_out, size(field_in%mode))
  field_out%mode = field_in%mode
  do i = 1, size(field_out%mode)
    mode => field_out%mode(i)
    if (associated(mode%map)) mode%map%n_link = mode%map%n_link + 1
    if (associated(mode%grid)) mode%grid%n_link = mode%grid%n_link + 1
  enddo
  return
endif

! field_in does not exist and field_out exists: Deallocate field_out.

if (.not. associated(field_in)) then
  call init_em_field (field_out, 0)
  return
endif

! Both field_in and field_out exist: If both point to the same memory then need
! to do nothing. Otherwise need to transfer the data.

call init_em_field (field_out, size(field_in%mode))

do i = 1, size(field_out%mode)

  mode_in => field_in%mode(i)
  mode_out => field_out%mode(i)

  if (associated(mode_in%map) .and. associated(mode_out%map)) then
    if (.not. associated(mode_in%map, mode_out%map)) then
      mode_out%map%n_link = mode_out%map%n_link - 1
      if (mode_out%map%n_link == 0) deallocate (mode_out%map)
      mode_out%map => mode_in%map
      mode_out%map%n_link = mode_out%map%n_link + 1
    endif
  elseif (associated(mode_out%map) .and. .not. associated(mode_in%map)) then 
    mode_out%map%n_link = mode_out%map%n_link - 1
    if (mode_out%map%n_link == 0) deallocate (mode_out%map)
  elseif (associated(mode_in%map) .and. .not. associated(mode_out%map)) then 
    mode_out%map => mode_in%map
    mode_out%map%n_link = mode_out%map%n_link + 1
  endif

  if (associated(mode_in%grid) .and. associated(mode_out%grid)) then
    if (.not. associated(mode_in%grid, mode_out%grid)) then
      mode_out%grid%n_link = mode_out%grid%n_link - 1
      if (mode_out%grid%n_link == 0) deallocate (mode_out%grid)
      mode_out%grid => mode_in%grid
      mode_out%grid%n_link = mode_out%grid%n_link + 1
    endif
  elseif (associated(mode_out%grid) .and. .not. associated(mode_in%grid)) then 
    mode_out%grid%n_link = mode_out%grid%n_link - 1
    if (mode_out%grid%n_link == 0) deallocate (mode_out%grid)
  elseif (associated(mode_in%grid) .and. .not. associated(mode_out%grid)) then 
    mode_out%grid => mode_in%grid
    mode_out%grid%n_link = mode_out%grid%n_link + 1
  endif

  mode_out = mode_in

enddo

end subroutine transfer_em_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_em_field (em_field, n_mode)
!
! Subroutine to initialize a em_field_struct pointer.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_mode   -- Integer: Size of %modes(:) to create. If 0, deallocate em_field
!
! Output:
!   em_field -- em_field_struct, pointer: Initialized structure.
!-

subroutine init_em_field (em_field, n_mode)

type (em_fields_struct), pointer :: em_field
type (em_field_mode_struct), pointer :: mode

integer n_mode

integer i

! Cases where nothing is to be done

if (n_mode < 1 .and. .not. associated(em_field)) return

if (n_mode > 0 .and. associated(em_field)) then
  if (size(em_field%mode) == n_mode) return
endif

! Must deallocate existing.

if (associated(em_field)) then
  do i = 1, size(em_field%mode)
    mode => em_field%mode(i)
    if (associated(mode%map)) then
      mode%map%n_link = mode%map%n_link - 1
      if (mode%map%n_link == 0) deallocate (mode%map)
    endif
    if (associated(mode%grid)) then
      mode%grid%n_link = mode%grid%n_link - 1
      if (mode%grid%n_link == 0) deallocate (mode%grid)
    endif
  enddo
  deallocate(em_field)
endif
  
if (n_mode < 1) return

! n_mode > 0 case.

allocate(em_field)
allocate(em_field%mode(n_mode))

end subroutine init_em_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_rf_wake (wake_in, wake_out)
!
! Subroutine to transfer the wake info from one struct to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   wake_in -- Wake_struct, pointer: Input wake.
!
! Output:
!   wake_out -- Wake_struct, pointer: Output wake.
!-

subroutine transfer_rf_wake (wake_in, wake_out)

implicit none

type (rf_wake_struct), pointer :: wake_in, wake_out
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

!

if (associated (wake_in)) then
  n_sr_table       = size(wake_in%sr_table)
  n_sr_mode_long   = size(wake_in%sr_mode_long)
  n_sr_mode_trans  = size(wake_in%sr_mode_trans)
  n_lr             = size(wake_in%lr)
  call init_wake (wake_out, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
  wake_out%sr_file        = wake_in%sr_file
  wake_out%lr_file        = wake_in%lr_file
  wake_out%z_sr_mode_max  = wake_in%z_sr_mode_max
  wake_out%sr_table       = wake_in%sr_table
  wake_out%sr_mode_long   = wake_in%sr_mode_long
  wake_out%sr_mode_trans  = wake_in%sr_mode_trans
  wake_out%lr             = wake_in%lr
else
  if (associated(wake_out)) call init_wake (wake_out, 0, 0, 0, 0)
endif

end subroutine transfer_rf_wake

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_wake (wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
!
! Subroutine to initialize a wake struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_sr_table      -- Integer: Number of terms: wake%sr_table(0:n_sr-1).
!   n_sr_mode_long  -- Integer: Number of terms: wake%nr(n_sr_mode_long).
!   n_sr_mode_trans -- Integer: Number of terms: wake%nr(n_sr_mode_trans).
!   n_lr            -- Integer: Number of terms: wake%nr(n_lr)
!
! Output:
!   wake -- Wake_struct, pointer: Initialized structure. 
!               If all inputs are 0 then wake is deallocated.
!-

subroutine init_wake (wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

implicit none

type (rf_wake_struct), pointer :: wake
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

! Deallocate wake if all inputs are zero.

if (n_sr_table == 0 .and. n_sr_mode_long == 0 .and. n_sr_mode_trans == 0 .and. n_lr == 0) then
  if (associated(wake)) then
    deallocate (wake%sr_table)
    deallocate (wake%sr_mode_long)
    deallocate (wake%sr_mode_trans)
    deallocate (wake%lr)
    deallocate (wake)
  endif
  return
endif

!

if (associated (wake)) then
  if (size(wake%sr_table) /= n_sr_table) then
    deallocate (wake%sr_table)
    allocate (wake%sr_table(0:n_sr_table-1))
  endif
  if (size(wake%sr_mode_long) /= n_sr_mode_long) then
    deallocate (wake%sr_mode_long)
    allocate (wake%sr_mode_long(n_sr_mode_long))
  endif
  if (size(wake%sr_mode_trans) /= n_sr_mode_trans) then
    deallocate (wake%sr_mode_trans)
    allocate (wake%sr_mode_trans(n_sr_mode_trans))
  endif
  if (size(wake%lr) /= n_lr) then
    deallocate (wake%lr)
    allocate (wake%lr(n_lr))
  endif

else
  allocate (wake)
  allocate (wake%sr_table(0:n_sr_table-1))
  allocate (wake%sr_mode_long(n_sr_mode_long))
  allocate (wake%sr_mode_trans(n_sr_mode_trans))
  allocate (wake%lr(n_lr))
endif

end subroutine init_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine calc_superimpose_key (ele1, ele2, ele3, create_em_field_slave)
!
! Function to decide what ele3%key and ele3%sub_key should be
! when two elements, ele1, and ele2, are superimposed.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1     -- Ele_struct:
!     %key
!     %sub_key
!   ele2     -- Ele_struct:
!     %key
!     %sub_key
!   create_em_field_slave
!            -- Logical, optional: Default is False. If True then ele3%key
!                     will be set to em_field.
!
! Output:
!   ele3 -- Ele_struct:
!     %key        -- Set to -1 if there is an error
!     %sub_key
!-

subroutine calc_superimpose_key (ele1, ele2, ele3, create_em_field_slave)

implicit none

type (ele_struct), target :: ele1, ele2, ele3
integer key1, key2
integer, pointer :: key3
logical, optional :: create_em_field_slave

!

key1 = ele1%key
key2 = ele2%key
key3 => ele3%key

key3 = -1  ! Default if no superimpse possible
ele3%sub_key = 0

! control elements cannot be superimposed.

if (any(key1 == [overlay$, group$, girder$])) return
if (any(key2 == [overlay$, group$, girder$])) return

! Superimposing two of like kind...

if (key1 == key2) then
  select case (key1)
  case (sbend$)
    ! Bad
  case (wiggler$, rfcavity$)
    key3 = em_field$
    ele3%sub_key = const_ref_energy$
  case (lcavity$)
    key3 = em_field$
    ele3%sub_key = nonconst_ref_energy$
  case (em_field$)
    key3 = em_field$
    if (ele1%sub_key == nonconst_ref_energy$ .or. ele2%sub_key == nonconst_ref_energy$) then
      ele3%sub_key = nonconst_ref_energy$
    else
      ele3%sub_key = const_ref_energy$
    endif
  case default
    key3 = key1
  end select
  return
endif

! If one element is a drift then key3 = key of other element.

if (key1 == drift$) then
  key3 = key2
  ele3%sub_key = ele2%sub_key
  return
endif

if (key2 == drift$) then
  key3 = key1
  ele3%sub_key = ele1%sub_key
  return
endif

! If one element is a pipe then key3 = key of other element.

if (any(key1 == [pipe$])) then
  key3 = key2
  return
endif

if (any(key2 == [pipe$])) then
  key3 = key1
  return
endif

! If one element is a rcollimator, monitor, or instrument then key3 = key of other element.

if (ele1%aperture_type == elliptical$ .or. ele2%aperture_type == elliptical$ ) ele3%aperture_type = elliptical$

if (any(key1 == [ecollimator$, rcollimator$, monitor$, instrument$])) then
  key3 = key2
  return
endif

if (any(key2 == [ecollimator$, rcollimator$, monitor$, instrument$])) then
  key3 = key1
  return
endif

! If one element is a kicker then key3 = key of other element.

if (any(key1 == [kicker$, hkicker$, vkicker$])) then
  if (any(key2 == [kicker$, hkicker$, vkicker$])) then
    key3 = kicker$
  else
    key3 = key2
  endif
  return
endif

if (any(key2 == [kicker$, hkicker$, vkicker$])) then
  key3 = key1
  return
endif

! General case...

! sbend elements are problematical due to the different reference orbit so cannot superimpose them.

if (key1 == sbend$ .or. key2 == sbend$) return

! em_field wanted

if (logic_option(.false., create_em_field_slave)) then
  key3 = em_field$
  if (key1 == lcavity$ .or. key2 == lcavity$) then
    ele3%sub_key = nonconst_ref_energy$
  elseif (key1 == em_field$) then
    ele3%sub_key = ele1%sub_key
  elseif (key2 == em_field$) then
    ele3%sub_key = ele2%sub_key
  else
    ele3%sub_key = const_ref_energy$
  endif
  return
endif

!

select case (key1)

case (quadrupole$,  solenoid$, sol_quad$) 
  select case (key2)
  case (quadrupole$);    key3 = sol_quad$
  case (solenoid$);      key3 = sol_quad$
  case (sol_quad$);      key3 = sol_quad$
  case (bend_sol_quad$); key3 = bend_sol_quad$
  case (sbend$);         key3 = bend_sol_quad$
  end select

case (bend_sol_quad$)
  select case (key2)
  case (quadrupole$);    key3 = bend_sol_quad$
  case (solenoid$);      key3 = bend_sol_quad$
  case (sol_quad$);      key3 = bend_sol_quad$
  case (sbend$);         key3 = bend_sol_quad$
  end select
end select

if (key3 /= -1) return  ! Have found something

! Only thing left is to use em_field type element.

key3 = em_field$
if (key1 == lcavity$ .or. key2 == lcavity$) then
  ele3%sub_key = nonconst_ref_energy$
elseif (key1 == em_field$) then
  ele3%sub_key = ele1%sub_key
elseif (key2 == em_field$) then
  ele3%sub_key = ele2%sub_key
else
  ele3%sub_key = const_ref_energy$
endif

end subroutine calc_superimpose_key

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

implicit none

type (ele_struct) ele
type (lat_param_struct) param
real(rp) grad_shift

! ele%value(grad_loss_sr_wake$) is an internal variable used with macroparticles.
! It accounts for the longitudinal short-range wakefields between macroparticles.
! It is continually being modified for each macroparticle

if (bmad_com%sr_wakes_on .and. ele%value(l$) /= 0) then
  grad_shift = ele%value(e_loss$) * param%n_part * abs(charge_of(param%particle)) * &
                                      e_charge / ele%value(l$) - ele%value(grad_loss_sr_wake$) 
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
!                      attribute_group$, control_group$, floor_position_group$,
!                      s_position_group$, ref_energy_group$, or mat6_group$, all_groups$
!   set_slaves    -- Logical, optional: If present and False then do not set
!                      the status for any slaves. Default is True.
!-

recursive subroutine set_ele_status_stale (ele, status_group, set_slaves)

implicit none

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
  call set_attributes
  call set_mat6

case (control_group$)
  call set_control

case (floor_position_group$)
  call set_floor_position
  call set_mat6

case (s_position_group$)
  call set_s_position
  call set_floor_position
  call set_mat6

case (ref_energy_group$)
  call set_ref_energy
  call set_mat6
  call set_attributes ! EG: k1 <--> b1_gradient calc needed 

case (mat6_group$)
  call set_mat6

case (rad_int_group$)
  call set_rad_int

case (all_groups$)
  call set_attributes
  call set_control
  call set_ref_energy
  call set_floor_position
  call set_s_position
  call set_mat6
  call set_rad_int

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

subroutine set_attributes
  if (ele%lord_status == overlay_lord$) return
  if (ele%lord_status == group_lord$) return
  ele%bookkeeping_state%attributes = stale$
  if (associated(state)) state%attributes = stale$
end subroutine set_attributes

!----------------------------------------------------------------------------
! contains

subroutine set_control
  if (ele%lord_status == not_a_lord$ .and. ele%slave_status == free$) return
  ele%bookkeeping_state%control = stale$
  if (associated(state)) state%control = stale$
end subroutine set_control

!----------------------------------------------------------------------------
! contains

subroutine set_floor_position
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%floor_position = stale$
  if (associated(state)) state%floor_position = stale$
end subroutine set_floor_position

!----------------------------------------------------------------------------
! contains

subroutine set_s_position
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%s_position = stale$
  if (associated(state)) state%s_position = stale$
end subroutine set_s_position

!----------------------------------------------------------------------------
! contains

subroutine set_ref_energy
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%ref_energy = stale$
  if (associated(state)) state%ref_energy = stale$
end subroutine set_ref_energy

!----------------------------------------------------------------------------
! contains

subroutine set_rad_int
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%rad_int = stale$
  if (associated(state)) state%rad_int = stale$
end subroutine set_rad_int

!----------------------------------------------------------------------------
! contains

! Ignore if the element does not have an associated linear transfer map
! Also the branch status does not get set since the transfer map calc
! must always check a branch due to possible reference orbit shifts.

subroutine set_mat6
  if (ele%lord_status == overlay_lord$) return
  if (ele%lord_status == group_lord$) return
  if (ele%lord_status == multipass_lord$) return
  ele%bookkeeping_state%mat6 = stale$
end subroutine set_mat6

end subroutine set_ele_status_stale 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_status_flags (bookkeeping_state, stat)
!
! Routine to set the bookkeeping status block.
!
! Input:
!   stat          -- Integer: bookkeeping status. ok$, stale$, etc.
!
! Output:
!   bookkeeping_state -- bookkeeping_state_struct: 
!-

subroutine set_status_flags (bookkeeping_state, stat)

implicit none

type (bookkeeping_state_struct) bookkeeping_state
integer stat

!

bookkeeping_state%control        = stat
bookkeeping_state%s_position     = stat
bookkeeping_state%floor_position = stat
bookkeeping_state%ref_energy     = stat
bookkeeping_state%attributes     = stat
bookkeeping_state%mat6           = stat
bookkeeping_state%rad_int        = stat

end subroutine set_status_flags

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

implicit none

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
! Function pointer_to_branch_given_name (branch_name, lat) result (branch_ptr)
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
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the nameed branch.
!                   Nullified if there is no such branch.
!-

function pointer_to_branch_given_name (branch_name, lat) result (branch_ptr)

implicit none

type (branch_struct), pointer :: branch_ptr
type (lat_struct), target :: lat

integer i, ib, ios
character(*) branch_name
character(32), parameter :: r_name = 'pointer_to_branch_given_name'

! Init in case of error

nullify(branch_ptr)

! Is index.

if (is_integer(trim(branch_name))) then
  read (branch_name, *, iostat = ios) ib
  if (ib < 0 .or. ib > ubound(lat%branch, 1)) then
    call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: ' // branch_name)
    return
  endif
  branch_ptr => lat%branch(ib)

! Is name.

else
  do i = lbound(lat%branch, 1), ubound(lat%branch, 1)
    if (lat%branch(i)%name == branch_name) then
      branch_ptr => lat%branch(ib)
      return
    endif
  enddo
  call out_io (s_error$, r_name, 'BRANCH NAME NOT FOUND: ' // branch_name)
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

implicit none

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

implicit none

type (ele_struct), target :: lord, slave
type (ele_struct), pointer :: ele
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch

integer ic_out, icon_out, ib, i, n

!

if (.not. associated(lord%branch, slave%branch)) call err_exit  ! Should not be

lat => lord%branch%lat

! Find lat%control(:) and lat%ic(:) elements associated with this link

do ic_out = slave%ic1_lord, slave%ic2_lord
  icon_out = lat%ic(ic_out)
  if (lat%control(icon_out)%ix_lord == lord%ix_ele) exit
enddo

if (icon_out == slave%ic2_lord + 1) call err_exit ! Should not be

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
    if (ele%ix2_slave >= icon_out) ele%ix2_slave = ele%ix2_slave - 1

    if (ele%ic1_lord >  ic_out) ele%ic1_lord = ele%ic1_lord - 1
    if (ele%ic2_lord >= ic_out) ele%ic2_lord = ele%ic2_lord - 1
  enddo
enddo

slave%n_lord = slave%n_lord - 1
if (slave%n_lord == 0) then
  slave%ic1_lord = 0
  slave%ic2_lord = -1
  slave%slave_status = free$
endif

lord%n_slave = lord%n_slave - 1
if (lord%n_slave == 0) then
  lord%ix1_slave = 0
  lord%ix2_slave = -1
  lord%lord_status = not_a_lord$
endif

end subroutine remove_lord_slave_link 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_slave (lord, ix_slave, ix_contrl) result (slave_ptr)
!
! Function to point to a slave of a lord.
! Also see:
!   pointer_to_lord
!   pointer_to_ele
!   pointer_to_field_ele
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   lord     -- Ele_struct: Lord element
!   ix_slave -- Integer: Index of the slave. ix_slave goes from 1 to lord%n_slave
!
! Output:
!   slave_ptr  -- Ele_struct, pointer: Pointer to the slave.
!                   Nullified if there is an error.
!   ix_control -- Integer, optional :: index of appropriate lat%control(:) element.
!                   Set to -1 is there is an error.
!-

function pointer_to_slave (lord, ix_slave, ix_control) result (slave_ptr)

implicit none

type (ele_struct), target :: lord
type (ele_struct), pointer :: slave_ptr
type (control_struct), pointer :: con
type (lat_struct), pointer :: lat

integer, optional :: ix_control
integer ix_slave, icon

!

if (ix_slave > lord%n_slave .or. ix_slave < 1) then
  nullify(slave_ptr)
  if (present(ix_control)) ix_control = -1
  return
endif

lat => lord%branch%lat
icon = lord%ix1_slave + ix_slave - 1
con => lat%control(icon)
slave_ptr => lat%branch(con%ix_branch)%ele(con%ix_slave)
if (present(ix_control)) ix_control = icon

end function pointer_to_slave

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_lord (slave, ix_lord, ix_control, ix_slave) result (lord_ptr)
!
! Function to point to a lord of a slave.
! Also see:
!   pointer_to_slave
!   pointer_to_ele
!   pointer_to_field_ele
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   slave      -- Ele_struct: Slave element.
!   ix_lord    -- Integer: Index of the lord. ix_lord goes from 1 to slave%n_lord
!
! Output:
!   lord_ptr   -- Ele_struct, pointer: Pointer to the lord.
!                   Nullified if there is an error.
!   ix_control -- Integer, optional: Index of appropriate lat%control(:) element.
!                   Set to -1 is there is an error or the slave is a slice_slave.
!   ix_slave   -- Integer, optional: Index of back to the slave. That is, 
!                   pointer_to_slave(lord_ptr, ix_slave) will point back to slave. 
!                   Set to -1 is there is an error or the slave is a slice_slave.
!-

function pointer_to_lord (slave, ix_lord, ix_control, ix_slave) result (lord_ptr)

implicit none

type (ele_struct), target :: slave
type (ele_struct), pointer :: lord_ptr
type (lat_struct), pointer :: lat

integer, optional :: ix_control, ix_slave
integer i, ix_lord, icon

! Case where there is no lord

if (ix_lord > slave%n_lord .or. ix_lord < 1) then
  nullify(lord_ptr)
  if (present(ix_control)) ix_control = -1
  if (present(ix_slave)) ix_slave = -1
  return
endif

! slice_ele stores info differently

lat => slave%branch%lat

if (slave%slave_status == slice_slave$) then
  lord_ptr => slave%lord
  if (present(ix_control)) ix_control = -1
  if (present(ix_slave)) ix_slave = -1
  return
endif

! Point to the lord

icon = lat%ic(slave%ic1_lord + ix_lord - 1)
lord_ptr => lat%ele(lat%control(icon)%ix_lord)

if (present(ix_control)) ix_control = icon

! There must be a corresponding ix_slave value such that
!   pointer_to_slave(lord_ptr, ix_slave) => slave 

if (present(ix_slave)) then

  do i = 1, lord_ptr%n_slave
    if (associated (pointer_to_slave(lord_ptr, i), slave)) then
      ix_slave = i
      return
    endif
  enddo

  ! If ix_slave not found then this is an error
  if (global_com%exit_on_error) call err_exit   

endif

end function pointer_to_lord

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function num_field_eles (ele_in) result (n_field_eles)
!
! Function to return the number of elements that have field information for
! the ele_in lattice_element.
!
! This routine is to be used in conjunction with pointer_to_field_ele.
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   ele_in     -- Ele_struct: Lattice element.
!
! Output:
!   n_field_eles -- Integer: Number of elements that have field information.
!-

function num_field_eles (ele_in) result (n_field_eles)

implicit none

type (ele_struct) ele_in
integer n_field_eles

!

select case (ele_in%slave_status)

case (super_slave$, multipass_slave$)
  n_field_eles = ele_in%n_lord

case default
  n_field_eles = 1

end select


end function num_field_eles

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_field_ele (ele_in, ix_field) result (field_ele_ptr)
!
! Routine to return a pointer to an element that contains field information
! for the ele_in lattice element.
!
! An element, if it is a super_slave or multipass_slave, will have field information 
! stored in the lord elements.
! In this case, pointer_to_field_ele is identical to pointer_to_lord.
!
! If the element is not a super_slave or a multipass_slave, the element will 
! itself contain the field information.
! In this case pointer_to_field_ele (ele_in, 1) will return a pointer to itself.
!
! Note: Use the function num_field_eles(ele) to obtain the number of elements where
! field information is stored.
!
! Also see:
!   pointer_to_slave
!   pointer_to_ele
!   pointer_to_lord
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   ele_in     -- Ele_struct: Lattice element.
!   ix_field   -- Integer: Index to select which of the field info elements to point to.
!                   ix_field should go from 1 to num_field_eles(ele).
!
! Output:
!   field_ele_ptr -- Ele_struct, pointer: Pointer to an element containing field information.
!                      Nullified if there is an error.
!-

function pointer_to_field_ele (ele_in, ix_field) result (field_ele_ptr)

implicit none

type (ele_struct), target :: ele_in
type (ele_struct), pointer :: field_ele_ptr
integer ix_field

!

select case (ele_in%slave_status)

case (super_slave$, multipass_slave$)
  field_ele_ptr => pointer_to_lord(ele_in, ix_field)

case default
  if (ix_field == 1) then
    field_ele_ptr => ele_in
  else
    nullify (field_ele_ptr)
  endif

end select

end function pointer_to_field_ele

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_next_ele (this_ele, dir) result (next_ele)
!
! Function to return a pointer to the next element in the lattice branch.
! next_ele will be nullified if there is a problem like this_ele is
! the last element of a linear lattice.
!
! Input:
!   this_ele  -- ele_struct: Starting element.
!   dir       -- integer, optional: If positive return next forward element.
!                   If negative then return previous element. Default = +1.
!
!   next_ele -- ele_struct, pointer: Element after this_ele.
!-

function pointer_to_next_ele (this_ele, dir) result (next_ele)

implicit none

type (ele_struct), target :: this_ele
type (ele_struct), pointer :: next_ele
type (ele_struct), pointer :: an_ele
type (branch_struct), pointer :: branch

integer, optional :: dir

!

next_ele => null()

if (.not. associated(this_ele%branch)) return

branch => this_ele%branch
if (this_ele%ix_ele < 0 .or. this_ele%ix_ele > branch%n_ele_max) return

if (this_ele%ix_ele > branch%n_ele_track) then
  if (this_ele%lord_status /= super_lord$) return
  an_ele => pointer_to_slave(this_ele, this_ele%n_slave)
else
  an_ele => this_ele
endif

if (integer_option(+1, dir) > 0) then
  if (an_ele%ix_ele == branch%n_ele_track) then
    if (branch%param%geometry == open$) return
    next_ele => branch%ele(0)
  else
   next_ele => branch%ele(an_ele%ix_ele+1)
  endif

else
  if (an_ele%ix_ele == 0) then
    if (branch%param%geometry == open$) return
    next_ele => branch%ele(branch%n_ele_track)
  else
   next_ele => branch%ele(an_ele%ix_ele-1)
  endif
endif

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

implicit none

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

implicit none

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

implicit none

type (ele_struct) ele
logical is_const

! Anything with longitudinal electric fields or anything
! where the "zero-orbit" is not a straight line down the middle
! has a varying ds/dt(ref).

select case (ele%key)
case (lcavity$, custom$, hybrid$, wiggler$, rfcavity$, em_field$)
  is_const = .false.
case default
  is_const = .true.
end select

end function ele_has_constant_ds_dt_ref

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function tracking_uses_end_drifts (ele) result (has_drifts)
!
! Function to determine if the tracking for an element uses a "hard edge model"
! where the tracking looks like (drift, model, drift). For example,
! RF cavity fields with ele%field_calc = bmad_standard$ use a hard edge model where
! the length of the cavity is c_light / (2 * freq).
!
! Module needed:
!   use bmad
!
! Input:
!   ele    -- ele_struct: Element.
!
! Output:
!   has_drifts -- Logical: True if tracking uses end drifts.
!-

function tracking_uses_end_drifts (ele) result (has_drifts)

implicit none

type (ele_struct) ele
logical has_drifts

!

has_drifts = .false.
if (.not. bmad_com%use_hard_edge_drifts) return

select case (ele%key)
case (lcavity$, rfcavity$, solenoid$)
    if (ele%field_calc == bmad_standard$) has_drifts = .true.
end select

end function tracking_uses_end_drifts

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function element_has_fringe_fields (ele) result (has_fringe)
!
! Function to determine if the element has fringe fields that must be accounted
! for when tracking through the element. For example, a solenoid has edge fields.
! Note: If an element has end drifts (see tracking_uses_end_drifts function), 
! Then the fringe_fields will not be internal to the element and not at the ends.
!
! Module needed:
!   use bmad
!
! Input:
!   ele    -- ele_struct: Element.
!
! Output:
!   has_fringe -- Logical: True if the element has fringe fields.
!-

function element_has_fringe_fields (ele) result (has_fringe)

implicit none

type (ele_struct) ele
logical has_fringe

!

has_fringe = .false.

select case (ele%key)
case (lcavity$, rfcavity$, solenoid$, sbend$, sol_quad$, bend_sol_quad$, e_gun$)
    if (ele%field_calc == bmad_standard$) has_fringe = .true.
end select

end function element_has_fringe_fields

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine calc_next_fringe_edge (track_ele, s_edge_track, hard_ele, s_edge_hard, hard_end)
!
! Routine to locate the next "hard edge" in an element when a hard edge model is being used. 
! This routine is used by integration tracking routines like Runge-Kutta.
! This routine is called repeatedly as the integration routine tracks through the element.
! If the element is a super_slave, there are potentially many hard edges.
!
! Rule: When track_ele is a super_slave, the edges of track_ele may be inside the field region.
! In this case, the hard edge is applied at the edges to make the particle look like it is
! outside the field. This is done to make the tracking symplectic from entrance edge to exit
! edge. [Remember that Bmad's coordinates are not canonical inside a field region.]
!
! Input:
!   track_ele     -- ele_struct: Element being tracked through.
!   hard_ele      -- ele_struct, pointer: Needs to be nullified at the start of tracking.
!
! Output:
!   s_edge_track -- Real(rp): S position of next hard edge in track_ele frame.
!                     If there are no more hard edges then s_pos will be set to ele%value(l$).
!   hard_ele     -- ele_struct, pointer: Points to element with the hard edge.
!                     Will be nullified if there is no hard edge.
!                     This will be track_ele unless track_ele is a super_slave.
!   s_edge_hard  -- Real(rp): S-position of next hard egde in hard_ele frame.
!   hard_end     -- Integer: Describes hard edge. Set to upstream_end$ or downstream_end$.
!-

subroutine calc_next_fringe_edge (track_ele, s_edge_track, hard_ele, s_edge_hard, hard_end)

implicit none

type (ele_struct), target :: track_ele
type (ele_struct), pointer :: hard_ele, lord

real(rp) s_edge_track, s_edge_hard
integer hard_end
integer i

! Init if needed.
! Keep track of things by setting: 
!   ele%ixx = 0  ! Init setting
!   ele%ixx = 1  ! About to track through entrance hard edge
!   ele%ixx = 2  ! About to track through exit hard edge

if (.not. associated(hard_ele)) then
  if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
    do i = 1, track_ele%n_lord
      lord => pointer_to_lord(track_ele, i)
      lord%ixx = 0
    enddo
  else
    track_ele%ixx = 0
  endif
endif

! Find next hard edge

s_edge_track = track_ele%value(l$)
nullify (hard_ele)

if (track_ele%slave_status == super_slave$ .or. track_ele%slave_status == slice_slave$) then
  do i = 1, track_ele%n_lord
    lord => pointer_to_lord(track_ele, i)
    call does_this_ele_contain_the_next_edge (lord)
  enddo
else
  call does_this_ele_contain_the_next_edge (track_ele)
endif

if (associated (hard_ele)) hard_ele%ixx = hard_ele%ixx + 1

!-------------------------------------------------------------------------
contains 

subroutine does_this_ele_contain_the_next_edge (this_ele)

type (ele_struct), target :: this_ele
real(rp) s_this_edge, s_hard_entrance, s_hard_exit, s_off
integer this_end

!

if (.not. element_has_fringe_fields (this_ele)) return
if (this_ele%ixx == 2) return

s_off = (this_ele%s - this_ele%value(l$)) - (track_ele%s - track_ele%value(l$))
s_hard_entrance = s_off + (this_ele%value(l$) - hard_edge_model_length(this_ele)) / 2 
s_hard_exit     = s_off + (this_ele%value(l$) + hard_edge_model_length(this_ele)) / 2 

if (s_hard_entrance < -bmad_com%significant_length .and. &
    s_hard_exit     < -bmad_com%significant_length) return

if (s_hard_entrance > track_ele%value(l$) + bmad_com%significant_length .and. &
    s_hard_exit     > track_ele%value(l$) + bmad_com%significant_length) return

if (this_ele%ixx == 0) then
  this_end = upstream_end$
  s_this_edge = max(0.0_rp, s_hard_entrance)

else   ! this_ele%ixx = 1
  this_end = downstream_end$
  s_this_edge = min(track_ele%value(l$), s_hard_exit)
endif

if (s_this_edge > s_edge_track) return

! This looks like the next hard edge

hard_ele => this_ele
hard_end = this_end
s_edge_track = s_this_edge
s_edge_hard = s_edge_track - s_off

end subroutine does_this_ele_contain_the_next_edge

end subroutine calc_next_fringe_edge

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

implicit none

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
drift_ele%value(l$)              = (ele_in%value(l$) - hard_edge_model_length(drift_ele)) / 2 
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

implicit none

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

if (ele%key == em_field$ .and. ele%slave_status == super_slave$) then
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
! Function particle_time (orbit, ele) result (time)
!
! Routine to return the current time for use, for example, 
! in calculations of time-dependent EM fields.
!
! The oscillations of such fields are synched relative to the absolute clock if 
! absolute time tracking is used or are synched relative to the reference particle
! if relative time tracking is used.
!
! Input:
!   orbit -- Coord_struct: Particle coordinates
!   ele   -- ele_struct: Element being tracked through.
!
! Ouput:
!   time  -- Real(rp): Current time.
!-

function particle_time (orbit, ele) result (time)

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) time
logical abs_time
character(16), parameter :: r_name = 'particle_time'

! Note: e_gun uses absolute time tracking to get around the problem when orbit%beta = 0.

if (absolute_time_tracking(ele)) then
  time = orbit%t 

else
  if (orbit%beta == 0) then
    call out_io (s_fatal$, r_name, 'PARTICLE IN NON E-GUN ELEMENT HAS VELOCITY = 0!')
    if (global_com%exit_on_error) call err_exit
    time = orbit%t  ! Just to keep on going
    return
  endif
  time = -orbit%vec(5) / (orbit%beta * c_light)
endif

end function particle_time

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function slave_time_offset (ele) result (time)
!
! Routine to return the reference time at the start of a slave element relative
! to the the reference time at the beginning of the lord.
! 
! Use of this routine is ambiguous if there are multiple lords.
!
! Input:
!   ele   -- ele_struct: Element being tracked through.
!
! Ouput:
!   time  -- Real(rp): Reference time as start of element relative to the lord.
!             Returns 0 if element is not a slave.
!-

function slave_time_offset (ele) result (time)

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (ele_struct), pointer :: lord
real(rp) time
logical abs_time
character(16), parameter :: r_name = 'slave_time_offset'

! If not a slave then return 0

if (ele%slave_status /= super_slave$ .and. ele%slave_status /= slice_slave$) then
  time = 0
  return
endif

!

lord => pointer_to_lord (ele, 1)
time = ele%value(ref_time_start$) - lord%value(ref_time_start$)

end function slave_time_offset

end module
