!+
! This file defines the interfaces for the BMAD subroutines
!-

module bmad_routine_interface

use bmad_struct

implicit none

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

interface assignment (=)
  module procedure ele_equal_ele
  module procedure ele_vec_equal_ele_vec
  module procedure lat_equal_lat 
  module procedure lat_vec_equal_lat_vec 
  module procedure branch_equal_branch
  module procedure taylor_equal_taylor
  module procedure taylors_equal_taylors
  module procedure em_taylor_equal_em_taylor
  module procedure em_taylors_equal_em_taylors
  module procedure complex_taylor_equal_complex_taylor
  module procedure complex_taylors_equal_complex_taylors
  module procedure bunch_equal_bunch
  module procedure beam_equal_beam
end interface

interface operator (+)
  module procedure em_field_plus_em_field
end interface

interface operator (*)
  module procedure map1_times_map1
end interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch
!
! Routine to return a pointer to the lattice branch associated with a given name
! or a given element.
!
! This routine is an overloaded name for:
!   pointer_to_branch_given_ele (ele) result (branch_ptr))
!   pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)
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
! Input:
!   ele                  -- ele_struct: Element contained in the branch.
!   branch_name          -- character(*): May be a branch name or a branch index.
!   lat                  -- lat_struct: Lattice to search.
!   parameter_is_branch0 -- logical, optional: If True, 'PARAMETER' is taken to be
!                             an alternative name for branch(0). Default is False.
!   blank_branch         -- integer, optional: Branch index if branch_name = ''. Default is blank is an error.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the branch.
!                   Nullified if there is no associated branch.
!-

interface pointer_to_branch
  function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)
    import
    implicit none
    type (branch_struct), pointer :: branch_ptr
    type (lat_struct), target :: lat
    integer, optional :: blank_branch
    logical, optional :: parameter_is_branch0
    character(*) branch_name
  end function

  recursive function pointer_to_branch_given_ele (ele) result (branch_ptr)
    import
    implicit none
    type (ele_struct), target :: ele
    type (branch_struct), pointer :: branch_ptr
  end function
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
!     Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
!     Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
!     Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)
!
! pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
! lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
! element in lat.
! 
! Note that using ele_name to locate an element is potentially dangerous if there
! are multiple elements that have the same name. Better in this case is to use:
!   lat_ele_locator
!
! Also see:
!   pointer_to_slave
!   pointer_to_lord
!
! Input:
!   lat           -- lat_struct: Lattice.
!   ix_ele        -- integer: Index of element in lat%branch(ix_branch).
!   ix_branch     -- integer: Index of the lat%branch(:) containing the element.
!   ix_nametable  -- integer: Nametable index. See above
!   ele_loc       -- lat_ele_loc_struct: Location identification.
!   ele_name      -- character(*): Name or index of element.
!   foreign_ele   -- ele_struct: Lattice element in another lattice.
!
! Output:
!   ele_ptr       -- ele_struct, pointer: Pointer to the element. Nullified if no match or error.
!-

interface pointer_to_ele
  module procedure pointer_to_ele1
  module procedure pointer_to_ele2
  module procedure pointer_to_ele3
  module procedure pointer_to_ele4
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
  subroutine reallocate_coord_n (coord, n_coord)
    import
    implicit none
    type (coord_struct), allocatable :: coord(:)
    integer, intent(in) :: n_coord
  end subroutine

  subroutine reallocate_coord_lat (coord, lat, ix_branch)
    import
    implicit none
    type (coord_struct), allocatable :: coord(:)
    type (lat_struct), target :: lat
    integer, optional :: ix_branch
  end subroutine

  subroutine reallocate_coord_array (coord_array, lat)
    import
    implicit none
    type (coord_array_struct), allocatable :: coord_array(:)
    type (lat_struct) lat
  end subroutine
end interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord (...)
!
! Routine to initialize a coord_struct. 
!
! This routine is an overloaded name for:
!   Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
!   Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
!   Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
!
! Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice), 
! or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6) 
! corresponding to the beginning_ele's value of ele%value(p0c_start$).
!
! Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead 
! of the standard which is to set orb%t from orb%vec(5).
!
! For photons:
!   orb%vec(5) is set depending upon where the photon is relative to the element.
!   If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle 
!       except if direction is set.
!
! Input:
!   orb_in       -- coord_struct: Input orbit.
!   vec(6)       -- real(rp), optional: Coordinate vector. If not present then taken to be zero.
!   ele          -- ele_struct, optional: Particle is initialized to start at element_end of this ele.
!   element_end  -- integer, optional: upstream_end$, downstream_end$, inside$, or start_end$.
!                     Must be present if ele argument is present.
!                     start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1.
!                     Default is upstream_end$. Note: If ele is the beginning element (index zero), the
!                     setting of element_end will not matter.
!   particle     -- integer, optional: Particle type (electron$, etc.). 
!                     If particle = not_set$ and orb_in is present, use orb_in%species instead.
!   dirction     -- integer, optional: +1 -> moving downstream +s direciton, -1 -> moving upstream.
!                     0 -> Ignore. Default is to not change orb%direction except for photons which get set
!                     according to orb%vec(6).
!   E_photon     -- real(rp), optional: Photon energy if particle is a photon. Ignored otherwise.
!   t_offset     -- real(rp), optional: Offset of the reference time. This is non-zero when
!                     there are multiple bunches and the reference time for a particular particle
!                     is pegged to the time of the center of the bunch.
!   shift_vec6   -- logical, optional: If present and False, prevent the shift of orb%vec(6).
!   spin(3)      -- real(rp), optional: Particle spin. Taken to be zero if not present.
!   s_pos        -- real(rp), optional: Particle s-position. Only relavent if element_end = inside$.
!   random_on    -- logical, optional: Default is True. Used only for photons being initalized with a photon_init
!                     element. If True, vary the photon coords using a random number generator. If False, the 
!                     photon coords will be centered within the distribution specified in the photon_init ele.
!
! Output:
!   orb -- Coord_struct: Initialized coordinate.
!                 Note: For photons, orb%vec(6) is computed as sqrt(1 - vec(2)^2 - vec(4)^2) if needed.
!-

interface init_coord
  subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
    import
    implicit none
    type (coord_struct) orb
    type (ele_struct), optional :: ele
    real(rp) :: vec(6)
    real(rp), optional :: t_offset, E_photon, spin(3), s_pos
    integer, optional :: element_end, particle, direction
    logical, optional :: shift_vec6, random_on
  end subroutine

  subroutine init_coord2 (orb_out, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
    import
    implicit none
    type (coord_struct) orb_out, orb_in
    type (ele_struct), optional, target :: ele
    real(rp), optional :: E_photon, t_offset, spin(3), s_pos
    integer, optional :: element_end, particle, direction
    logical, optional :: shift_vec6, random_on
  end subroutine

  subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)
    import
    implicit none
    type (coord_struct) orb
    type (ele_struct), optional :: ele
    real(rp), optional :: t_offset, E_photon, spin(3)
    integer, optional :: element_end, particle, direction
    logical, optional :: shift_vec6
  end subroutine
end interface

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_attribute (...)
!
! Routine to mark an element or lattice as modified for use with "intelligent" bookkeeping.
! Also will do some dependent variable bookkeeping when a particular attribute has 
! been altered.
!
! This routine should be called after the attribute has been set.
!
! set_flags_for_changed_attribute is an overloaded name for:
!   set_flags_for_changed_lat_attribute (lat, set_dependent)
!   set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
!   set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
!   set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
!   set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)
!
! The set_flags_for_changed_lat_attribute (lat) routine is used when one
! does not know what has changed and wants a complete bookkeeping done.
!
! NOTE: The attribute argument MUST be the component that was changed. For example:
!     ele%value(x_offset$) = off_value
!     call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
! And NOT:
!     call set_flags_for_changed_attribute (ele, off_value)  ! WRONG
!
! Input:
!   lat           -- lat_struct: Lattice being modified.
!   ele           -- ele_struct, Element being modified.
!   real_attrib   -- real(rp), optional: Attribute that has been changed.
!                      For example: ele%value(hkick$).
!                      If not present then assume everything has potentially changed.
!   int_attrib    -- integer: Attribute that has been changed.
!                      For example: ele%mat6_calc_method.
!   logic_attrib  -- logical; Attribute that has been changed.
!                      For example: ele%is_on.
!   all_attrib    -- all_pointer_struct: Pointer to attribute.
!   set_dependent -- logical, optional: If False then dependent parameter bookkeeping will not be done.
!                     False is used, for example, during parsing when dependent bookkeepin is not wanted.
!                     Default is True. Do not set False unless you know what you are doing.
!
! Output:
!   lat  -- lat_struct: Lattice with appropriate changes.
!-

interface set_flags_for_changed_attribute
  subroutine set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)
    import
    implicit none
    type (ele_struct), target :: ele
    type (all_pointer_struct) all_attrib
    logical, optional :: set_dependent
  end subroutine

  subroutine set_flags_for_changed_integer_attribute (ele, attrib, set_dependent)
    import
    implicit none
    type (ele_struct), target :: ele
    integer, target :: attrib
    logical, optional :: set_dependent
  end subroutine

  subroutine set_flags_for_changed_logical_attribute (ele, attrib, set_dependent)
    import
    implicit none
    type (ele_struct), target :: ele
    logical, target :: attrib
    logical, optional :: set_dependent
  end subroutine

  subroutine set_flags_for_changed_lat_attribute (lat, set_dependent)
    import
    implicit none
    type (lat_struct), target :: lat
    logical, optional :: set_dependent
  end subroutine

  subroutine set_flags_for_changed_real_attribute (ele, attrib, set_dependent)
    import
    implicit none
    type (ele_struct), target :: ele
    real(rp), optional, target :: attrib
    logical, optional :: set_dependent
  end subroutine
end interface


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

interface

function absolute_time_tracking (ele) result (is_abs_time)
  import
  implicit none
  type (ele_struct), target :: ele
  logical is_abs_time
end function

function ac_kicker_amp(ele, orbit, true_time) result (ac_amp)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  real(rp), optional :: true_time
  real(rp) ac_amp
end function

subroutine add_lattice_control_structs (ele, n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field, add_at_end)
  import
  implicit none
  type (ele_struct) ele
  integer, optional :: n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field
  logical, optional :: add_at_end
end subroutine

subroutine allocate_branch_array (lat, upper_bound)
  import
  implicit none
  type (lat_struct), target :: lat
  integer :: upper_bound
end subroutine

subroutine allocate_element_array (ele, upper_bound)
  import
  implicit none
  type (ele_struct), pointer :: ele(:)
  integer, optional :: upper_bound
end subroutine

subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch, do_ramper_slave_setup)
  import
  implicit none
  type (lat_struct), target :: lat
  integer, optional :: upper_bound
  integer, optional :: ix_branch
  logical, optional :: do_ramper_slave_setup
end subroutine

subroutine aml_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

function angle_between_polars (polar1, polar2) result (angle)
  import
  implicit none
  type (spin_polar_struct), intent(in) :: polar1, polar2
  real(rp) :: angle
end function

subroutine angle_to_canonical_coords (orbit, coord_type)
  import
  implicit none
  type (coord_struct) orbit
  character(*), optional :: coord_type
end subroutine

subroutine apply_all_rampers (lat, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  logical err_flag
end subroutine

subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix, rf_time, apply_sol_fringe)
  import
  implicit none
  type (coord_struct) orb
  type (fringe_field_info_struct) fringe_info
  type (ele_struct), target :: track_ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6), rf_time
  logical, optional :: make_matrix, apply_sol_fringe
  logical track_spin
end subroutine

subroutine apply_energy_kick (dE, orbit, ddE_dr, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  real(rp) dE, ddE_dr(2)
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

recursive subroutine apply_rampers_to_slave (slave, err_flag)
  import
  implicit none
  type (ele_struct), target :: slave
  logical err_flag
end subroutine

function at_this_ele_end (now_at, where_at) result (is_at_this_end)
  import
  implicit none
  integer now_at, where_at
  logical is_at_this_end
end function

subroutine attribute_bookkeeper (ele, force_bookkeeping)
  import
  implicit none
  type (ele_struct), target :: ele
  logical, optional :: force_bookkeeping
end subroutine

subroutine attribute_set_bookkeeping (ele, attrib_name, err_flag, attrib_ptr)
  import
  implicit none
  type (ele_struct) ele
  type (all_pointer_struct), optional :: attrib_ptr
  character(*) attrib_name
  logical err_flag
end subroutine

subroutine autoscale_phase_and_amp(ele, param, err_flag, scale_phase, scale_amp, call_bookkeeper)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  logical, optional :: scale_phase, scale_amp, call_bookkeeper
end subroutine

function average_twiss (frac1, twiss1, twiss2) result (ave_twiss)
  import
  implicit none
  type (twiss_struct) twiss1, twiss2, ave_twiss
  real(rp) frac1
end function

subroutine bbi_kick (x, y, sigma, nk, dnk)
  import
  implicit none
  real(rp) x, y, sigma(2), nk(2), dnk(2,2)
end subroutine

subroutine bbi_slice_calc (ele, n_slice, z_slice)
  import
  implicit none
  type (ele_struct) ele
  integer :: n_slice
  real(rp) z_slice(:)
end subroutine

subroutine bend_exact_multipole_field (ele, param, orbit, local_ref_frame, field, calc_dfield, calc_potential)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  type (em_field_struct) field
  logical local_ref_frame
  logical, optional :: calc_dfield, calc_potential
end subroutine

function bend_length_has_been_set(ele) result (is_set)
  import
  implicit none
  type (ele_struct) ele
  logical is_set
end function

function bend_shift (position1, g, delta_s, w_mat, ref_tilt) result(position2)
  import
  implicit none
  type (floor_position_struct) :: position1, position2
  real(rp) :: g, delta_s, S_mat(3,3), L_vec(3), tlt, angle
  real(rp), optional :: w_mat(3,3), ref_tilt
end function bend_shift

subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag, parse_lat)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  type (lat_struct), optional :: parse_lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine bmad_parser2 (in_file, lat, orbit, make_mats6, err_flag, parse_lat)
  import
  implicit none
  character(*) in_file
  type (lat_struct), target :: lat
  type (lat_struct), optional :: parse_lat
  type (coord_struct), optional :: orbit(0:)
  logical, optional :: make_mats6, err_flag
end subroutine

function branch_name(branch) result (name)
  import
  implicit none
  type (branch_struct), target :: branch
  character(40) name
end function

subroutine ramper_slave_setup(lat, do_setup)
  import
  implicit none
  type (lat_struct), target :: lat
  logical, optional :: do_setup
end subroutine

function ramper_value (ramper, r1, err_flag) result (value)
  import
  implicit none
  type (ele_struct) ramper
  type (control_ramp1_struct) r1
  real(rp) value
  logical err_flag
end function

subroutine remove_dead_from_bunch(bunch_in, bunch_out)
  import
  implicit none
  type (bunch_struct) bunch_in, bunch_out
end subroutine

function c_multi (n, m, no_n_fact, c_full) result (c_out)
  import
  implicit none
  integer, intent(in) :: n, m
  real(rp) c_out
  real(rp), optional :: c_full(0:n_pole_maxx, 0:n_pole_maxx)
  logical, optional :: no_n_fact
end function

subroutine c_to_cbar (ele, cbar_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp) cbar_mat(2,2)
end subroutine

subroutine calc_next_fringe_edge (track_ele, s_edge_body, fringe_info, orbit, init_needed, time_tracking)
  import
  type (ele_struct), target :: track_ele
  type (fringe_field_info_struct) fringe_info
  type (coord_struct) :: orbit
  real(rp) s_edge_body
  logical, optional :: init_needed, time_tracking
end subroutine

subroutine calc_super_slave_key (lord1, lord2, slave, create_jumbo_slave)
  import
  implicit none
  type (ele_struct), target :: lord1, lord2, slave
  logical, optional :: create_jumbo_slave
end subroutine

subroutine calc_z_tune (branch)
  import
  implicit none
  type (branch_struct), target :: branch
end subroutine

subroutine canonical_to_angle_coords (orbit, coord_type)
  import
  implicit none
  type (coord_struct) orbit
  character(*), optional :: coord_type
end subroutine

subroutine cbar_to_c (cbar_mat, a, b, c_mat)
  import
  implicit none
  real(rp) cbar_mat(2,2), c_mat(2,2)
  type (twiss_struct) a, b
end subroutine

recursive subroutine check_aperture_limit (orb, ele, particle_at, param, old_orb, check_momentum)
  import
  implicit none
  type (coord_struct) :: orb
  type (ele_struct) :: ele
  type (lat_param_struct), intent(inout) :: param
  integer particle_at
  type (coord_struct), optional :: old_orb
  logical, optional :: check_momentum
end subroutine

function classical_radius (species) result (radius)
  import
  integer species
  real(rp) radius
end function

function coords_body_to_rel_exit (body_position, ele, w_mat, calculate_angles) result(rel_exit)
  import
  implicit none
  type (floor_position_struct) :: body_position, rel_exit
  type (ele_struct) :: ele
  real(rp), optional :: w_mat(3,3)
  logical, optional :: calculate_angles
end function 

function coords_body_to_local (body_position, ele, w_mat, calculate_angles) result(local_position)
  import
  implicit none
  type (floor_position_struct) :: body_position, local_position
  type (ele_struct) :: ele
  real(rp), optional :: w_mat(3,3)
  logical, optional :: calculate_angles
end function coords_body_to_local

function coords_relative_to_floor (floor0, dr, theta, phi, psi) result (floor1)
  import
  implicit none
  type (floor_position_struct) floor0, floor1
  real(rp) dr(3)
  real(rp), optional :: theta, phi, psi
end function coords_relative_to_floor

function coords_floor_to_relative (floor0, global_position, calculate_angles, is_delta_position) result (local_position)
  import
  implicit none
  type (floor_position_struct) floor0, global_position, local_position
  logical, optional :: calculate_angles, is_delta_position
end function coords_floor_to_relative

function coords_floor_to_local_curvilinear (global_position, ele, status, w_mat, relative_to) result(local_position)
  import
  implicit none
  type (floor_position_struct) :: global_position, local_position
  type (ele_struct)   :: ele
  real(rp), optional :: w_mat(3,3)
  integer :: status
  integer, optional :: relative_to
end function coords_floor_to_local_curvilinear

function coords_floor_to_curvilinear (floor_coords, ele0, ele1, status, w_mat) result (local_coords)
  import
  implicit none
  type (floor_position_struct) floor_coords, local_coords
  type (ele_struct), target :: ele0
  type (ele_struct), pointer :: ele1
  integer status
  real(rp), optional :: w_mat(3,3)
end function coords_floor_to_curvilinear

function coords_local_curvilinear_to_body (local_position, ele, w_mat, calculate_angles) result (body_position)
  import
  implicit none
  type (floor_position_struct) :: local_position, body_position, p, floor0
  type (ele_struct) :: ele
  real(rp), optional :: w_mat(3,3)
  logical, optional :: calculate_angles
end function

function coords_local_curvilinear_to_floor (local_position, ele, in_body_frame, &
                                                w_mat, calculate_angles, relative_to) result (global_position)
  import
  implicit none
  type (floor_position_struct) :: local_position, global_position
  type (ele_struct), target :: ele
  real(rp), optional :: w_mat(3,3)
  logical, optional :: in_body_frame
  logical, optional :: calculate_angles
  integer, optional :: relative_to
end function coords_local_curvilinear_to_floor

function coords_curvilinear_to_floor (xys, branch, err_flag) result (global)
  import
  implicit none
  type (branch_struct), target :: branch
  type (floor_position_struct) global, local
  real(rp) xys(3)
  logical err_flag
end function coords_curvilinear_to_floor

subroutine check_controller_controls (ele_key, contrl, name, err)
  import
  implicit none
  type (control_struct), target :: contrl(:)
  integer ele_key
  logical err
  character(*) name
end subroutine

subroutine check_if_s_in_bounds (branch, s, err_flag, translated_s, print_err)
  import
  implicit none
  type (branch_struct) branch
  real(rp) s
  real(rp), optional :: translated_s
  logical err_flag
  logical, optional :: print_err
end subroutine

subroutine choose_quads_for_set_tune (branch, dk1, eles, mask, err_flag)
  import
  implicit none
  type (branch_struct), target :: branch
  type (ele_pointer_struct), allocatable :: eles(:)
  character(*), optional :: mask
  real(rp), allocatable :: dk1(:)
  logical, optional :: err_flag
end subroutine

subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
                       pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch, orb0)
  import
  implicit none
  type (lat_struct), target :: lat
  type (lat_struct), optional, target :: low_E_lat, high_E_lat
  type (coord_struct), allocatable, optional, target :: low_E_orb(:), high_E_orb(:)
  type (coord_struct), optional :: orb0
  real(rp) delta_e
  real(rp) chrom_x
  real(rp) chrom_y
  real(rp), optional :: pz
  logical, optional, intent(out) :: err_flag
  integer, optional :: ix_branch
end subroutine

subroutine chrom_tune (lat, delta_e, chrom_x, chrom_y, err_tol, err_flag)
  import
  implicit none
  type (lat_struct) lat
  real(rp) delta_e
  real(rp) chrom_x
  real(rp) chrom_y
  real(rp) err_tol
  logical err_flag
end subroutine

subroutine clear_lat_1turn_mats (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine clear_taylor_maps_from_elements (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag, print_err)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), allocatable, target :: closed_orb(:)
  integer, optional :: direction, ix_branch, i_dim
  logical, optional, intent(out) :: err_flag
  logical, optional, intent(in) :: print_err
end subroutine

subroutine closed_orbit_from_tracking (lat, closed_orb, i_dim, eps_rel, eps_abs, init_guess, err_flag)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable :: closed_orb(:)
  type (coord_struct), optional :: init_guess
  real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)
  integer i_dim
  logical, optional :: err_flag
end subroutine

subroutine combine_consecutive_elements (lat, error)
  import
  implicit none
  type (lat_struct), target :: lat
  logical error
end subroutine

subroutine control_bookkeeper (lat, ele, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  type (ele_struct), optional :: ele
  logical, optional :: err_flag
end subroutine

subroutine convert_particle_coordinates_s_to_t (particle, s_body, orientation)
  import
  implicit none
  type (coord_struct), intent(inout), target :: particle
  real(rp) s_body
  integer :: orientation
end subroutine

subroutine convert_particle_coordinates_t_to_s (particle, ele, s_body, use_downstream_p0c)
  import
  implicit none
  type (coord_struct), target :: particle
  type (ele_struct) :: ele
  real(rp), optional :: s_body
  logical, optional :: use_downstream_p0c
end subroutine

subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, beta1, err_flag, print_err)
  import
  implicit none
  real(rp), intent(in) :: E_tot
  real(rp), intent(out), optional :: pc, kinetic, beta, brho, gamma, beta1
  integer, intent(in) :: particle
  logical, optional :: err_flag, print_err
end subroutine

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, beta1, err_flag)
  import
  implicit none
  real(rp), intent(in) :: pc
  real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma, beta1
  integer, intent(in) :: particle
  logical, optional :: err_flag
end subroutine

subroutine convert_bend_exact_multipole (g, out_type, an, bn)
  import
  implicit none
  real(rp) g, an(0:n_pole_maxx), bn(0:n_pole_maxx)
  integer out_type
end subroutine

subroutine create_feedback(lord, input, output, err)
  import
  implicit none
  type (ele_struct), target :: lord
  character(*) input(:), output(:)
  logical err
end subroutine

recursive subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset, &
                       param, include_upstream_end, include_downstream_end, err_flag, old_slice, orb_in)
  import
  implicit none
  type (ele_struct), target :: sliced_ele, ele_in
  type (ele_struct), optional :: old_slice
  type (coord_struct), optional :: orb_in
  type (lat_param_struct) param
  real(rp) l_slice, offset
  logical include_upstream_end, include_downstream_end, err_flag
end subroutine

subroutine create_field_overlap (lat, lord_name, slave_name, err_flag)
  import
  implicit none
  type (lat_struct) lat
  character(*) lord_name, slave_name
  logical err_flag
end subroutine

subroutine create_girder (lat, ix_ele, con, init_ele, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  type (ele_struct) :: init_ele
  type (control_struct) con(:)
  integer, intent(in) :: ix_ele
  logical err_flag
end subroutine

subroutine create_group (lord, con, err)
  import
  implicit none
  type (ele_struct), target :: lord
  type (control_struct) con(:)
  logical err
end subroutine

subroutine create_lat_ele_nametable (lat, nametable)
  import
  implicit none
  type (lat_struct), target :: lat
  type (nametable_struct), target :: nametable
end subroutine

subroutine create_overlay (lord, contl, err)
  import
  implicit none
  type (ele_struct), target :: lord
  type (control_struct) contl(:)
  logical err
end subroutine

subroutine create_ramper (lord, contl, err)
  import
  implicit none
  type (ele_struct), target :: lord
  type (control_struct), target :: contl(:)
  logical err
end subroutine

subroutine create_unique_ele_names (lat, key, suffix)
  import
  implicit none
  type (lat_struct), target :: lat
  integer key
  character(*) suffix
end subroutine

subroutine create_wiggler_cartesian_map (ele, cart_map)
  import
  implicit none
  type (ele_struct) ele
  type (cartesian_map_struct), target :: cart_map
end subroutine

subroutine crystal_attribute_bookkeeper (ele)
  import
  implicit none
  type (ele_struct), target :: ele
end subroutine

subroutine convert_coords (in_type_str, coord_in, ele, out_type_str, coord_out, err_flag)
  import
  implicit none
  character(*) in_type_str
  character(*) out_type_str
  type (coord_struct) coord_in
  type (coord_struct) coord_out
  type (ele_struct) ele
  logical, optional :: err_flag
end subroutine

subroutine deallocate_ele_array_pointers (eles)
  import
  implicit none
  type (ele_struct), pointer :: eles(:)
end subroutine

subroutine deallocate_ele_pointers (ele, nullify_only, nullify_branch, dealloc_poles)
  import
  implicit none
  type (ele_struct), target :: ele
  logical, optional, intent(in) :: nullify_only, nullify_branch, dealloc_poles
end subroutine

subroutine deallocate_lat_pointers (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

function default_tracking_species (param) result (species)
  import
  implicit none
  type (lat_param_struct) param
  integer species
end function

function diffraction_plate_or_mask_hit_spot (ele, orbit) result (ix_section)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  integer :: ix_section
end function

recursive function distance_to_aperture (orbit, particle_at, ele, no_aperture_here) result (dist)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) dist
  integer particle_at
  logical no_aperture_here
end function

subroutine do_mode_flip (ele, err_flag)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: err_flag
end subroutine

function e_accel_field (ele, voltage_or_gradient, bmad_standard_tracking) result (field)
  import
  implicit none
  type (ele_struct) ele
  real(rp) field
  integer voltage_or_gradient
  logical, optional :: bmad_standard_tracking
end function

recursive subroutine ele_compute_ref_energy_and_time (ele0, ele, param, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele0, ele
  type (lat_param_struct) param
  real(rp) e_tot_start, p0c_start, ref_time_start
  logical err_flag
end subroutine

recursive subroutine ele_geometry (floor_start, ele, floor_end, len_scale, ignore_patch_err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (floor_position_struct), optional, target :: floor_end
  type (floor_position_struct) :: floor_start
  real(rp), optional :: len_scale
  logical, optional :: ignore_patch_err
end subroutine ele_geometry

function ele_geometry_with_misalignments (ele, len_scale) result (floor)
  import
  implicit none
  type (ele_struct), target :: ele
  type (floor_position_struct) floor
  real(rp), optional :: len_scale
end function ele_geometry_with_misalignments

function ele_has_constant_ds_dt_ref (ele) result (is_const)
  import
  implicit none
  type (ele_struct) ele
  logical is_const
end function

function ele_has_nonzero_kick (ele) result (has_kick)
  import
  implicit none
  type (ele_struct) ele
  logical has_kick
end function

function ele_has_nonzero_offset (ele) result (has_offset)
  import
  implicit none
  type (ele_struct) ele
  logical has_offset
end function

function ele_loc_name (ele, show_branch0, parens) result (str)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: show_branch0
  character(2), optional :: parens
  character(10) str
end function

function ele_full_name (ele, template) result (str)
  import
  implicit none
  type (ele_struct) ele
  character(*), optional :: template
  character(:), allocatable :: str
end function

subroutine ele_misalignment_L_S_calc (ele, L_mis, S_mis)
  import
  implicit none
  type(ele_struct) :: ele 
  real(rp) :: L_mis(3), S_mis(3,3)
end subroutine ele_misalignment_L_S_calc

function ele_loc (ele) result (loc)
  import
  implicit none
  type (ele_struct) ele
  type (lat_ele_loc_struct) loc
end function

function ele_nametable_index(ele) result(ix_nt)
  import
  implicit none
  type (ele_struct), target :: ele
  integer ix_nt
end function

subroutine ele_order_calc (lat, order)
  import
  implicit none
  type (lat_struct), target :: lat
  type (lat_ele_order_struct) order
end subroutine

subroutine ele_rad_int_cache_calc (ele)
  import
  implicit none
  type (ele_struct) ele
end subroutine

subroutine ele_to_fibre (ele, ptc_fibre, use_offsets, err_flag, integ_order, steps, for_layout, ref_in)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct), optional :: ref_in
  type (fibre), pointer :: ptc_fibre
  integer, optional :: integ_order, steps
  logical use_offsets, err_flag
  logical, optional :: for_layout
end subroutine

subroutine ele_to_spin_taylor(ele, param, orb0)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orb0
end subroutine

subroutine ele_to_taylor (ele, orb0, taylor_map_includes_offsets, include_damping, orbital_taylor, spin_taylor)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct), optional, intent(in) :: orb0
  type (taylor_struct), optional, target :: orbital_taylor(6), spin_taylor(0:3)
  logical, optional :: taylor_map_includes_offsets, include_damping
end subroutine

function ele_unique_name (ele, order) result (unique_name)
  import
  implicit none
  type (ele_struct) ele
  type (lat_ele_order_struct) order
  character(40) unique_name
end function ele_unique_name

function ele_value_has_changed (ele, list, abs_tol, set_old) result (has_changed)
  import
  implicit none
  type (ele_struct) ele
  integer list(:)
  real(rp) abs_tol(:)
  logical set_old, has_changed
end function

subroutine elec_multipole_field (a, b, n, coord, Ex, Ey, dE, compute_dE)
  import
  implicit none
  type (coord_struct)  coord
  real(rp) a, b
  real(rp), optional :: dE(2,2)
  real(rp) Ex, Ey
  integer n
  logical, optional :: compute_dE
end subroutine

subroutine element_slice_iterator (ele, param, i_slice, n_slice_tot, sliced_ele, s_start, s_end)
  import
  implicit none
  type (ele_struct) ele, sliced_ele
  type (lat_param_struct) param
  integer i_slice, n_slice_tot
  real(rp), optional :: s_start, s_end
end subroutine

recursive subroutine em_field_calc (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, err_flag, &
             calc_potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles, print_err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) :: orbit
  type (em_field_struct) :: field
  type (ele_pointer_struct), allocatable, optional :: used_eles(:)
  real(rp) s_pos
  real(rp), optional :: rf_time
  logical :: local_ref_frame
  logical, optional :: calc_dfield, calc_potential, err_flag, use_overlap, grid_allow_s_out_of_bounds, print_err
end subroutine

function entering_element(orbit, particle_at) result (is_entering)
  import
  implicit none
  type (coord_struct) orbit
  integer particle_at
  logical is_entering
end function

function equivalent_taylor_attributes (ele_taylor, ele2) result (equiv)
  import
  implicit none
  type (ele_struct) :: ele_taylor, ele2
  logical equiv
end function

subroutine fibre_to_ele (ptc_fibre, branch, ix_ele, err_flag, from_mad)
  import
  implicit none
  type (fibre), target :: ptc_fibre
  type (branch_struct) branch
  integer ix_ele
  logical err_flag
  logical, optional :: from_mad
end subroutine

subroutine find_element_ends (ele, ele1, ele2, ix_multipass)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: ele1, ele2
  integer, optional :: ix_multipass
end subroutine

subroutine find_matching_fieldmap (file_name, ele, t_type, match_ele, ix_field, ignore_slaves)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: match_ele
  integer t_type, ix_field
  logical, optional :: ignore_slaves
  character(*) file_name
end subroutine

subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat, w_mat_inv)
  import
  implicit none
  real(rp), optional :: w_mat(3,3), w_mat_inv(3,3)
  real(rp) theta, phi, psi
end subroutine floor_angles_to_w_mat

subroutine floor_w_mat_to_angles (w_mat, theta, phi, psi, floor0)
  import
  implicit none
  type (floor_position_struct), optional :: floor0
  real(rp) theta, phi, psi, w_mat(3,3)
end subroutine floor_w_mat_to_angles

function fringe_here (ele, orbit, particle_at) result (is_here)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orbit
  integer particle_at
  logical is_here
end function

subroutine g_bending_strength_from_em_field (ele, param, s_rel, orbit, local_ref_frame, g, dg)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  real(rp), intent(in) :: s_rel
  real(rp), intent(out) :: g(3)
  real(rp), optional :: dg(3,3)
  logical local_ref_frame
end subroutine

function gamma_ref(ele) result (gamma)
  import
  implicit none
  type (ele_struct) ele
  real(rp) gamma
end function

subroutine gen_grad_at_s_to_em_taylor (ele, gen_grad, s_pos, em_taylor)
  import
  implicit none
  type (ele_struct) ele
  type (gen_grad_map_struct), target :: gen_grad
  type (em_taylor_struct), target :: em_taylor(3)
  real(rp) s_pos
end subroutine

subroutine gen_grad1_to_em_taylor (ele, gen_grad, iz, em_taylor)
  import
  implicit none
  type (ele_struct) ele
  type (gen_grad_map_struct), target :: gen_grad
  type (em_taylor_struct), target :: em_taylor(3)
  integer iz
end subroutine

subroutine get_slave_list (lord, slaves, n_slave)
  import
  implicit none
  type (ele_struct), target :: lord
  type (ele_pointer_struct), allocatable :: slaves(:)
  integer n_slave
end subroutine

function gradient_shift_sr_wake (ele, param) result (grad_shift)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp) grad_shift
end function

subroutine hdf5_read_beam (file_name, beam, error, ele, pmd_header, print_p0c_shift_warning, conserve_momentum)
  import
  implicit none
  type (beam_struct), target :: beam
  type (ele_struct), optional :: ele
  type (pmd_header_struct), optional :: pmd_header
  logical, optional :: print_p0c_shift_warning, conserve_momentum
  logical error
  character(*) file_name
end subroutine

subroutine hdf5_read_grid_field (file_name, ele, g_field, err_flag, pmd_header, combine)
  import
  implicit none
  type (grid_field_struct), pointer :: g_field(:)
  type (ele_struct) ele
  type (pmd_header_struct), optional :: pmd_header
  logical, optional :: combine
  logical err_flag
  character(*) file_name
end subroutine

subroutine hdf5_write_beam (file_name, bunches, append, error, lat, alive_only)
  import
  implicit none
  type (bunch_struct), target :: bunches(:)
  type (lat_struct), optional :: lat
  logical error, append
  logical, optional :: alive_only
  character(*) file_name
end subroutine

subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)
  import
  implicit none
  type (grid_field_struct), target :: g_field(:)
  type (ele_struct) ele
  logical err_flag
  character(*) file_name
end subroutine

subroutine g_integrals_calc (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine init_photon_from_a_photon_init_ele (ele, param, orbit, random_on)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  logical, optional :: random_on
end subroutine

subroutine init_bmad()
  import
  implicit none
end subroutine

subroutine init_bmad_parser_common(lat)
  import
  implicit none
  type (lat_struct), optional:: lat
end subroutine

subroutine init_custom (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine init_ele (ele, key, sub_key, ix_ele, branch)
  import
  implicit none
  type (ele_struct)  ele
  type (branch_struct), optional, target :: branch
  integer, optional :: key, sub_key
  integer, optional :: ix_ele
end subroutine

subroutine init_fringe_info (fringe_info, ele, orbit, leng_sign)
  import
  implicit none
  type (fringe_field_info_struct) fringe_info
  type (ele_struct) ele
  type (coord_struct), optional :: orbit
  integer, optional :: leng_sign
end subroutine

subroutine init_lat (lat, n, init_beginning_ele)
  import
  implicit none
  type (lat_struct), target :: lat
  integer, optional :: n
  logical, optional :: init_beginning_ele
end subroutine

subroutine init_multipole_cache(ele)
  import
  implicit none
  type (ele_struct) ele
end subroutine

subroutine init_taylor_series (bmad_taylor, n_term, save_old)
  import
  implicit none
  type (taylor_struct) bmad_taylor
  integer n_term
  logical, optional :: save_old
end subroutine

subroutine init_wake (wake, n_sr_long, n_sr_trans, n_sr_z, n_lr_mode, always_allocate)
  import
  implicit none
  type (wake_struct), pointer :: wake
  integer n_sr_long, n_sr_trans, n_sr_z, n_lr_mode
  logical, optional :: always_allocate
end subroutine

subroutine insert_element (lat, insert_ele, insert_index, ix_branch, orbit)
  import
  implicit none
  type (lat_struct), target :: lat
  type (ele_struct) insert_ele
  integer insert_index
  integer, optional :: ix_branch
  type (coord_struct), optional, allocatable :: orbit(:)
end subroutine

subroutine ion_kick (orbit, r_beam, n_beam_part, a_twiss, b_twiss, sig_ee, kick)
  import
  implicit none
  type (coord_struct) orbit
  type (twiss_struct) a_twiss, b_twiss
  real(rp) r_beam(2), n_beam_part, sig_ee, kick(3)
end subroutine

function key_name_to_key_index (key_str, abbrev_allowed) result (key_index)
  import
  implicit none
  character(*) key_str
  logical, optional :: abbrev_allowed
  integer key_index
end function

subroutine kill_ptc_layouts (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine kill_taylor (bmad_taylor)
  import
  implicit none
  type (taylor_struct) :: bmad_taylor(:)
end subroutine

function knot_interpolate (x_knot, y_knot, x_pt, interpolation, err_flag) result (y_pt)
  import
  implicit none
  real(rp) x_knot(:), y_knot(:), x_pt, y_pt
  integer interpolation
  logical err_flag
end function

function knots_to_string (x_knot, y_knot) result (str)
  import
  implicit none
  real(rp) x_knot(:), y_knot(:)
  character(:), allocatable :: str
end function

subroutine lat_compute_ref_energy_and_time (lat, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  logical err_flag
end subroutine

subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch, order_by_index)
  import
  implicit none
  character(*) loc_str
  type (lat_struct), target :: lat
  type (ele_pointer_struct), allocatable :: eles(:)
  integer n_loc
  logical, optional :: above_ubound_is_err, err, order_by_index
  integer, optional :: ix_dflt_branch
end subroutine

subroutine lat_geometry (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine lat_geometry

recursive subroutine lat_make_mat6 (lat, ix_ele, ref_orb, ix_branch, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), optional :: ref_orb(0:)
  integer, optional :: ix_ele, ix_branch
  logical, optional :: err_flag
end subroutine

subroutine lat_sanity_check (lat, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  logical, intent(out) :: err_flag
end subroutine

subroutine lat_to_ptc_layout (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine lattice_bookkeeper (lat, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  logical, optional :: err_flag
end subroutine

subroutine linear_to_spin_taylor(q_map, spin_taylor)
  import
  type (taylor_struct) spin_taylor(0:3)
  real(rp) q_map(0:3, 0:6)
end subroutine

function lord_edge_aligned (slave, slave_edge, lord) result (is_aligned)
  import
  implicit none
  type (ele_struct), target :: slave, lord
  integer slave_edge
  logical is_aligned
end function

function low_energy_z_correction (orbit, ele, ds, mat6, make_matrix) result (dz)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) ds, dz
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end function

subroutine make_g_mats (ele, g_mat, g_inv_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp) g_mat(4,4)
  real(rp) g_inv_mat(4,4)
end subroutine

subroutine make_g2_mats (twiss, g2_mat, g2_inv_mat)
  import
  implicit none
  type (twiss_struct) twiss
  real(rp) g2_mat(2,2), g2_inv_mat(2,2)
end subroutine

subroutine make_hybrid_lat (r_in, r_out, use_taylor, orb0_arr)
  import
  implicit none
  type (lat_struct), target :: r_in
  type (lat_struct), target :: r_out
  logical, optional :: use_taylor
  type (coord_array_struct), optional :: orb0_arr(0:)
end subroutine

function map1_inverse (map1) result (inv_map1)
  import
  implicit none
  type (spin_orbit_map1_struct) map1, inv_map1
end function

subroutine map1_make_unit(map1)
  import
  implicit none
  type (spin_orbit_map1_struct) map1
end subroutine

recursive subroutine make_mat6 (ele, param, start_orb, end_orb, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct), optional :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err_flag
end subroutine

subroutine make_mat6_taylor (ele, start_orb, end_orb, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  logical, optional :: err_flag
end subroutine

subroutine make_mat6_bmad (ele, param, start_orb, end_orb, err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err
end subroutine

subroutine make_mat6_bmad_photon (ele, param, start_orb, end_orb, err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: err
end subroutine

subroutine make_mat6_runge_kutta (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_mat6_symp_lie_ptc (ele, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
end subroutine

subroutine make_mat6_tracking (ele, param, start_orb, end_orb, err_flag, spin_only)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical err_flag
  logical, optional :: spin_only
end subroutine

subroutine make_v_mats (ele, v_mat, v_inv_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp), optional :: v_mat(4,4)
  real(rp), optional :: v_inv_mat(4,4)
end subroutine

subroutine map_to_angle_coords (t_canon, t_angle)
  import
  implicit none
  type (taylor_struct) t_canon(6), t_angle(6)
end subroutine

function master_parameter_value (master_parameter, ele) result (value)
  import
  implicit none
  type (ele_struct) ele
  real(rp) value
  integer master_parameter
end function

subroutine mat6_add_offsets (ele, param)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
end subroutine

subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)
  import
  implicit none
  real(rp) mat6(6,6), x_pitch_tot, y_pitch_tot
  integer orientation
end subroutine

subroutine mat_symp_decouple(t0, stat, U, V, Ubar, Vbar, G,  twiss1, twiss2, gamma, type_out)
  import
  implicit none
  real(rp) t0(4,4), U(4,4), V(4,4)
  integer stat
  real(rp) Ubar(4,4), Vbar(4,4), G(4,4)
  type (twiss_struct)  twiss1, twiss2
  real(rp) gamma
  logical type_out
end subroutine

subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)
  import
  implicit none
  real(rp) knl, tilt
  integer n
  type (coord_struct) orbit
  real(rp) kick_mat(4,4)
end subroutine

subroutine match_ele_to_mat6 (ele, start_orb, mat6, vec0, err_flag, include_delta_time, set_trombone)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) start_orb
  real(rp) mat6(6,6), vec0(6)
  logical :: err_flag
  logical, optional :: include_delta_time, set_trombone
end subroutine

function mexp (x, m) result (this_exp)
  import
  implicit none
  real(rp) x, this_exp
  integer m
end function

function momentum_compaction(branch) result (mom_comp)
  import
  implicit none
  type (branch_struct), target :: branch
  real(rp) mom_comp
end function

subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele, stable, growth_rate, chi, err_flag)
  import
  implicit none
  type (coord_struct), intent(in) :: track(:)
  type (coord_struct), intent(out) :: track0
  type (ele_struct) :: ele
  real(rp), intent(out) :: growth_rate, chi
  integer, intent(in) :: i_dim
  logical, intent(out) :: stable, err_flag
end subroutine

subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, map0, track0, chi)
  import
  implicit none
  type (coord_struct), intent(in), target :: track(:)
  type (coord_struct), intent(out) :: track0
  real(rp), intent(out) :: mat1(:,:), map0(:)
  real(rp), intent(out) :: chi
  integer, intent(in) :: i_dim
end subroutine

subroutine multipass_all_info (lat, info)
  import
  implicit none
  type (lat_struct), target :: lat
  type (multipass_all_info_struct), target :: info
end subroutine

subroutine multipass_chain (ele, ix_pass, n_links, chain_ele, use_super_lord)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_pointer_struct), allocatable, optional :: chain_ele(:)
  integer :: n_links, ix_pass
  logical, optional :: use_super_lord
end subroutine

subroutine multipole1_ab_to_kt (an, bn, n, knl, tn)
  import
  implicit none
  real(rp) an, bn
  real(rp) knl, knsl, tn
  integer n
end subroutine

subroutine multipole1_kt_to_ab (knl, knsl, tn, n, an, bn)
  import
  implicit none
  real(rp) an, bn
  real(rp) knl, knsl, tn
  integer n
end subroutine

subroutine multipole_ab_to_kt (an, bn, knl, tn)
  import
  implicit none
  real(rp) an(0:), bn(0:)
  real(rp) knl(0:), tn(0:)
end subroutine

recursive subroutine multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks, b1)
  import
  implicit none
  type (ele_struct), target :: ele
  real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rp), optional :: b1
  integer ix_pole_max
  integer, optional :: pole_type, include_kicks
  integer include_kck
  logical use_ele_tilt
end subroutine

subroutine multipole_ele_to_kt (ele, use_ele_tilt, ix_pole_max, knl, tilt, pole_type, include_kicks)
  import
  implicit none
  type (ele_struct), target :: ele
  real(rp) knl(0:), tilt(0:)
  integer ix_pole_max
  integer, optional :: pole_type, include_kicks
  logical use_ele_tilt
end subroutine

subroutine multipole_kt_to_ab (knl, knsl, tn, an, bn)
  import
  implicit none
  real(rp) an(0:), bn(0:)
  real(rp) knl(0:), knsl(0:), tn(0:)
end subroutine

subroutine multipole_init (ele, who, zero)
  import
  implicit none
  type (ele_struct) ele
  integer who
  logical, optional :: zero
end subroutine

subroutine multipole_kick_mat (knl, tilt, ref_species, ele, orbit, factor, mat6)
  import
  implicit none
  real(rp) knl(0:), tilt(0:)
  integer ref_species
  type (ele_struct) ele
  type (coord_struct) orbit
  real(rp) mat6(6,6), factor
end subroutine

subroutine multipole_spin_tracking (ele, param, orbit)
  import
  implicit none
  type (ele_struct) :: ele
  type (lat_param_struct) param
  type (coord_struct) orbit
end subroutine

subroutine new_control (lat, ix_ele, ele_name)
  import
  implicit none
  type (lat_struct) lat
  integer ix_ele
  character(*), optional :: ele_name
end subroutine

subroutine normal_mode_dispersion(ele, reverse) 
  import
  type (ele_struct) ele
  logical, optional :: reverse
end subroutine

function num_field_eles (ele) result (n_field_ele)
  import
  implicit none
  type (ele_struct) ele
  integer n_field_ele
end function

function num_lords (slave, lord_type) result (num)
  import
  implicit none
  type (ele_struct), target :: slave
integer lord_type, num
end function

subroutine offset_particle (ele, set, coord, set_tilt, set_hvkicks, drift_to_edge, &
                                        s_pos, s_out, set_spin, mat6, make_matrix, spin_qrot, time)
  import
  implicit none
  type (ele_struct) :: ele
  type (coord_struct), intent(inout) :: coord
  integer, optional :: drift_to_edge
  logical, intent(in) :: set
  logical, optional, intent(in) :: set_tilt, set_hvkicks, set_spin
  real(rp), optional :: s_pos, mat6(6,6), s_out, spin_qrot(0:3), time
  logical, optional :: make_matrix
end subroutine

subroutine offset_photon (ele, orbit, set, offset_position_only, rot_mat)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct), target :: orbit
  logical :: set
  logical, optional :: offset_position_only
  real(rp), optional :: rot_mat(3,3)
end subroutine

subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)
  import
  implicit none
  type (ele_struct) ele
  real(rp) phi_a
  real(rp) phi_b
  real(rp) mat4(4,4)
end subroutine

subroutine orbit_amplitude_calc (ele, orb, amp_a, amp_b, amp_na, amp_nb)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orb
  real(rp), optional :: amp_a, amp_b, amp_na, amp_nb
end subroutine

function orbit_to_floor_phase_space (orbit, ele) result (floor_phase_space)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) floor_phase_space(6)
end function

function orbit_to_local_curvilinear (orbit, ele, z_direction, relative_to) result (local_position)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (floor_position_struct) local_position
  integer, optional :: z_direction, relative_to
end function

function orbit_too_large (orbit, param, check_momentum) result (is_too_large)
  import
  implicit none
  type (coord_struct) orbit
  type (lat_param_struct), optional :: param
  logical, optional :: check_momentum
  logical is_too_large
end function

subroutine order_super_lord_slaves (lat, ix_lord)
  import
  implicit none
  type (lat_struct), target :: lat
  integer ix_lord
end subroutine

subroutine converter_distribution_parser (ele, delim, delim_found, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  character(*) delim
  logical delim_found, err_flag
end subroutine

function particle_is_moving_backwards (orbit) result (is_moving_backwards)
  import
  implicit none
  type (coord_struct) orbit
  logical is_moving_backwards
end function

function particle_is_moving_forward (orbit, dir) result (is_moving_forward)
  import
  implicit none
  type (coord_struct) orbit
  integer, optional :: dir
  logical is_moving_forward
end function

function particle_rf_time (orbit, ele, reference_active_edge, s_rel, time_coords, rf_freq, rf_clock_harmonic, abs_time) result (time)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  real(rp), optional :: s_rel, rf_freq
  real(rp) time
  integer, optional :: rf_clock_harmonic
  logical, optional :: reference_active_edge, time_coords, abs_time
end function

function patch_flips_propagation_direction (x_pitch, y_pitch) result (is_flip)
  import
  implicit none
  real(rp) x_pitch, y_pitch
  logical is_flip
end function patch_flips_propagation_direction

function patch_length (patch, ref_coords) result (length)
  import
  implicit none
  type (ele_struct) patch
  real(rp) length
  integer, optional :: ref_coords
end function

subroutine phase_space_fit (x, xp, twiss, tune, emit, x_0, xp_0, chi, tol)
  import
  implicit none
  type (twiss_struct) twiss
  real(rp), optional :: tol
  real(rp) x(:), xp(:)
  real(rp) tune, emit
  real(rp) x_0, xp_0, chi
end subroutine

function physical_ele_end (track_end, orbit, ele_orientation, return_stream_end) result (physical_end)
  import
  implicit none
  type (coord_struct) orbit
  integer track_end, ele_orientation, physical_end
  logical, optional :: return_stream_end
end function

subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib, do_unlink)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) a_ptr
  character(*) attrib_name
  logical err_flag
  logical do_allocation
  logical, optional :: err_print_flag, do_unlink
  integer, optional :: ix_attrib
end subroutine

subroutine pointer_to_ele_multipole (ele, a_pole, b_pole, ksl_pole, pole_type)
  import
  implicit none
  type (ele_struct), target :: ele
  real(rp), pointer :: a_pole(:), b_pole(:), ksl_pole(:)
  integer, optional :: pole_type
end subroutine

function pointer_to_fibre(ele) result (assoc_fibre)
  import
  implicit none
  type (ele_struct), target :: ele
  type (fibre), pointer :: assoc_fibre
end function

function pointer_to_field_ele(ele, ix_field_ele, dz_offset) result (field_ele)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: field_ele
  integer ix_field_ele
  real(rp), optional :: dz_offset
end function

function pointer_to_girder(ele, ix_slave_back) result (girder)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: girder
  integer, optional :: ix_slave_back
end function

subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, a_ptr, err_flag, err_print_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) :: a_ptr
  integer :: ix_attrib
  logical err_flag, do_allocation
  logical, optional :: err_print_flag
end subroutine

function pointer_to_lord (slave, ix_lord, control, ix_slave_back, lord_type, ix_control, ix_ic) result (lord_ptr)
  import
  implicit none
  type (ele_struct), target :: slave
  type (control_struct), pointer, optional :: control
  type (ele_struct), pointer :: lord_ptr
  integer ix_lord
  integer, optional :: ix_slave_back, lord_type, ix_control, ix_ic
end function

function pointer_to_multipass_lord (ele, ix_pass, super_lord) result (multi_lord)
  import
  implicit none
  type (ele_struct), target :: ele
  integer, optional :: ix_pass
  type (ele_struct), pointer :: multi_lord
  type (ele_struct), pointer, optional :: super_lord
end function

function pointer_to_next_ele (this_ele, offset, skip_beginning, follow_fork) result (next_ele)
  import
  implicit none
  type (ele_struct), target :: this_ele
  type (ele_struct), pointer :: next_ele
  integer, optional :: offset
  logical, optional :: skip_beginning, follow_fork
end function

function pointer_to_super_lord (slave, control, ix_slave_back, ix_control, ix_ic) result (lord_ptr)
  import
  implicit none
  type (ele_struct), target :: slave
  type (control_struct), pointer, optional :: control
  type (ele_struct), pointer :: lord_ptr
  integer, optional :: ix_slave_back, ix_control, ix_ic
end function

function pointer_to_wake_ele (ele, delta_s) result (wake_ele)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: wake_ele
  real(rp), optional :: delta_s
end function

subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                  ptr_array, err_flag, err_print_flag, eles, ix_attrib, do_unlink)
  import
  implicit none
  type (lat_struct), target :: lat
  type (all_pointer_struct), allocatable :: ptr_array(:)
  character(*) ele_name, attrib_name
  logical err_flag
  logical do_allocation
  logical, optional :: err_print_flag, do_unlink
  type (ele_pointer_struct), optional, allocatable :: eles(:)
  integer, optional :: ix_attrib
end subroutine

function polar_to_spinor (polar) result (spinor)
  import
  implicit none
  type (spin_polar_struct) polar
  complex(rp) :: spinor(2)
end function

function polar_to_vec (polar) result (vec)
  import
  implicit none
  type (spin_polar_struct) polar
  real(rp) vec(3)
end function

subroutine ptc_bookkeeper (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine ptc_linear_isf_calc (branch, ele_isf)
  import
  implicit none
  type (branch_struct), target :: branch
  type (linear_ele_isf_struct), allocatable, target :: ele_isf(:)
end subroutine

subroutine ptc_ran_seed_put (iseed)
  implicit none
  integer iseed
end subroutine

subroutine ptc_read_flat_file (flat_file, err_flag, lat, create_end_marker, from_mad)
  import
  implicit none
  type (lat_struct), optional, target :: lat
  character(*) flat_file(:)
  logical err_flag
  logical, optional :: create_end_marker, from_mad
end subroutine

subroutine ptc_set_rf_state_for_c_normal (nocavity)
  implicit none
  logical nocavity
end subroutine

subroutine ptc_spin_matching_calc (branch, match_info)
  import
  implicit none
  type (branch_struct), target :: branch
  type (spin_matching_struct), allocatable, target :: match_info(:)
end subroutine

subroutine ptc_transfer_map_with_spin (branch, t_map, s_map, orb0, err_flag, ix1, ix2, one_turn, unit_start)
  import
  implicit none
  type (branch_struct) :: branch
  type (taylor_struct) :: t_map(6), s_map(4)
  type (coord_struct) orb0
  integer, optional :: ix1, ix2
  logical err_flag
  logical, optional :: one_turn, unit_start
end subroutine

subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)
  import
  implicit none
  real(rp) k1, length, rel_p
  real(rp) mat2(:,:)
  real(rp), optional :: z_coef(3), dz_dpz_coef(3)
end subroutine

subroutine radiation_integrals (lat, orb, mode, ix_cache, ix_branch, rad_int_by_ele)
  import
  implicit none
  type (lat_struct), target :: lat
  type (rad_int_all_ele_struct), optional, target :: rad_int_by_ele
  type (coord_struct), target :: orb(0:)
  type (normal_modes_struct) mode
  integer, optional :: ix_cache, ix_branch
end subroutine

subroutine re_allocate_eles (eles, n, save_old, exact)
  import
  implicit none
  type (ele_pointer_struct), allocatable :: eles(:)
  integer n
  logical, optional :: save_old, exact
end subroutine

subroutine reallocate_control (lat, n)
  import
  implicit none
  type (lat_struct) lat
  integer, intent(in) :: n
end subroutine

subroutine reallocate_expression_stack (stack, n, exact)
  import
  implicit none
  type (expression_atom_struct), allocatable :: stack(:)
  integer n
  logical, optional :: exact
end subroutine

subroutine reference_energy_correction (ele, orbit, particle_at, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) :: ele
  type (coord_struct) :: orbit
  real(rp), optional :: mat6(6,6)
  integer particle_at
  logical, optional :: make_matrix
end subroutine

function rel_tracking_charge_to_mass (orbit, ref_species) result (rel_charge)
  import
  implicit none
  type (coord_struct) orbit
  integer ref_species
  real(rp) rel_charge
end function

subroutine remove_eles_from_lat (lat, check_sanity)
  import
  implicit none
  type (lat_struct), target :: lat
  logical, optional :: check_sanity
end subroutine

subroutine remove_lord_slave_link (lord, slave)
  import
  implicit none
  type (ele_struct), target :: lord, slave
end subroutine

subroutine read_digested_bmad_file (in_file_name, lat, version, err_flag, parser_calling, lat_files)
  import
  implicit none
  type (lat_struct), target, intent(inout) :: lat
  integer version
  character(*) in_file_name
  logical, optional :: err_flag, parser_calling
  character(*), optional, allocatable :: lat_files(:)
end subroutine

subroutine reallocate_beam (beam, n_bunch, n_particle, extend)
  import
  implicit none
  type (beam_struct) beam
  integer n_bunch
  integer, optional :: n_particle
  logical, optional :: extend
end subroutine

subroutine reallocate_bunch (bunch, n_particle, save)
  import
  implicit none
  type (bunch_struct) bunch
  integer n_particle
  logical, optional :: save
end subroutine

function relative_mode_flip (ele1, ele2)
  import
  implicit none
  logical relative_mode_flip
  type (ele_struct) ele1
  type (ele_struct) ele2
end function

subroutine reverse_lat (lat_in, lat_rev, track_antiparticle)
  import
  implicit none
  type (lat_struct), target :: lat_in, lat_rev
  logical, optional :: track_antiparticle
end subroutine

function rf_clock_setup (branch, n_rf_included, n_rf_excluded) result (ok)
  import
  implicit none
  type (branch_struct), target :: branch
  integer n_rf_included, n_rf_excluded
  logical ok
end function

subroutine rf_coupler_kick (ele, param, particle_at, phase, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  real(rp) phase
  real(rp), optional :: mat6(6,6)
  integer particle_at
  logical, optional :: make_matrix
end subroutine

function rf_is_on (branch, ix_ele1, ix_ele2) result (is_on)
  import
  implicit none
  type (branch_struct), target :: branch
  integer, optional :: ix_ele1, ix_ele2
  logical is_on
end function

function rf_ref_time_offset (ele) result (time)
  import
  implicit none
  type (ele_struct) ele
  real(rp) time
end function

subroutine rotate_for_curved_surface (ele, orbit, set, rot_mat)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  real(rp) rot_mat(3,3)
  logical set
end subroutine

subroutine rotate_spin (rot_vec, spin, qrot)
  import
  implicit none
  real(rp) :: spin(3), rot_vec(3)
  real(rp), optional :: qrot(0:3)
end subroutine

subroutine rotate_spin_a_step (orbit, field, ele, ds)
  import
  implicit none
  type (coord_struct) orbit
  type (em_field_struct) field
  type (ele_struct) ele
  real(rp) ds
end subroutine

subroutine rotate_spin_given_field (orbit, sign_z_vel, BL, EL, qrot)
  import
  implicit none
  type (coord_struct) orbit
  real(rp), optional :: BL(3), EL(3), qrot(0:3)
  integer sign_z_vel
end subroutine

function s_body_calc (orbit, ele) result (s_body)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) s_body
end function

subroutine s_calc (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

subroutine save_a_beam_step (ele, beam, bunch_tracks, s_body, is_time_coords)
  import
  implicit none
  type (ele_struct), target :: ele
  type (beam_struct) beam
  type (bunch_track_struct), optional, target :: bunch_tracks(:)
  real(rp), optional :: s_body
  logical, optional :: is_time_coords
end subroutine

subroutine save_a_bunch_step (ele, bunch, bunch_track, s_body, is_time_coords)
  import
  implicit none
  type (ele_struct), target :: ele
  type (bunch_struct) bunch
  type (bunch_track_struct), optional, target :: bunch_track
  real(rp), optional :: s_body
  logical, optional :: is_time_coords
end subroutine

subroutine save_a_step (track, ele, param, local_ref_frame, orb, s_rel, save_field, mat6, make_matrix, rf_time, strong_beam)
  import
  implicit none
  type (track_struct), target :: track
  type (ele_struct), target :: ele
  type (lat_param_struct), intent(in) :: param
  type (coord_struct) orb
  type (strong_beam_struct), optional :: strong_beam
  real(rp) s_rel
  real(rp), optional :: mat6(6,6), rf_time
  logical local_ref_frame
  logical, optional :: save_field, make_matrix
end subroutine

subroutine sbend_body_with_k1_map (ele, dg, k_1, param, n_step, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  integer n_step
  real(rp) dg, k_1
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag, set_lords, err_id)
  import
  implicit none
  type (ele_struct), target :: ele
  integer, optional :: err_id
  logical, optional :: err_print_flag, set_lords
  logical err_flag
  character(*) set_string
end subroutine

subroutine set_ele_name(ele, name)
  import
  implicit none
  type (ele_struct) ele
  character(*) name
end subroutine

subroutine set_ele_real_attribute (ele, attrib_name, value, err_flag, err_print_flag)
  import
  implicit none
  type (ele_struct) ele
  real(rp) value
  logical, optional :: err_print_flag
  logical err_flag
  character(*) attrib_name
end subroutine

recursive subroutine set_ele_status_stale (ele, status_group, set_slaves)
  import
  implicit none
  type (ele_struct), target :: ele
  integer status_group
  logical, optional :: set_slaves
end subroutine

function set_emit_from_beam_init (beam_init_in, ele, species, modes, err_flag) result (beam_init_set)
  import
  implicit none
  type (beam_init_struct), target :: beam_init_set, beam_init_in
  type (ele_struct) ele
  type (normal_modes_struct), optional :: modes
  integer species
  logical, optional :: err_flag
end function

subroutine set_fringe_on_off (fringe_at, ele_end, on_or_off) 
  import
  implicit none
  integer ele_end, on_or_off
  real(rp) fringe_at
end subroutine

recursive subroutine set_lords_status_stale (ele, stat_group, control_bookkeeping, flag)
  import
  implicit none
  type (ele_struct) ele
  integer stat_group
  logical, optional :: control_bookkeeping
  integer, optional :: flag
end subroutine

subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch, saved_values, attribute, set_val)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), optional :: orb(0:)
  real(rp), optional, allocatable :: saved_values(:)
  integer :: key, switch
  integer, optional :: ix_branch, set_val
  logical, optional :: use_ref_orb
  character(*), optional :: attribute
end subroutine

subroutine set_orbit_to_zero (orbit, n1, n2, ix_noset)
  import
  implicit none
  type (coord_struct) orbit(0:)
  integer n1, n2
  integer, optional :: ix_noset
end subroutine

subroutine set_particle_from_rf_time (rf_time, ele, reference_active_edge, orbit)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  real(rp) rf_time
  logical reference_active_edge
end subroutine

subroutine set_ptc (e_tot, particle, taylor_order, integ_order, n_step, no_cavity, force_init) 
  import
  implicit none
  integer, optional :: integ_order, particle, n_step, taylor_order
  real(rp), optional :: e_tot
  logical, optional :: no_cavity, force_init
end subroutine

subroutine set_ptc_base_state (component, set_val, old_val)
  implicit none
  character(*) component
  logical set_val
  logical, optional :: old_val
end subroutine

subroutine set_status_flags (bookkeeping_state, stat)
  import
  implicit none
  type (bookkeeping_state_struct) bookkeeping_state
  integer stat
end subroutine

subroutine set_twiss(branch, twiss_ele, ix_ele, match_deta_ds, err_flag, print_err)
  import
  implicit none
  type (branch_struct), target :: branch
  type (ele_struct) twiss_ele
  integer ix_ele
  logical match_deta_ds, err_flag
  logical, optional :: print_err
end subroutine

subroutine set_z_tune (branch, z_tune, ok, print_err)
  import
  implicit none
  type (branch_struct), target :: branch
  real(rp) :: z_tune
  logical, optional :: ok, print_err
end subroutine

subroutine set_on (key, lat, on_switch, orb)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), optional :: orb(0:)
  integer key
  logical on_switch
end subroutine

subroutine set_ele_defaults (ele, do_allocate)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: do_allocate
end subroutine

function set_tune (phi_a_set, phi_b_set, dk1, eles, branch, orb, print_err) result (ok)
  import
  implicit none
  type (branch_struct), target :: branch
  type (coord_struct), allocatable :: orb(:)
  type (ele_pointer_struct) eles(:)
  real(rp) phi_a_set
  real(rp) phi_b_set
  real(rp) dk1(:)
  logical, optional :: print_err
  logical ok
end function

function set_tune_via_group_knobs (phi_set, branch, group_knobs, orb, print_err) result (ok)
  import
  implicit none
  type (branch_struct), target :: branch
  type (coord_struct), allocatable :: orb(:)
  character(*) group_knobs(2)
  real(rp) phi_set(2)
  logical, optional :: print_err
  logical ok
end function

function significant_difference (value1, value2, abs_tol, rel_tol) result (is_different)
  import
  implicit none
  real(rp), intent(in) :: value1, value2
  real(rp), intent(in), optional :: abs_tol, rel_tol
  logical is_different
end function

subroutine slice_lattice (lat, ele_list, error, do_bookkeeping)
  import
  implicit none
  type (lat_struct), target :: lat
  character(*) ele_list
  logical error
  logical, optional :: do_bookkeeping
end subroutine

subroutine sol_quad_mat6_calc (ks, k1, tilt, s_len, ele, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orbit
  real(rp) ks, k1, tilt, s_len
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine solenoid_track_and_mat (ele, length, param, start_orb, end_orb, mat6, make_matrix)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) start_orb, end_orb
  real(rp) length
  real(rp), optional :: mat6(:,:)
  logical, optional :: make_matrix
end subroutine

subroutine spin_concat_linear_maps (err_flag, mat1, branch, n1, n2, mat1_ele, orbit, excite_zero)
  import
  implicit none
  type (spin_orbit_map1_struct) mat1
  type (spin_orbit_map1_struct), optional :: mat1_ele(0:)
  type (branch_struct), target :: branch
  type (coord_struct), optional :: orbit(0:)
  logical err_flag
  integer n1, n2
  character(*), optional :: excite_zero(3)
end subroutine

function spin_depolarization_rate (branch, match_info, rad_int_by_ele) result (depol_rate)
  import
  type (branch_struct), target :: branch
  type (spin_matching_struct) match_info(0:)
  type (rad_int_all_ele_struct) rad_int_by_ele
  real(rp) depol_rate
end function

function spin_dn_dpz_from_mat8 (mat_1turn, dn_dpz_partial, error) result (dn_dpz)
  import
  implicit none
  real(rp) mat_1turn(8,8), dn_dpz(3)
  real(rp), optional :: dn_dpz_partial(3,3)
  logical error
end function

function spin_dn_dpz_from_qmap (orb_mat, q_map, dn_dpz_partial, dn_dpz_partial2, error, n0) result (dn_dpz)
  import
  implicit none
  real(rp), optional :: n0(3)
  real(rp) orb_mat(6,6), q_map(0:3,0:6), dn_dpz(3)
  real(rp) :: dn_dpz_partial(3,3), dn_dpz_partial2(3,3)
  logical error
end function

subroutine spin_map1_normalize (spin1)
  import
  implicit none
  real(rp) spin1(0:3,0:6)
end subroutine

subroutine spin_mat_to_eigen (orb_mat, spin_map, eigen_val, orb_evec, n0, spin_evec, error)
  import
  implicit none
  real(rp) orb_mat(6,6), spin_map(0:3,0:6), n0(3)
  complex(rp) eigen_val(6), orb_evec(6,6), spin_evec(6,3)
  logical error
end subroutine

subroutine spin_mat8_resonance_strengths (orb_evec, mat8, xi_sum, xi_diff)
  import
  implicit none
  real(rp) mat8(6,6), xi_sum, xi_diff
  complex(rp) orb_evec(6)
end subroutine

function spin_omega (field, orbit, sign_z_vel, phase_space_coords) result (omega)
  import
  implicit none
  type (em_field_struct) :: field
  type (coord_struct) :: orbit
  integer sign_z_vel
  logical, optional :: phase_space_coords
  real(rp) omega(3)
end function

subroutine spin_quat_resonance_strengths (orb_evec, spin_q, xi_sum, xi_diff)
  import
  implicit none
  real(rp) spin_q(0:3,0:6), xi_sum, xi_diff
  complex(rp) orb_evec(6)
end subroutine

function spin_taylor_to_linear (spin_taylor, normalize, dref_orb, is_on) result (spin_map1)
  import
  implicit none
  type (taylor_struct), target :: spin_taylor(0:3)
  logical normalize, is_on
  real(rp) dref_orb(6), spin_map1(0:3,0:6)
end function

function spinor_to_polar (spinor) result (polar)
  import
  implicit none
  complex(rp) spinor(2)
  type (spin_polar_struct) ::  polar
end function

function spinor_to_vec (spinor) result (vec)
  import
  implicit none
  complex(rp) spinor(2)
  real(rp) vec(3)
end function

subroutine spline_fit_orbit (start_orb, end_orb, spline_x, spline_y)
  import
  implicit none
  type (coord_struct) start_orb, end_orb
  real(rp) spline_x(0:3), spline_y(0:3)
end subroutine

subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_sanity, &
                                                 save_null_drift, err_flag, choose_max, ix_insert)
  import
  implicit none
  type (lat_struct), target :: lat
  real(rp) s_split
  integer ix_branch
  integer ix_split
  integer, optional :: ix_insert
  logical split_done
  logical, optional :: add_suffix, check_sanity, save_null_drift, err_flag, choose_max
end subroutine

subroutine sprint_spin_taylor_map (ele, start_orbit)
  import
  implicit none
  type (ele_struct) ele
  real(rp), optional :: start_orbit(6)
end subroutine

subroutine start_branch_at (lat, ele_start, move_end_marker, error)
  import
  implicit none
  type (lat_struct), target :: lat
  character(*) ele_start
  logical move_end_marker, error
end subroutine

function stream_ele_end (physical_end, ele_orientation) result (stream_end)
  import
  implicit none
  integer stream_end, ele_orientation, physical_end
end function

subroutine strong_beam_sigma_calc (ele, s_pos, sigma, bbi_const, dsigma_ds)
  import
  implicit none
  type (ele_struct) ele
  real(rp) s_pos, bbi_const, sigma(2), dsigma_ds(2)
end subroutine

function strong_beam_strength (ele) result (strength)
  import
  implicit none
  type (ele_struct) ele
  real(rp) strength
end function

subroutine symp_lie_bmad (ele, param, orbit, track, mat6, make_matrix, offset_ele)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: orbit
  type (lat_param_struct)  param
  type (track_struct), optional, target :: track
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix, offset_ele
end subroutine

subroutine taper_mag_strengths (lat, ref_lat, except, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  type (lat_struct), optional, target :: ref_lat
  character(*), optional :: except
  logical, optional :: err_flag
end subroutine

subroutine tilt_coords (tilt_val, coord, mat6, make_matrix)
  import
  implicit none
  real(rp) tilt_val
  real(rp) coord(:)
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine tilt_coords_photon (tilt_val, coord, w_mat)
  import
  implicit none
  real(rp) tilt_val, coord(:)
  real(rp), optional :: w_mat(3,3)
end subroutine

subroutine tilt_mat6 (mat6, tilt)
  import
  implicit none
  real(rp) tilt, mat6(6,6)
end subroutine

subroutine to_surface_coords (lab_orbit, ele, surface_orbit)
  import
  implicit none
  type (coord_struct) lab_orbit, surface_orbit
  type (ele_struct) ele
end subroutine

subroutine track_a_beambeam (orbit, ele, param, track, mat6, make_matrix)
  import
  implicit none
  type (coord_struct), target :: orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (track_struct), optional :: track
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_bend (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_converter (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_crab_cavity (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_drift (orb, length, mat6, make_matrix, ele_orientation, include_ref_motion, time)
  import
  implicit none
  type (coord_struct) orb
  real(rp) length
  real(rp), optional :: mat6(6,6), time
  integer, optional :: ele_orientation
  logical, optional :: make_matrix, include_ref_motion
end subroutine

subroutine track_a_drift_photon (orb, length, phase_relative_to_ref)
  import
  implicit none
  type (coord_struct) orb
  real(rp) length
  logical phase_relative_to_ref
end subroutine

subroutine track_a_gkicker (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_mask (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_match (orbit, ele, param, err_flag, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix, err_flag
end subroutine

subroutine track_a_pickup (orbit, ele, param, err_flag, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix, err_flag
end subroutine

subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, track_spin, mat6, make_matrix)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) orbit
  real(rp), optional :: mat6(6,6), s_ent, ds_ref
  logical, optional :: drift_to_exit, track_spin, make_matrix
end subroutine

subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_sad_mult (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_sol_quad (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

recursive subroutine track_a_foil (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_wiggler (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_zero_length_element (orbit, ele, param, err_flag, track)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine track_all (lat, orbit, ix_branch, track_state, err_flag, orbit0)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), allocatable, target :: orbit(:)
  type (coord_struct), optional, allocatable, target :: orbit0(:)
  integer, optional :: ix_branch, track_state
  logical, optional :: err_flag
end subroutine

subroutine track_bunch_time (bunch, branch, t_end, s_end, dt_step, extra_field)
  import
  implicit none
  type (bunch_struct), target :: bunch
  type (branch_struct), target :: branch
  real(rp) t_end, s_end
  real(rp), optional :: dt_step(:)
  type (em_field_struct), optional :: extra_field(:)
end subroutine

subroutine track_from_s_to_s (lat, s_start, s_end, orbit_start, orbit_end, all_orb, ix_branch, track_state, ix_ele_end)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct) orbit_start, orbit_end
  type (coord_struct), optional, allocatable :: all_orb(:)
  real(rp) s_start, s_end
  integer, optional :: ix_branch, track_state, ix_ele_end
end subroutine

subroutine track_many (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct)  orbit(0:)
  integer ix_start
  integer ix_end
  integer direction
  integer, optional :: ix_branch, track_state
end subroutine

subroutine track_to_surface (ele, orbit, param, w_surface)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orbit
  type (lat_param_struct) param
  real(rp) :: w_surface(3,3)
end subroutine

recursive subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, &
                                                        ignore_radiation, make_map1, init_to_edge)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct)   :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical, optional :: err_flag, ignore_radiation, init_to_edge
  logical, optional :: make_map1
end subroutine

subroutine track1_bmad (orbit, ele, param, err_flag, track, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical, optional :: err_flag
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track1_bmad_photon (orbit, ele, param, err_flag)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
end subroutine

subroutine track1_bunch_space_charge (bunch, ele, err, track_to_same_s, bunch_track)
  import
  implicit none
  type (bunch_struct), target :: bunch
  type (ele_struct), target :: ele
  type (bunch_track_struct), optional :: bunch_track
  logical err
  logical, optional :: track_to_same_s
end subroutine

subroutine track1_linear (orbit, ele, param)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_runge_kutta (orbit, ele, param, err_flag, track, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track1_spin (start_orb, ele, param, end_orb, make_quaternion)
  import
  implicit none
  type (coord_struct) :: start_orb, end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: make_quaternion
end subroutine

subroutine track1_spin_integration (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (coord_struct) :: end_orb
end subroutine

subroutine track1_spin_magnus (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb, end_orb
  type (ele_struct) ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_spin_taylor (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb, end_orb
  type (ele_struct) ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_symp_lie_ptc (orbit, ele, param, track)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
end subroutine

subroutine track1_taylor (orbit, ele, taylor, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) :: orbit
  type (ele_struct), target :: ele
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
  type (taylor_struct), optional, target :: taylor(6)
end subroutine

subroutine track1_time_runge_kutta (orbit, ele, param, err_flag, track, t_end, dt_step)
  import
  implicit none
  real(rp), optional :: t_end, dt_step
  type (coord_struct) :: orbit
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine tracking_rad_map_setup (ele, tollerance, ref_edge, rad_map, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (rad_map_struct) rad_map
  real(rp) tollerance
  integer ref_edge
  logical err_flag
end subroutine

subroutine transfer_ac_kick (ac_kick_in, ac_kick_out)
  import
  implicit none
  type (ac_kicker_struct), pointer :: ac_kick_in, ac_kick_out
end subroutine transfer_ac_kick

subroutine transfer_branch (branch1, branch2)
  import
  implicit none
  type (branch_struct) :: branch1
  type (branch_struct) :: branch2
end subroutine

subroutine transfer_branch_parameters (branch_in, branch_out)
  import
  implicit none
  type (branch_struct), intent(in) :: branch_in
  type (branch_struct) :: branch_out
end subroutine

subroutine transfer_branches (branch1, branch2)
  import
  implicit none
  type (branch_struct) :: branch1(:)
  type (branch_struct) :: branch2(:)
end subroutine

subroutine transfer_ele (ele1, ele2, nullify_pointers)
  import
  implicit none
  type (ele_struct), target :: ele1
  type (ele_struct) :: ele2
  logical, optional :: nullify_pointers
end subroutine

subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)
  import
  implicit none
  type (ele_struct) ele_in, ele_out
  integer, optional :: taylor_order
end subroutine

subroutine transfer_eles (ele1, ele2)
  import
  implicit none
  type (ele_struct), intent(inout) :: ele1(:)
  type (ele_struct), intent(inout) :: ele2(:)
end subroutine

subroutine transfer_fieldmap (ele_in, ele_out, who)
  import
  implicit none
  type (ele_struct) :: ele_in, ele_out
  integer who
end subroutine

subroutine transfer_lat (lat1, lat2)
  import
  implicit none
  type (lat_struct), intent(in) :: lat1
  type (lat_struct), intent(out) :: lat2
end subroutine

subroutine transfer_lat_parameters (lat_in, lat_out)
  import
  implicit none
  type (lat_struct), intent(in) :: lat_in
  type (lat_struct) :: lat_out
end subroutine

subroutine transfer_map_calc (lat, t_map, err_flag, ix1, ix2, ref_orb, ix_branch, one_turn, unit_start, concat_if_possible)
  import
  implicit none
  type (lat_struct), target :: lat
  type (taylor_struct) :: t_map(:)
  type (coord_struct), optional :: ref_orb
  integer, intent(in), optional :: ix1, ix2, ix_branch
  logical err_flag
  logical, optional :: one_turn, unit_start, concat_if_possible
end subroutine

subroutine transfer_mat_from_twiss (ele1, ele2, orb1, orb2, m)
  import
  implicit none
  type (ele_struct) ele1, ele2
  real(rp) orb1(6), orb2(6)
  real(rp) m(6,6)
end subroutine

subroutine transfer_mat2_from_twiss (twiss1, twiss2, mat)
  import
  implicit none
  type (twiss_struct) twiss1, twiss2
  real(rp) mat(2,2)
end subroutine

subroutine transfer_matrix_calc (lat, xfer_mat, xfer_vec, ix1, ix2, ix_branch, one_turn)
  import
  implicit none
  type (lat_struct), target :: lat
  real(rp) :: xfer_mat(6,6)
  real(rp), optional :: xfer_vec(6)
  integer, optional :: ix1, ix2, ix_branch
  logical, optional :: one_turn
end subroutine

subroutine transfer_twiss (ele_in, ele_out, reverse)
  import
  implicit none
  type (ele_struct) ele_in, ele_out
  logical, optional :: reverse
end subroutine

subroutine transfer_wake (wake_in, wake_out)
  import
  implicit none
  type (wake_struct), pointer :: wake_in, wake_out
end subroutine

subroutine transfer_wall3d (wall3d_in, wall3d_out)
  import
  implicit none
  type (wall3d_struct), pointer :: wall3d_in(:), wall3d_out(:)
end subroutine

subroutine twiss_and_track_from_s_to_s (branch, orbit_start, s_end, orbit_end, &
                                          ele_start, ele_end, err, compute_floor_coords, compute_twiss)
  import
  implicit none
  type (coord_struct) :: orbit_start, orbit_end
  type (ele_struct), optional :: ele_start, ele_end
  type (branch_struct), target :: branch
  real(rp) s_end
  logical, optional, intent(inout) :: err
  logical, optional :: compute_floor_coords, compute_twiss
end subroutine

recursive subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, track_downstream_end, &
                 orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords, compute_twiss, reuse_ele_end)
  import
  implicit none
  type (coord_struct), optional :: orbit_start, orbit_end
  type (ele_struct), optional, target :: ele_start, ele_end
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  real(rp) l_start, l_end
  logical track_upstream_end, track_downstream_end
  logical, optional :: err, compute_floor_coords, reuse_ele_end, compute_twiss
end subroutine

recursive subroutine twiss_at_element (ele, start_ele, end_ele, average)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), optional :: start_ele
  type (ele_struct), optional :: end_ele
  type (ele_struct), optional :: average
end subroutine

subroutine twiss_at_start (lat, status, ix_branch, type_out)
  import
  implicit none
  type (lat_struct), target :: lat
  integer, optional, intent(in) :: ix_branch
  integer, optional, intent(out) :: status
  logical, optional :: type_out
end subroutine

subroutine twiss1_propagate (twiss1, mat2, ele_type, length, twiss2, err)
  import
  implicit none
  type (twiss_struct) twiss1, twiss2
  integer ele_type
  real(rp) mat2(2,2), length
  logical err
end subroutine

subroutine twiss_from_mat2 (mat_in, twiss, stat, type_out)
  import
  implicit none
  type (twiss_struct)  twiss
  real(rp) mat_in(:, :)
  integer stat
  logical type_out
end subroutine

subroutine twiss_from_mat6 (mat6, map0, ele, stable, growth_rate, status, type_out)
  import
  implicit none
  type (ele_struct) :: ele
  real(rp) :: mat6(:,:), map0(:)
  real(rp) :: growth_rate
  logical :: stable, type_out
  integer :: status
end subroutine

subroutine twiss_from_tracking (lat, ref_orb0, symp_err, err_flag, d_orb) 
  import
  implicit none
  type (lat_struct), intent(inout) :: lat
  type (coord_struct), intent(in) :: ref_orb0
  real(rp), intent(in), optional :: d_orb(:)   
  real(rp), intent(out) :: symp_err
  logical err_flag
end subroutine

subroutine twiss_propagate1 (ele1, ele2, err)
  import
  implicit none
  type (ele_struct), target :: ele1, ele2
  logical, optional :: err
end subroutine

subroutine twiss_propagate_all (lat, ix_branch, err_flag, ie_start, ie_end)
  import
  implicit none
  type (lat_struct), target :: lat
  integer, optional :: ix_branch, ie_start, ie_end
  logical, optional :: err_flag
end subroutine

subroutine twiss_to_1_turn_mat (twiss, phi, mat2)
  import
  implicit none
  type (twiss_struct) twiss
  real(rp) phi, mat2(2,2)
end subroutine

subroutine type_coord (coord)
  import
  implicit none
  type (coord_struct) coord
end subroutine

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, type_control, &
      type_wake, type_floor_coords, type_field, type_wall, type_rad_kick, lines, n_lines)
  import
  implicit none
  type (ele_struct), target :: ele
  integer, optional, intent(in) :: type_mat6
  integer, optional, intent(out) :: n_lines
  integer, optional, intent(in) :: twiss_out, type_field, type_control
  logical, optional, intent(in) :: type_taylor, type_floor_coords
  logical, optional, intent(in) :: type_zero_attrib, type_wake, type_wall, type_rad_kick
  character(*), optional, allocatable :: lines(:)
end subroutine

subroutine type_taylors (bmad_taylor, max_order, lines, n_lines, file_id, out_style, clean, out_var_suffix)
  import
  implicit none
  type (taylor_struct), intent(in), target :: bmad_taylor(:)
  integer, optional, intent(out) :: n_lines
  integer, optional :: max_order, file_id
  logical, optional :: clean
  character(*), optional :: out_style, out_var_suffix
  character(*), optional, allocatable :: lines(:)
end subroutine

subroutine type_twiss (ele, frequency_units, compact_format, lines, n_lines)
  import
  implicit none
  type (ele_struct) ele
  integer, optional :: frequency_units
  integer, optional :: n_lines
  character(*), optional :: lines(:)
  logical, optional :: compact_format
end subroutine

subroutine unlink_fieldmap (cartesian_map, cylindrical_map, gen_grad_map, grid_field)
  import
  implicit none
  type (cartesian_map_struct), pointer, optional :: cartesian_map(:)
  type (cylindrical_map_struct), pointer, optional :: cylindrical_map(:)
  type (gen_grad_map_struct), pointer, optional :: gen_grad_map(:)
  type (grid_field_struct), pointer, optional :: grid_field(:)
end subroutine

subroutine unlink_wall3d (wall3d)
  import
  implicit none
  type (wall3d_struct), pointer :: wall3d(:)
end subroutine

subroutine update_floor_angles (floor, floor0)
  import
  implicit none
  type(floor_position_struct) :: floor
  type(floor_position_struct), optional :: floor0
end subroutine update_floor_angles

subroutine update_fibre_from_ele (ele, survey_needed)
  import
  implicit none
  type (ele_struct), target :: ele
  logical survey_needed
end subroutine

function valid_field_calc (ele, field_calc) result (is_valid)
  import
  implicit none
  type (ele_struct) ele
  integer field_calc
  logical is_valid
end function

function valid_fringe_type (ele, fringe_type) result (is_valid)
  import
  implicit none
  type (ele_struct) ele
  integer fringe_type
  logical is_valid
end function

function valid_mat6_calc_method (ele, species, mat6_calc_method) result (is_valid)
  import
  implicit none
  type (ele_struct), target :: ele
  integer mat6_calc_method, species
  logical is_valid
end function

function valid_spin_tracking_method (ele, spin_tracking_method) result (is_valid)
  import
  implicit none
  type (ele_struct) ele
  integer spin_tracking_method
  logical is_valid
end function

function valid_tracking_method (ele, species, tracking_method) result (is_valid)
  import
  implicit none
  type (ele_struct), target :: ele
  integer tracking_method, species
  logical is_valid
end function

function value_of_attribute (ele, attrib_name, err_flag, err_print_flag, err_value) result (value)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) a_ptr
  real(rp) value
  real(rp), optional :: err_value
  character(*) attrib_name
  logical, optional :: err_print_flag, err_flag
end function

function vec_to_polar (vec, phase) result (polar)
  import
  implicit none
  real(rp) vec(3)
  real(rp), optional :: phase
  type (spin_polar_struct) :: polar
end function

function vec_to_spinor (vec, phase) result (spinor)
  import
  implicit none
  real(rp) vec(3)
  real(rp), optional :: phase
  complex(rp) :: spinor(2)
end function

function w_mat_for_bend_angle (angle, ref_tilt, r_vec) result (w_mat)
  import
  implicit none
  real(rp) angle, ref_tilt, w_mat(3,3), t_mat(3,3)
  real(rp), optional :: r_vec(3)
end function w_mat_for_bend_angle

function w_mat_for_x_pitch (x_pitch, return_inverse) result (w_mat)
  import
  implicit none
  real(rp) x_pitch, c_ang, s_ang
  real(rp) :: w_mat(3,3)
  logical, optional :: return_inverse
end function w_mat_for_x_pitch

function w_mat_for_y_pitch (y_pitch, return_inverse) result (w_mat)
  import
  implicit none
  real(rp) y_pitch, c_ang, s_ang
  real(rp) :: w_mat(3,3)
  logical, optional :: return_inverse
end function w_mat_for_y_pitch

function w_mat_for_tilt (tilt, return_inverse) result (w_mat)
  import
  implicit none
  real(rp) tilt, c_ang, s_ang
  real(rp) :: w_mat(3,3)
  logical, optional :: return_inverse
end function w_mat_for_tilt

subroutine write_beam_floor_positions (file_name, beam, ele, new_file)
  import
  implicit none
  type (beam_struct), target :: beam
  type (ele_struct) :: ele
  character(*) file_name
  logical, optional :: new_file
end subroutine

subroutine write_bmad_lattice_file (bmad_file, lat, err, output_form, orbit0)
  import
  implicit none
  character(*) bmad_file
  type (lat_struct), target :: lat
  logical, optional :: err
  integer, optional :: output_form
  type (coord_struct), optional :: orbit0
end subroutine

subroutine write_digested_bmad_file (digested_name, lat,  n_files, file_names, extra, err_flag)
  import
  implicit none
  type (lat_struct), target, intent(in) :: lat
  integer, intent(in), optional :: n_files
  character(*) digested_name
  character(*), optional :: file_names(:)
  type (extra_parsing_info_struct), optional :: extra
  logical, optional :: err_flag
end subroutine

subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, use_matrix_model, &
                                           include_apertures, dr12_drift_max, ix_branch, err)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), allocatable, optional :: ref_orbit(:)
  real(rp), optional :: dr12_drift_max
  integer, optional :: ix_branch
  character(*) out_type, out_file_name
  logical, optional :: use_matrix_model, include_apertures, err
end subroutine

subroutine write_lattice_in_mad_format (out_type, out_file_name, lat, ref_orbit, use_matrix_model, &
                                           include_apertures, dr12_drift_max, ix_branch, err)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), allocatable, optional :: ref_orbit(:)
  real(rp), optional :: dr12_drift_max
  integer, optional :: ix_branch
  character(*) out_type, out_file_name
  logical, optional :: use_matrix_model, include_apertures, err
end subroutine

subroutine write_lattice_in_elegant_format (out_file_name, lat, ref_orbit, use_matrix_model, &
                                           include_apertures, dr12_drift_max, ix_branch, err)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), allocatable, optional :: ref_orbit(:)
  real(rp), optional :: dr12_drift_max
  integer, optional :: ix_branch
  character(*) out_file_name
  logical, optional :: use_matrix_model, include_apertures, err
end subroutine

subroutine write_lattice_in_sad_format (out_file_name, lat, include_apertures, ix_branch, err)
  import
  implicit none
  type (lat_struct), target :: lat
  integer, optional :: ix_branch
  character(*) out_file_name
  logical, optional :: include_apertures, err
end subroutine

subroutine write_lattice_in_julia (julia_name, lat, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  character(*) :: julia_name
  logical, optional :: err_flag
end subroutine

subroutine xsif_parser (xsif_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) xsif_file
  character(*), optional :: use_line
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  end subroutine

subroutine zero_ele_kicks (ele)
  import
  implicit none
  type (ele_struct) ele
end subroutine

subroutine zero_ele_offsets (ele)
  import
  implicit none
  type (ele_struct) ele
end subroutine

end interface

! Custom and hook routines

! Hook and custom abstract definitions 

abstract interface

subroutine apply_element_edge_kick_hook_def (orb, fringe_info, track_ele, param, finished, mat6, make_matrix, rf_time)
  import
  implicit none
  type (coord_struct) orb
  type (fringe_field_info_struct) fringe_info
  type (ele_struct) track_ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6), rf_time
  logical, optional :: make_matrix
  logical finished
end subroutine

subroutine check_aperture_limit_custom_def (orb, ele, particle_at, param, err_flag)
  import
  implicit none
  type (coord_struct) :: orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  integer particle_at
  logical err_flag
end subroutine

function distance_to_aperture_custom_def (orbit, particle_at, ele, no_aperture_here) result (dist)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp) dist
  integer particle_at
  logical no_aperture_here
end function

subroutine ele_geometry_hook_def (floor0, ele, floor, finished, len_scale)
  import
  implicit none
  type (ele_struct) ele
  type (floor_position_struct) floor0, floor
  real(rp) len_scale
  logical finished
end subroutine

recursive subroutine em_field_custom_def (ele, param, s_rel, orbit, local_ref_frame, field, calc_dfield, err_flag, &
                                                  calc_potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct), intent(in) :: orbit
  type (ele_pointer_struct), allocatable, optional :: used_eles(:)
  real(rp), intent(in) :: s_rel
  real(rp), optional :: rf_time
  logical local_ref_frame
  type (em_field_struct) :: field
  logical, optional :: err_flag, grid_allow_s_out_of_bounds
  logical, optional :: calc_dfield, calc_potential, use_overlap
end subroutine

subroutine ele_to_fibre_hook_def (ele, ptc_fibre, use_offsets, err_flag)
  import
  implicit none
  type (ele_struct) ele
  type (fibre) ptc_fibre
  logical use_offsets, err_flag
end subroutine

recursive subroutine lat_make_mat6_hook_def (finished, lat, ix_ele, ref_orb, ix_branch, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), optional :: ref_orb(0:)
  integer ix_ele, ix_branch
  logical finished, err_flag
end subroutine

subroutine radiation_integrals_custom_def (lat, ir, orb, rad_int1, err_flag)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct) orb(0:)
  type (rad_int1_struct) rad_int1
  integer ir
  logical err_flag
end subroutine

subroutine init_custom_def (ele, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  logical err_flag
end subroutine

subroutine make_mat6_custom_def (ele, param, start_orb, end_orb, err_flag)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical err_flag, finished
end subroutine

subroutine time_runge_kutta_periodic_kick_hook_def (orbit, ele, param, stop_time, init_needed)
  import
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp) stop_time
  integer :: init_needed
end subroutine

subroutine track1_bunch_hook_def (bunch, ele, err, centroid, direction, finished, bunch_track)
  import
  implicit none
  type (bunch_struct), target :: bunch
  type (ele_struct), target :: ele
  type (coord_struct), optional :: centroid(0:)
  type (bunch_track_struct), optional :: bunch_track
  integer, optional :: direction
  logical err, finished
end subroutine

subroutine track1_custom_def (start_orb, ele, param, err_flag, finished, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical err_flag, finished, radiation_included
end subroutine

subroutine track_many_hook_def (finished, lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  import
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct)  orbit(0:)
  integer ix_start
  integer ix_end
  integer direction
  integer, optional :: ix_branch, track_state
  logical finished
end subroutine

subroutine track1_postprocess_def (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_preprocess_def (start_orb, ele, param, err_flag, finished, radiation_included, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct), target :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical err_flag, finished, radiation_included
end subroutine

subroutine track1_spin_custom_def (start_orb, ele, param, end_orb, err_flag, make_quaternion)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical err_flag
  logical, optional :: make_quaternion
end subroutine

subroutine track1_wake_hook_def (bunch, ele, finished)
  import
  implicit none
  type (bunch_struct) bunch
  type (ele_struct) ele
  logical finished
end subroutine

subroutine wall_hit_handler_custom_def (orb, ele, s)
  import
  implicit none
  type (coord_struct) :: orb
  type (ele_struct) :: ele
  real(rp) s
end subroutine

end interface

! Function pointers

procedure(apply_element_edge_kick_hook_def), pointer :: apply_element_edge_kick_hook_ptr => null()
procedure(check_aperture_limit_custom_def), pointer :: check_aperture_limit_custom_ptr => null()
procedure(distance_to_aperture_custom_def), pointer :: distance_to_aperture_custom_ptr => null()
procedure(ele_geometry_hook_def), pointer :: ele_geometry_hook_ptr => null()
procedure(em_field_custom_def), pointer :: em_field_custom_ptr => null()
procedure(ele_to_fibre_hook_def), pointer :: ele_to_fibre_hook_ptr => null()
procedure(radiation_integrals_custom_def), pointer :: radiation_integrals_custom_ptr => null()
procedure(init_custom_def), pointer :: init_custom_ptr => null()
procedure(lat_make_mat6_hook_def), pointer :: lat_make_mat6_hook_ptr => null()

procedure(make_mat6_custom_def), pointer :: make_mat6_custom_ptr => null()
procedure(time_runge_kutta_periodic_kick_hook_def), pointer :: time_runge_kutta_periodic_kick_hook_ptr => null()
procedure(track1_bunch_hook_def), pointer :: track1_bunch_hook_ptr => null()
procedure(track1_custom_def), pointer :: track1_custom_ptr => null()
procedure(track_many_hook_def), pointer :: track_many_hook_ptr => null()
procedure(track1_postprocess_def), pointer :: track1_postprocess_ptr => null()
procedure(track1_preprocess_def), pointer :: track1_preprocess_ptr => null()
procedure(track1_spin_custom_def), pointer :: track1_spin_custom_ptr => null()
procedure(track1_wake_hook_def), pointer :: track1_wake_hook_ptr => null()
procedure(wall_hit_handler_custom_def), pointer :: wall_hit_handler_custom_ptr => null()

contains

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

integer :: ix_ele, ib, ixe
integer, optional :: ix_branch

!

ele_ptr => null()

! ix_ele, ix_branch given

if (present(ix_branch)) then
  if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) return
  if (ix_ele < 0 .or. ix_ele > lat%branch(ix_branch)%n_ele_max) return
  ele_ptr => lat%branch(ix_branch)%ele(ix_ele)

! ix_ele = Nametable index

else
  if (ix_ele < 0) return
  ixe = ix_ele
  do ib = 0, ubound(lat%branch, 1)
    if (ixe > lat%branch(ib)%n_ele_max) then
      ixe = ixe - lat%branch(ib)%n_ele_max - 1
      cycle
    else
      ele_ptr => lat%branch(ib)%ele(ixe)
      return
    endif
  enddo
endif

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
! Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
!
! Function to return a pointer to an element in a lattice.
! This routine is overloaded by pointer_to_ele.
! See pointer_to_ele for more details.
!-

function pointer_to_ele3 (lat, ele_name) result (ele_ptr)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_ptr
type (ele_pointer_struct), allocatable :: eles(:)

integer n_loc
logical err

character(*) ele_name

!

ele_ptr => null()

call lat_ele_locator (ele_name, lat, eles, n_loc, err)
if (n_loc == 0) return

ele_ptr => eles(1)%ele

end function pointer_to_ele3

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)
!
! Function to return a pointer to an element in a lattice.
! This routine is overloaded by pointer_to_ele.
! See pointer_to_ele for more details.
!-

function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

type (lat_struct), target :: lat
type (ele_struct) foreign_ele
type (ele_struct), pointer :: ele_ptr

!

ele_ptr => lat%branch(foreign_ele%ix_branch)%ele(foreign_ele%ix_ele)

end function pointer_to_ele4

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function map1_times_map1 (map2, map1) result (map_out)
!
! Routine to concatenate two spin orbital linear maps.
!     map_out = map2(map1)
! Order is like applying matrices. map1 is before map2.
!
! Input:
!   map2      -- spin_orbit_map1_struct: Second map.
!   map1      -- spin_orbit_map1_struct: First map.
!
! Output:
!   map_out   -- spin_orbit_map1_struct: Concatenated map.
!-

function map1_times_map1 (map2, map1) result (map_out)

type (spin_orbit_map1_struct), intent(in) :: map1, map2
type (spin_orbit_map1_struct) map_out
real(rp) m2q1(0:3,6)
integer i

!

map_out%orb_mat = matmul(map2%orb_mat, map1%orb_mat)
map_out%vec0 = matmul(map2%orb_mat, map1%vec0) + map2%vec0

map_out%spin_q(:,0) = quat_mul(map2%spin_q(:,0), map1%spin_q(:,0))

m2q1 = matmul(map2%spin_q(:,1:6), map1%orb_mat)
do i = 1, 6
  map_out%spin_q(:,i) = quat_mul(map2%spin_q(:,0), map1%spin_q(:,i)) + quat_mul(m2q1(:,i), map1%spin_q(:,0))
enddo

end function map1_times_map1

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function em_field_plus_em_field (field1, field2) result (field_tot)
!
! Routine to add fields.
!
! Note: This subroutine is called by the overloaded plus sign:
!		field_tot = field1 + field2 
!
! Input:
!   field1 -- em_field_struct: Input field
!   field2 -- em_field_struct: Input field
!
! Output:
!   field_tot -- em_field_struct: Combined field.
!-

function em_field_plus_em_field (field1, field2) result (field_tot)

type (em_field_struct), intent(in) :: field1, field2
type (em_field_struct) field_tot

!

field_tot%e = field1%e + field2%e
field_tot%b = field1%b + field2%b

field_tot%de = field1%de + field2%de
field_tot%db = field1%db + field2%db

end function em_field_plus_em_field 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equal_ele (ele_out, ele_in)
!
! Subroutine that is used to set one element equal to another. 
! This routine takes care of the pointers in ele_out. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele_out = ele_in
!
! Input:
!   ele_in  -- Ele_struct: Input element.
!
! Output:
!   ele_out -- Ele_struct: Output element.
!-

subroutine ele_equal_ele (ele_out, ele_in)

implicit none
	
type (ele_struct), intent(inout), target :: ele_out
type (ele_struct), intent(in), target :: ele_in

call ele_equals_ele(ele_out, ele_in, .true.)

end subroutine ele_equal_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equals_ele (ele_out, ele_in, update_nametable)
!
! Subroutine that is used to set an element equal to another.
! Note: Use ele_equal_ele instead unless you know what you are doing.
! 
!
! Input:
!   ele_in            -- Ele_struct: Input element.
!   update_nametable  -- logical: If true, update the nametable. If false, do not.
!                         Note: nametable updates can take time if this routine
!                         is called a many times. See remove_eles_from_lat as an example.
!
! Output:
!   ele_out -- Ele_struct: Output element.
!-

subroutine ele_equals_ele (ele_out, ele_in, update_nametable)

implicit none
	
type (ele_struct), intent(inout), target :: ele_out
type (ele_struct), intent(in), target :: ele_in
type (ele_struct) ele_save
type (nametable_struct), pointer :: nt
type (converter_sub_distribution_struct), pointer :: sd_in, sd_out

integer i, j, n, n1, n2, lb(2), ub(2), ub1
logical update_nametable, comensurate

! 1) Save ele_out pointers in ele_save
! 2) Set ele_out = ele_in.

call transfer_ele (ele_out, ele_save)
call transfer_ele (ele_in, ele_out)

! ele_out%ix_ele and ele_out%ix_branch should not change.
! ele_out%branch should not change if ele_out is a component of a lat_struct.
!   Otherwise ele_out%lat should point to ele_in%lat (For cases where ele_out 
!   represents a sliced piece of ele_in)

ele_out%ix_ele    = ele_save%ix_ele    ! This should not change.
ele_out%ix_branch = ele_save%ix_branch ! This should not change.

if (update_nametable .and. ele_out%ix_ele > -1) then          ! If part of a lattice...
  ele_out%branch => ele_save%branch    !   then ele_out%branch should not change.
  if (associated(ele_out%branch)) then
    n = ele_nametable_index(ele_out)
    nt => ele_out%branch%lat%nametable
    ! During parsing the nametable may not have yet been updated so do not modify the
    ! nametable if this is the case.
    if (n <= nt%n_max .and. allocated(nt%name)) then
      if (nt%name(n) /= ele_out%name) call nametable_change1(ele_out%branch%lat%nametable, ele_out%name, n)
    endif
  endif
endif

! Transfer pointer info.
! When finished ele_out's pointers will be pointing to a different memory
! location from ele_in's so that the elements are separate.
! Exceptions: %em_field%mode%cylindrical_map, %em_field%mode%grid.

! %cartesian_map

ele_out%cartesian_map => ele_save%cartesian_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, cartesian_map$)

! %ac_kick

ele_out%ac_kick => ele_save%ac_kick  ! reinstate
call transfer_ac_kick (ele_in%ac_kick, ele_out%ac_kick)

! %cylindrical_map, etc.

ele_out%cylindrical_map => ele_save%cylindrical_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, cylindrical_map$)

ele_out%gen_grad_map => ele_save%gen_grad_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, gen_grad_map$)

ele_out%grid_field => ele_save%grid_field ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, grid_field$)

! %rad_map

if (associated(ele_in%rad_map)) then
  if (associated (ele_save%rad_map)) then
      ele_out%rad_map => ele_save%rad_map
  else
    allocate (ele_out%rad_map)
  endif
  ele_out%rad_map = ele_in%rad_map
else
  if (associated (ele_save%rad_map)) deallocate (ele_save%rad_map)
endif

! %r

if (associated(ele_in%r)) then
  if (associated (ele_save%r)) then
    if (all(lbound(ele_save%r) == lbound(ele_in%r)) .and. &
        all(ubound(ele_save%r) == ubound(ele_in%r)) ) then
      ele_out%r => ele_save%r
    else
      deallocate (ele_save%r)
      allocate (ele_out%r(lbound(ele_in%r,1):ubound(ele_in%r,1), &
                       lbound(ele_in%r,2):ubound(ele_in%r,2), &
                       lbound(ele_in%r,3):ubound(ele_in%r,3)))
    endif
  else
    allocate (ele_out%r(lbound(ele_in%r,1):ubound(ele_in%r,1), &
                     lbound(ele_in%r,2):ubound(ele_in%r,2), &
                     lbound(ele_in%r,3):ubound(ele_in%r,3)))
  endif
  ele_out%r = ele_in%r
else
  if (associated (ele_save%r)) deallocate (ele_save%r)
endif

! %photon

if (associated(ele_in%photon)) then
  ele_out%photon => ele_save%photon  ! reinstate
  if (.not. associated(ele_out%photon)) allocate(ele_out%photon)

  if (allocated (ele_in%photon%segmented%pt)) then
    lb = lbound(ele_in%photon%segmented%pt)
    ub = ubound(ele_in%photon%segmented%pt)
    if (allocated (ele_out%photon%segmented%pt)) then
      if (any(ub /= ubound(ele_out%photon%segmented%pt))) deallocate (ele_out%photon%segmented%pt)
    endif
    if (.not. allocated (ele_out%photon%segmented%pt)) allocate (ele_out%photon%segmented%pt(lb(1):ub(1), lb(2):ub(2)))
  else
    if (allocated(ele_out%photon%segmented%pt)) deallocate (ele_out%photon%segmented%pt)
  endif

  if (allocated (ele_in%photon%displacement%pt)) then
    lb = lbound(ele_in%photon%displacement%pt)
    ub = ubound(ele_in%photon%displacement%pt)
    if (allocated (ele_out%photon%displacement%pt)) then
      if (any(ub /= ubound(ele_out%photon%displacement%pt))) deallocate (ele_out%photon%displacement%pt)
    endif
    if (.not. allocated (ele_out%photon%displacement%pt)) allocate (ele_out%photon%displacement%pt(lb(1):ub(1), lb(2):ub(2)))
  else
    if (allocated(ele_out%photon%displacement%pt)) deallocate (ele_out%photon%displacement%pt)
  endif

  if (allocated (ele_in%photon%h_misalign%pt)) then
    lb = lbound(ele_in%photon%h_misalign%pt)
    ub = ubound(ele_in%photon%h_misalign%pt)
    if (allocated (ele_out%photon%h_misalign%pt)) then
      if (any(ub /= ubound(ele_out%photon%h_misalign%pt))) deallocate (ele_out%photon%h_misalign%pt)
    endif
    if (.not. allocated (ele_out%photon%h_misalign%pt)) allocate (ele_out%photon%h_misalign%pt(lb(1):ub(1), lb(2):ub(2)))
  else
    if (allocated(ele_out%photon%h_misalign%pt)) deallocate (ele_out%photon%h_misalign%pt)
  endif

  if (allocated (ele_in%photon%pixel%pt)) then
    lb = lbound(ele_in%photon%pixel%pt)
    ub = ubound(ele_in%photon%pixel%pt)
    if (allocated (ele_out%photon%pixel%pt)) then
      if (any(ub /= ubound(ele_out%photon%pixel%pt))) deallocate (ele_out%photon%pixel%pt)
    endif
    if (.not. allocated (ele_out%photon%pixel%pt)) allocate (ele_out%photon%pixel%pt(lb(1):ub(1), lb(2):ub(2)))
  else
    if (allocated(ele_out%photon%pixel%pt)) deallocate (ele_out%photon%pixel%pt)
  endif

  ele_out%photon = ele_in%photon
else
  if (associated (ele_save%photon)) deallocate (ele_save%photon)
endif

! %control
! Note: control%var will not be allocated for ramper slaves.

if (associated(ele_in%control)) then
  n1 = -1
  if (allocated(ele_in%control%var)) n1 = size(ele_in%control%var)
  n2 = -1
  if (allocated(ele_in%control%ramp)) n2 = size(ele_in%control%ramp)

  ele_out%control => ele_save%control   ! reinstate

  if (associated(ele_out%control)) then
    if (allocated(ele_out%control%var)) then
      if (size(ele_out%control%var) /= n1)  deallocate(ele_out%control%var)
    endif
    if (allocated(ele_out%control%ramp)) then
      if (size(ele_out%control%ramp) /= n2) deallocate(ele_out%control%ramp)
    endif
  endif

  if (.not. associated(ele_out%control)) allocate(ele_out%control)
  if (.not. allocated(ele_out%control%var) .and. n1 > -1) allocate(ele_out%control%var(n1))
  if (.not. allocated(ele_out%control%ramp) .and. n2 > -1) allocate(ele_out%control%ramp(n2))

  ele_out%control = ele_in%control

else
  if (associated (ele_save%control)) deallocate(ele_save%control)
endif

! %converter

if (associated(ele_in%converter)) then
  n = size(ele_in%converter%dist)
  ele_out%converter => ele_save%converter   ! reinstate
  if (associated(ele_out%converter)) then
    comensurate = .false.
    if (size(ele_out%converter%dist) /= n) goto 100
    do i = 1, n
      if (size(ele_in%converter%dist(i)%sub_dist) /= size(ele_out%converter%dist(i)%sub_dist)) goto 100
      do j = 1, size(ele_in%converter%dist(i)%sub_dist)
        sd_in => ele_in%converter%dist(i)%sub_dist(j)
        sd_out => ele_out%converter%dist(i)%sub_dist(j)
        if (.not. all(sd_in%prob_pc_r%prob == sd_out%prob_pc_r%prob)) goto 100
        if (size(sd_in%dir_out%beta%fit_1d_r) /= size(sd_out%dir_out%beta%fit_1d_r)) goto 100
        if (size(sd_in%dir_out%alpha_x%fit_1d_r) /= size(sd_out%dir_out%alpha_x%fit_1d_r)) goto 100
        if (size(sd_in%dir_out%alpha_y%fit_1d_r) /= size(sd_out%dir_out%alpha_y%fit_1d_r)) goto 100
      enddo
    enddo
    comensurate = .true.
    100 continue
    if (.not. comensurate) deallocate(ele_out%converter)
  endif

  if (.not. associated(ele_out%converter)) then
    allocate(ele_out%converter)
    allocate(ele_out%converter%dist(n))
  endif
  ele_out%converter = ele_in%converter

else
  if (associated (ele_save%converter)) deallocate(ele_save%converter)
endif

! %foil

if (associated(ele_in%foil)) then
  n = size(ele_in%foil%material)
  ele_out%foil => ele_save%foil   ! reinstate
  if (associated(ele_out%foil)) then
    comensurate = .false.
    if (size(ele_out%foil%material) /= n) deallocate(ele_out%foil)
  endif

  if (.not. associated(ele_out%foil)) then
    allocate(ele_out%foil)
    allocate(ele_out%foil%material(n))
  endif
  ele_out%foil = ele_in%foil

else
  if (associated (ele_save%foil)) deallocate(ele_save%foil)
endif

! %taylor

do i = 1, 6
  ele_out%taylor(i)%term => ele_save%taylor(i)%term ! reinstate
  ele_out%taylor(i) = ele_in%taylor(i)      ! use overloaded taylor_equal_taylor
enddo

! %spin_taylor

do i = 0, 3
  ele_out%spin_taylor(i)%term => ele_save%spin_taylor(i)%term ! reinstate
  ele_out%spin_taylor(i) = ele_in%spin_taylor(i)      ! use overloaded taylor_equal_taylor
enddo

! %wall3d

ele_out%wall3d => ele_save%wall3d        ! reinstate
call transfer_wall3d (ele_in%wall3d, ele_out%wall3d)

! %a_pole, and %b_pole

if (associated(ele_in%a_pole)) then
  ele_out%a_pole => ele_save%a_pole   ! reinstate
  ele_out%b_pole => ele_save%b_pole   ! reinstate
  call multipole_init (ele_out, magnetic$)
  ele_out%a_pole = ele_in%a_pole
  ele_out%b_pole = ele_in%b_pole
else
  if (associated (ele_save%a_pole)) deallocate (ele_save%a_pole, ele_save%b_pole)
endif

! %a_pole_elec, and %b_pole_elec

if (associated(ele_in%a_pole_elec)) then
  ele_out%a_pole_elec => ele_save%a_pole_elec   ! reinstate
  ele_out%b_pole_elec => ele_save%b_pole_elec   ! reinstate
  call multipole_init (ele_out, electric$)
  ele_out%a_pole_elec = ele_in%a_pole_elec
  ele_out%b_pole_elec = ele_in%b_pole_elec
else
  if (associated (ele_save%a_pole_elec)) deallocate (ele_save%a_pole_elec, ele_save%b_pole_elec)
endif

! %custom

if (associated(ele_in%custom)) then
  if (associated (ele_save%custom)) then
    if (size(ele_save%custom) == size(ele_in%custom)) then
      ele_out%custom => ele_save%custom
    else
      deallocate (ele_save%custom)
      allocate (ele_out%custom(size(ele_in%custom)))
    endif
  else
    allocate (ele_out%custom(size(ele_in%custom)))
  endif
  ele_out%custom = ele_in%custom
else
  if (associated (ele_save%custom)) deallocate (ele_save%custom)
endif

! %descrip

if (associated(ele_in%descrip)) then
  if (associated (ele_save%descrip)) then
    ele_out%descrip => ele_save%descrip
  else
    allocate (ele_out%descrip)
  endif
  ele_out%descrip = ele_in%descrip
else
  if (associated (ele_save%descrip)) deallocate (ele_save%descrip)
endif

! %mode3

if (associated(ele_in%mode3)) then
  if (associated (ele_save%mode3)) then
    ele_out%mode3 => ele_save%mode3
  else
    allocate (ele_out%mode3)
  endif
  ele_out%mode3 = ele_in%mode3
else
  if (associated (ele_save%mode3)) deallocate (ele_save%mode3)
endif

! %high_energy_space_charge

if (associated(ele_in%high_energy_space_charge)) then
  if (associated (ele_save%high_energy_space_charge)) then
    ele_out%high_energy_space_charge => ele_save%high_energy_space_charge
  else
    allocate (ele_out%high_energy_space_charge)
  endif
  ele_out%high_energy_space_charge = ele_in%high_energy_space_charge
else
  if (associated (ele_save%high_energy_space_charge)) deallocate (ele_save%high_energy_space_charge)
endif

! %wake

ele_out%wake => ele_save%wake  ! reinstate
call transfer_wake (ele_in%wake, ele_out%wake)

end subroutine ele_equals_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_vec_equal_ele_vec (ele1, ele2)
!
! Subroutine that is used to set one ele vector equal to another. 
! This routine takes care of the pointers in ele1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele1(:) = ele2(:)
!
! Input:
!   ele2(:) -- ele_struct: Input ele vector.
!
! Output:
!   ele1(:) -- ele_struct: Output ele vector.
!-

subroutine ele_vec_equal_ele_vec (ele1, ele2)

implicit none
	
type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(in) :: ele2(:)

integer i

! error check

if (size(ele1) /= size(ele2)) then
  print *, 'ERROR IN ele_vec_equal_ele_vec: ARRAY SIZES ARE NOT THE SAME!'
  if (global_com%exit_on_error) call err_exit
endif

! transfer

do i = 1, size(ele1)
  call ele_equal_ele (ele1(i), ele2(i))
enddo

end subroutine ele_vec_equal_ele_vec 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine lat_equal_lat (lat_out, lat_in)
!
! Subroutine that is used to set one lat equal to another. 
! This routine takes care of the pointers in lat_in. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		lat_out = lat_in
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat.
!-

subroutine lat_equal_lat (lat_out, lat_in)

implicit none

type (lat_struct), intent(inout), target :: lat_out
type (lat_struct), intent(in), target :: lat_in
type (branch_struct), pointer :: branch_out
type (ele_struct), pointer :: ele
type (control_struct), pointer :: c_in, c_out
integer i, n, nb, ne, ie, n_out, n_in
logical do_alloc

! Kill allociated PTC layouts if they exist

call kill_ptc_layouts(lat_out)

! If the element arrays have not been initialized in lat_in then deallocate lat_out.

if (.not. allocated (lat_in%branch) .or. .not. associated(lat_in%ele)) then
  call deallocate_lat_pointers (lat_out)
  return
endif

! Care must be taken here since lat%ele points to the same memory as lat%branch(0).
! First take care of the branch lines.

nb = ubound(lat_in%branch, 1)
call allocate_branch_array (lat_out, nb)

do i = 0, nb
  do_alloc = .true.
  if (associated(lat_out%branch(i)%ele)) then
    if (ubound(lat_out%branch(i)%ele, 1) >= lat_in%branch(i)%n_ele_max) do_alloc = .false.
  endif
  if (do_alloc) then
    ne = min(lat_in%branch(i)%n_ele_max+10, ubound(lat_in%branch(i)%ele, 1))
    call allocate_lat_ele_array (lat_out, ne, i)
  endif
  branch_out => lat_out%branch(i)
  branch_out = lat_in%branch(i)
  branch_out%lat => lat_out
  do ie = 0, ubound(branch_out%ele, 1)
    ele => branch_out%ele(ie)
    ele%ix_ele = ie
    ele%ix_branch = i
    ele%branch => branch_out
    if (associated(ele%ptc_fibre)) nullify(ele%ptc_fibre)
  enddo
enddo

lat_out%ele_init = lat_in%ele_init
nullify(lat_out%ele_init%branch)

! non-pointer transfer

call transfer_lat_parameters (lat_in, lat_out)

! super_ok$ is used to signal this routine to not do any ramper bookkeeping.
if (lat_in%ramper_slave_bookkeeping /= super_ok$) call ramper_slave_setup(lat_out, .true.)


end subroutine lat_equal_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine lat_vec_equal_lat_vec (lat1, lat2)
!
! Subroutine that is used to set one lat vector equal to another. 
! This routine takes care of the pointers in lat1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		lat1(:) = lat2(:)
!
! Input:
!   lat2(:) -- lat_struct: Input lat vector.
!
! Output:
!   lat1(:) -- lat_struct: Output lat vector.
!-

subroutine lat_vec_equal_lat_vec (lat1, lat2)

implicit none
	
type (lat_struct), intent(inout) :: lat1(:)
type (lat_struct), intent(in) :: lat2(:)

integer i

! error check

if (size(lat1) /= size(lat2)) then
  print *, 'ERROR IN lat_vec_equal_lat_vec: ARRAY SIZES ARE NOT THE SAME!'
  if (global_com%exit_on_error) call err_exit
endif

! transfer

do i = 1, size(lat1)
  call lat_equal_lat (lat1(i), lat2(i))
enddo

end subroutine lat_vec_equal_lat_vec 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine branch_equal_branch (branch1, branch2)
!
! Subroutine that is used to set one branch equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		branch1 = branch2
!
! Input:
!   branch2 -- branch_struct: Input branch.
!
! Output:
!   branch1 -- branch_struct: Output branch.
!-

subroutine branch_equal_branch (branch1, branch2)

implicit none
	
type (branch_struct), intent(inout) :: branch1
type (branch_struct), intent(in) :: branch2
integer i, n
logical do_alloc

!

do_alloc = .true.
if (associated(branch1%ele)) then
  if (ubound(branch1%ele, 1) >= branch2%n_ele_max) do_alloc = .false.
endif

if (do_alloc) then
  n = min(branch2%n_ele_max+10, ubound(branch2%ele, 1))
  call allocate_element_array (branch1%ele, n)
endif

do i = 0, branch2%n_ele_max
  branch1%ele(i)  = branch2%ele(i)
enddo

call transfer_wall3d (branch2%wall3d, branch1%wall3d)

call transfer_branch_parameters(branch2, branch1)

end subroutine branch_equal_branch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine coord_equal_coord (coord1, coord2)
!
! Subroutine that is used to set one coord equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		coord1 = coord2
!
! Input:
!   coord2 -- coord_struct: Input coord.
!
! Output:
!   coord1 -- coord_struct: Output coord.
!-

elemental subroutine coord_equal_coord (coord1, coord2)

implicit none
	
type (coord_struct), intent(inout) :: coord1
type (coord_struct), intent(in) :: coord2

!

coord1%vec = coord2%vec
coord1%spin = coord2%spin
 
end subroutine coord_equal_coord

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine taylor_equal_taylor (taylor1, taylor2)
!
! Subroutine that is used to set one taylor equal to another. 
! This routine takes care of the pointers in taylor1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		taylor1 = taylor2
!
! Input:
!   taylor2 -- Taylor_struct: Input taylor.
!
! Output:
!   taylor1 -- Taylor_struct: Output taylor.
!-

subroutine taylor_equal_taylor (taylor1, taylor2)

implicit none
	
type (taylor_struct), intent(inout) :: taylor1
type (taylor_struct), intent(in) :: taylor2

!

taylor1%ref = taylor2%ref

if (associated(taylor2%term)) then
  call init_taylor_series (taylor1, size(taylor2%term))
  taylor1%term = taylor2%term
  taylor1%ref = taylor2%ref
else
  if (associated (taylor1%term)) deallocate (taylor1%term)
endif

end subroutine taylor_equal_taylor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine taylors_equal_taylors (taylor1, taylor2)
!
! Subroutine to transfer the values from one taylor map to another:
!     Taylor1 <= Taylor2
!
! Input:
!   taylor2(:) -- Taylor_struct: Taylor map.
!
! Output:
!   taylor1(:) -- Taylor_struct: Taylor map. 
!-

subroutine taylors_equal_taylors (taylor1, taylor2)

implicit none

type (taylor_struct), intent(inout) :: taylor1(:)
type (taylor_struct), intent(in)    :: taylor2(:)

integer i

!

do i = 1, size(taylor1)
  taylor1(i) = taylor2(i)
enddo

end subroutine taylors_equal_taylors

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine em_taylor_equal_em_taylor (em_taylor1, em_taylor2)
!
! Subroutine that is used to set one em_taylor equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		em_taylor1 = em_taylor2
!
! Input:
!   em_taylor2 -- Em_taylor_struct: Input em_taylor.
!
! Output:
!   em_taylor1 -- Em_taylor_struct: Output em_taylor.
!-

subroutine em_taylor_equal_em_taylor (em_taylor1, em_taylor2)

implicit none
	
type (em_taylor_struct), intent(inout) :: em_taylor1
type (em_taylor_struct), intent(in) :: em_taylor2
integer n2

!

em_taylor1%ref = em_taylor2%ref

if (allocated(em_taylor2%term)) then
  n2 = size(em_taylor2%term)
  if (allocated(em_taylor1%term)) then
    if (size(em_taylor1%term) /= n2) then
      deallocate(em_taylor1%term)
      allocate (em_taylor1%term(n2))
    endif
  else
    allocate (em_taylor1%term(n2))
  endif
  em_taylor1%term = em_taylor2%term

else
  if (allocated(em_taylor1%term)) deallocate (em_taylor1%term)
endif

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine em_taylors_equal_em_taylors (em_taylor1, em_taylor2)
!
! Subroutine to transfer the values from one em_taylor map to another:
!     Em_taylor1 <= Em_taylor2
!
! Input:
!   em_taylor2(:) -- Em_taylor_struct: Em_taylor map.
!
! Output:
!   em_taylor1(:) -- Em_taylor_struct: Em_taylor map. 
!-

subroutine em_taylors_equal_em_taylors (em_taylor1, em_taylor2)

implicit none

type (em_taylor_struct), intent(inout) :: em_taylor1(:)
type (em_taylor_struct), intent(in)    :: em_taylor2(:)

integer i

!

do i = 1, size(em_taylor1)
  em_taylor1(i) = em_taylor2(i)
enddo

end subroutine em_taylors_equal_em_taylors

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_em_taylor_series (em_taylor, n_term, save_old)
!
! Subroutine to initialize a Bmad Em_taylor series (6 of these series make
! a Em_taylor map). Note: This routine does not zero the structure. The calling
! routine is responsible for setting all values.
!
! Input:
!   em_taylor   -- Em_taylor_struct: Old structure.
!   n_term      -- Integer: Number of terms to allocate. 
!                   n_term < 0 => em_taylor%term pointer will be disassociated.
!   save_old    -- Logical, optional: If True then save any old terms when
!                   em_taylor is resized. Default is False.
!
! Output:
!   em_taylor -- Em_taylor_struct: Initalized structure.
!-

subroutine init_em_taylor_series (em_taylor, n_term, save_old)

implicit none

type (em_taylor_struct) em_taylor
type (em_taylor_term_struct), allocatable :: term(:)
integer n_term
integer n
logical, optional :: save_old

!

if (n_term < 0) then
  if (allocated(em_taylor%term)) deallocate(em_taylor%term)
  return
endif

if (.not. allocated (em_taylor%term)) then
  allocate (em_taylor%term(n_term))
  return
endif

if (size(em_taylor%term) == n_term) return

!

if (logic_option (.false., save_old) .and. n_term > 0 .and. size(em_taylor%term) > 0) then
  n = min (n_term, size(em_taylor%term))
  call move_alloc(em_taylor%term, term)
  allocate (em_taylor%term(n_term))
  em_taylor%term(1:n) = term(1:n)
  deallocate (term)

else
  deallocate (em_taylor%term)
  allocate (em_taylor%term(n_term))
endif

end subroutine init_em_taylor_series

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine complex_taylor_equal_complex_taylor (complex_taylor1, complex_taylor2)
!
! Subroutine that is used to set one complex_taylor equal to another. 
! This routine takes care of the pointers in complex_taylor1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		complex_taylor1 = complex_taylor2
!
! Input:
!   complex_taylor2 -- complex_taylor_struct: Input complex_taylor.
!
! Output:
!   complex_taylor1 -- complex_taylor_struct: Output complex_taylor.
!-

subroutine complex_taylor_equal_complex_taylor (complex_taylor1, complex_taylor2)

implicit none
	
type (complex_taylor_struct), intent(inout) :: complex_taylor1
type (complex_taylor_struct), intent(in) :: complex_taylor2

!

complex_taylor1%ref = complex_taylor2%ref

if (associated(complex_taylor2%term)) then
  call init_complex_taylor_series (complex_taylor1, size(complex_taylor2%term))
  complex_taylor1%term = complex_taylor2%term
else
  if (associated (complex_taylor1%term)) deallocate (complex_taylor1%term)
endif

end subroutine complex_taylor_equal_complex_taylor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine complex_taylors_equal_complex_taylors (complex_taylor1, complex_taylor2)
!
! Subroutine to transfer the values from one complex_taylor map to another:
!     complex_taylor1 <= complex_taylor2
!
! Input:
!   complex_taylor2(:) -- complex_taylor_struct: complex_taylor map.
!
! Output:
!   complex_taylor1(:) -- complex_taylor_struct: complex_taylor map. 
!-

subroutine complex_taylors_equal_complex_taylors (complex_taylor1, complex_taylor2)

implicit none

type (complex_taylor_struct), intent(inout) :: complex_taylor1(:)
type (complex_taylor_struct), intent(in)    :: complex_taylor2(:)

integer i

!

do i = 1, size(complex_taylor1)
  complex_taylor1(i) = complex_taylor2(i)
enddo

end subroutine complex_taylors_equal_complex_taylors

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_complex_taylor_series (complex_taylor, n_term, save)
!
! Subroutine to initialize a Bmad complex_taylor series (6 of these series make
! a complex_taylor map). Note: This routine does not zero the structure. The calling
! routine is responsible for setting all values.
!
! Input:
!   complex_taylor -- complex_taylor_struct: Old structure.
!   n_term      -- Integer: Number of terms to allocate. 
!                   n_term < 1 => complex_taylor%term pointer will be disassociated.
!   save        -- Logical, optional: If True then save any old terms when
!                   complex_taylor is resized. Default is False.
!
! Output:
!   complex_taylor -- complex_taylor_struct: Initalized structure.
!-

subroutine init_complex_taylor_series (complex_taylor, n_term, save)

implicit none

type (complex_taylor_struct) complex_taylor
type (complex_taylor_term_struct), allocatable :: term(:)
integer n_term
integer n
logical, optional :: save

!

if (n_term < 1) then
  if (associated(complex_taylor%term)) deallocate(complex_taylor%term)
  return
endif

if (.not. associated (complex_taylor%term)) then
  allocate (complex_taylor%term(n_term))
  return
endif

if (size(complex_taylor%term) == n_term) return

!

if (logic_option (.false., save) .and. n_term > 0 .and. size(complex_taylor%term) > 0) then
  n = min (n_term, size(complex_taylor%term))
  allocate (term(n))
  term = complex_taylor%term(1:n)
  deallocate (complex_taylor%term)
  allocate (complex_taylor%term(n_term))
  complex_taylor%term(1:n) = term
  deallocate (term)

else
  deallocate (complex_taylor%term)
  allocate (complex_taylor%term(n_term))
endif

end subroutine init_complex_taylor_series

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine bunch_equal_bunch (bunch1, bunch2)
!
! Subroutine to set one particle bunch equal to another.
!
! Note: This subroutine is called by the overloaded equal sign:
!    bunch1 = bunch2
!
! Input: 
!   bunch2 -- bunch_struct: Input bunch
!
! Output
!   bunch1 -- bunch_struct: Output bunch
!-

subroutine bunch_equal_bunch (bunch1, bunch2)

implicit none

type (bunch_struct), intent(inout) :: bunch1
type (bunch_struct), intent(in)    :: bunch2

integer i, np

!

if (.not. allocated(bunch2%particle)) then
  if (allocated(bunch1%particle)) deallocate (bunch1%particle, bunch1%ix_z)
else
  np = size(bunch2%particle)
  if (.not. allocated(bunch1%particle)) allocate(bunch1%particle(np), bunch1%ix_z(np))
  if (size(bunch1%particle) /= size(bunch2%particle)) then
    deallocate (bunch1%particle, bunch1%ix_z)
    allocate (bunch1%particle(np), bunch1%ix_z(np))
  endif
  bunch1%particle = bunch2%particle
  bunch1%ix_z     = bunch2%ix_z
endif

bunch1%charge_tot  = bunch2%charge_tot
bunch1%charge_live = bunch2%charge_live
bunch1%t_center    = bunch2%t_center
bunch1%ix_ele      = bunch2%ix_ele
bunch1%ix_bunch    = bunch2%ix_bunch
bunch1%n_live      = bunch2%n_live

end subroutine bunch_equal_bunch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine beam_equal_beam (beam1, beam2)
!
! Subroutine to set one particle beam equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!    beam1 = beam2
!
! Input: 
!  beam2 -- beam_struct: Input beam
!
! Output
!  beam1 -- beam_struct: Output beam
!
!-

subroutine beam_equal_beam (beam1, beam2)

implicit none

type (beam_struct), intent(inout) :: beam1
type (beam_struct), intent(in)    :: beam2

integer i, j, n_bun, n_particle
logical allocate_this

! The following rule must be observed: If beam%bunch is allocated then
! beam%bunch%particle must be also.

n_bun = size(beam2%bunch)

allocate_this = .true.
if (allocated(beam1%bunch)) then
  if (size(beam1%bunch) /= size(beam2%bunch)) then
    do i = 1, size(beam1%bunch)
      deallocate (beam1%bunch(i)%particle)
    enddo
    deallocate (beam1%bunch)
  else
    allocate_this = .false.
  endif
endif

if (allocate_this) allocate (beam1%bunch(n_bun))

do i = 1, n_bun
  beam1%bunch(i) = beam2%bunch(i)
enddo

end subroutine beam_equal_beam

end module
