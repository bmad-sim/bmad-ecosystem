!+
! This file defines the interfaces for the BMAD subroutines
!-

module basic_bmad_interface

use bmad_struct

interface

subroutine add_lattice_control_structs (ele, n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field, add_at_end)
  import
  implicit none
  type (ele_struct) ele
  integer, optional :: n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field
  logical, optional :: add_at_end
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

subroutine apply_element_edge_kick (orb, fringe_info, track_ele, param, track_spin, mat6, make_matrix, rf_time)
  import
  implicit none
  type (coord_struct) orb
  type (fringe_edge_info_struct) fringe_info
  type (ele_struct) hard_ele, track_ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6), rf_time
  logical, optional :: make_matrix
  logical track_spin
end subroutine

subroutine attribute_bookkeeper (ele, param, force_bookkeeping)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  logical, optional :: force_bookkeeping
end subroutine

subroutine bbi_kick (x, y, r, kx, ky)
  import
  implicit none
  real(rp) x, y, r, kx, ky
end subroutine

subroutine bbi_kick_matrix (ele, param, orb, s_pos, mat6)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orb
  real(rp) s_pos
  real(rp) mat6(6,6)
end subroutine

subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) lat_file
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  character(*), optional :: use_line
end subroutine

subroutine bmad_parser2 (in_file, lat, orbit, make_mats6, err_flag)
  import
  implicit none
  character(*) in_file
  type (lat_struct), target :: lat
  type (coord_struct), optional :: orbit(0:)
  logical, optional :: make_mats6, err_flag
end subroutine

function branch_name(branch) result (name)
  import
  implicit none
  type (branch_struct), target :: branch
  character(40) name
end function

subroutine c_to_cbar (ele, cbar_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp) cbar_mat(2,2)
end subroutine

subroutine cbar_to_c (cbar_mat, a, b, c_mat)
  import
  implicit none
  real(rp) cbar_mat(2,2), c_mat(2,2)
  type (twiss_struct) a, b
end subroutine

subroutine calc_z_tune (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine crystal_attribute_bookkeeper (ele)
  import
  type (ele_struct) ele
end subroutine

subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
                       pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch)
  import
  implicit none
  type (lat_struct) lat
  type (lat_struct), optional, target :: low_E_lat, high_E_lat
  type (coord_struct), allocatable, optional, target :: low_E_orb(:), high_E_orb(:)
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

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable, target :: closed_orb(:)
  integer, optional :: direction, ix_branch, i_dim
  logical, optional, intent(out) :: err_flag
end subroutine

subroutine closed_orbit_from_tracking (lat, closed_orb, i_dim, &
                                     eps_rel, eps_abs, init_guess, err_flag)
  import
  type (lat_struct) lat
  type (coord_struct), allocatable :: closed_orb(:)
  type (coord_struct), optional :: init_guess
  real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)
  integer i_dim
  logical, optional :: err_flag
end subroutine

subroutine combine_consecutive_elements (lat)
  import
  type (lat_struct), target :: lat
end subroutine

subroutine create_field_overlap (lat, lord_name, slave_name, err_flag)
  import
  implicit none
  type (lat_struct) lat
  character(*) lord_name, slave_name
  logical err_flag
end subroutine

subroutine create_girder (lat, ix_ele, con, init_ele)
  import
  implicit none
  type (lat_struct) lat
  type (ele_struct), optional :: init_ele
  type (control_struct) con(:)
  integer, intent(in) :: ix_ele
end subroutine

subroutine create_group (lord, con, err, err_print_flag)
  import
  implicit none
  type (ele_struct) lord
  type (control_struct) con(:)
  logical err
  logical, optional :: err_print_flag
end subroutine

subroutine create_overlay (lord, contl, err, err_print_flag)
  import
  implicit none
  type (ele_struct) lord
  type (control_struct) contl(:)
  logical err
  logical, optional :: err_print_flag
end subroutine

subroutine create_uniform_element_slice (ele, param, i_slice, n_slice_tot, sliced_ele, s_start, s_end)
  import
  implicit none
  type (ele_struct) ele, sliced_ele
  type (lat_param_struct) param
  integer i_slice, n_slice_tot
  real(rp), optional :: s_start, s_end
end subroutine

subroutine create_unique_ele_names (lat, key, suffix)
  import
  type (lat_struct), target :: lat
  integer key
  character(*) suffix
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

function default_tracking_species (param) result (species)
  import
  implicit none
  type (lat_param_struct) param
  integer species
end function

subroutine do_mode_flip (ele, err_flag)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: err_flag
end subroutine

subroutine ele_compute_ref_energy_and_time (ele0, ele, param, err_flag)
  import
  type (ele_struct) ele0, ele
  type (lat_param_struct) param
  real(rp) e_tot_start, p0c_start, ref_time_start
  logical err_flag
end subroutine

subroutine exact_bend_multipole_field (ele, param, orbit, local_ref_frame, field, potential, calc_dfield)
  import
  implicit none
  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  type (em_field_struct) field
  type (em_potential_struct) potential
  logical local_ref_frame
  logical, optional :: calc_dfield
end subroutine

subroutine fibre_to_ele (ptc_fibre, branch, ix_ele, err_flag, from_mad)
  import
  implicit none
  type (fibre) ptc_fibre
  type (branch_struct) branch
  integer ix_ele
  logical err_flag
  logical, optional :: from_mad
end subroutine

subroutine find_element_ends (ele, ele1, ele2, ix_multipass)
  import
  implicit none
  type (ele_struct) ele
  type (ele_struct), pointer :: ele1, ele2
  integer, optional :: ix_multipass
end subroutine

subroutine find_matching_fieldmap (file_name, ele, t_type, match_ele, ix_field)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: match_ele
  integer t_type, ix_field
  character(*) file_name
end subroutine

subroutine init_a_photon_from_a_photon_init_ele (ele, param, orbit)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
end subroutine

subroutine init_exact_bend_multipole_coefs (ele, param, local_ref_frame, has_nonzero_pole)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  logical local_ref_frame
  logical has_nonzero_pole
end subroutine

subroutine init_wake (wake, n_sr_long, n_sr_trans, n_lr_mode, n_lr_spline, always_allocate)
  import
  implicit none
  type (wake_struct), pointer :: wake
  integer n_sr_long, n_sr_trans, n_lr_mode, n_lr_spline
  logical, optional :: always_allocate
end subroutine

subroutine insert_element (lat, insert_ele, insert_index, ix_branch, orbit)
  import
  implicit none
  type (lat_struct) lat
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
  type (lat_struct) lat
end subroutine

subroutine lat_compute_ref_energy_and_time (lat, err_flag)
  import
  type (lat_struct) lat
  logical err_flag
end subroutine

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

subroutine make_g_mats (ele, g_mat, g_inv_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp) g_mat(4,4)
  real(rp) g_inv_mat(4,4)
end subroutine

subroutine make_hybrid_lat (r_in, r_out, use_taylor, orb0)
  import
  implicit none
  type (lat_struct), target :: r_in
  type (lat_struct), target :: r_out
  logical, optional :: use_taylor
  type (coord_array_struct), optional :: orb0(0:)
end subroutine

recursive subroutine make_mat6 (ele, param, start_orb, end_orb, end_in, err_flag)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct), optional :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: end_in
  logical, optional :: err_flag
end subroutine

subroutine make_mat6_taylor (ele, param, start_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb
  type (lat_param_struct) param
end subroutine

subroutine make_mat6_bmad (ele, param, start_orb, end_orb, end_in, err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: end_in
  logical, optional :: err
end subroutine

subroutine make_mat6_bmad_photon (ele, param, start_orb, end_orb, end_in, err)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
  logical, optional :: end_in
  logical, optional :: err
end subroutine

subroutine make_mat6_runge_kutta (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_mat6_symp_lie_ptc (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_mat6_tracking (ele, param, start_orb, end_orb)
  import
  implicit none
  type (ele_struct), target :: ele
  type (coord_struct) :: start_orb, end_orb
  type (lat_param_struct) param
end subroutine

subroutine make_v_mats (ele, v_mat, v_inv_mat)
  import
  implicit none
  type (ele_struct) ele
  real(rp), optional :: v_mat(4,4)
  real(rp), optional :: v_inv_mat(4,4)
end subroutine

subroutine mat6_add_offsets (ele, param)
  import
  type (ele_struct) ele
  type (lat_param_struct) param
end subroutine

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

subroutine multipole_spin_precession (ele, param, orbit)
  import
  implicit none
  type (ele_struct) :: ele
  type (lat_param_struct) param
  type (coord_struct) orbit
end subroutine

subroutine name_to_list (lat, ele_names)
  import
  implicit none
  type (lat_struct) lat
  character(*) ele_names(:)
end subroutine

subroutine new_control (lat, ix_ele)
  import
  implicit none
  type (lat_struct) lat
  integer ix_ele
end subroutine

subroutine offset_particle (ele, param, set, coord, set_tilt, set_multipoles, set_hvkicks, set_z_offset, ds_pos, set_spin)
  import
  implicit none
  type (ele_struct) :: ele
  type (lat_param_struct) param
  type (coord_struct), intent(inout) :: coord
  integer particle
  logical, intent(in) :: set
  logical, optional, intent(in) :: set_multipoles, set_tilt, set_hvkicks, set_z_offset, set_spin
  real(rp), optional, intent(in) :: ds_pos
end subroutine

subroutine offset_photon (ele, coord, set, offset_position_only, rot_mat)
  import
  implicit none
  type (ele_struct) :: ele
  type (coord_struct) :: coord
  logical :: set
  logical, optional :: offset_position_only
  real(rp), optional :: rot_mat(3,3)
end subroutine

subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)
  import
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

subroutine order_super_lord_slaves (lat, ix_lord)
  import
  implicit none
  type (lat_struct), target :: lat
  integer ix_lord
end subroutine

function particle_rf_time (orbit, ele, apply_hard_edge_offset, s_rel) result (time)
  import
  type (coord_struct) orbit
  type (ele_struct) ele
  real(rp), optional :: s_rel
  real(rp) time
  logical apply_hard_edge_offset
end function

subroutine ptc_bookkeeper (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine ptc_read_flat_file (flat_file, err_flag, lat, create_end_marker, from_mad)
  import
  implicit none
  type (lat_struct), optional :: lat
  character(*) flat_file(:)
  logical err_flag
  logical, optional :: create_end_marker, from_mad
end subroutine

subroutine phase_space_fit (x, xp, twiss, tune, emit, x_0, xp_0, chi, tol)
  import
  implicit none
  type (twiss_struct) twiss
  real(rp), optional :: tol
  real(rp) x(:), xp(:)
  real(rp) tune, emit
  real(rp) x_0, xp_0, chi
end subroutine

subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) a_ptr
  character(*) attrib_name
  logical err_flag
  logical do_allocation
  logical, optional :: err_print_flag
  integer, optional :: ix_attrib
end subroutine

function pointer_to_lord (slave, ix_lord, control, ix_slave, field_overlap_ptr) result (lord_ptr)
  import
  implicit none
  type (ele_struct), target :: slave
  type (control_struct), pointer, optional :: control
  type (ele_struct), pointer :: lord_ptr
  integer, optional :: ix_slave
  integer ix_lord
  logical, optional :: field_overlap_ptr
end function

function pointer_to_slave (lord, ix_slave, control, field_overlap_ptr) result (slave_ptr)
  import
  implicit none
  type (ele_struct), target :: lord
  type (control_struct), pointer, optional :: control
  type (ele_struct), pointer :: slave_ptr
  integer ix_slave
  logical, optional :: field_overlap_ptr
end function

function pointer_to_wake_ele (ele, delta_s) result (wake_ele)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: wake_ele
  real(rp), optional :: delta_s
end function

subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                  ptr_array, err_flag, err_print_flag, eles, ix_attrib)
  import
  implicit none
  type (lat_struct) lat
  type (all_pointer_struct), allocatable :: ptr_array(:)
  character(*) ele_name, attrib_name
  logical err_flag
  logical do_allocation
  logical, optional :: err_print_flag
  type (ele_pointer_struct), optional, allocatable :: eles(:)
  integer, optional :: ix_attrib
end subroutine

subroutine quad_beta_ave (ele, beta_a_ave, beta_b_ave)
  import
  implicit none
  type (ele_struct) ele
  real(rp) beta_a_ave
  real(rp) beta_b_ave
end subroutine

subroutine radiation_integrals (lat, orb, mode, ix_cache, ix_branch, rad_int_by_ele)
  import
  implicit none
  type (lat_struct), target :: lat
  type (rad_int_all_ele_struct), optional :: rad_int_by_ele
  type (coord_struct), target :: orb(0:)
  type (normal_modes_struct) mode
  integer, optional :: ix_cache, ix_branch
end subroutine

subroutine reference_energy_correction (ele, orbit, particle_at)
  import
  implicit none
  type (ele_struct) :: ele
  type (coord_struct) :: orbit
  integer particle_at
end subroutine

subroutine remove_eles_from_lat (lat, check_sanity)
  import
  implicit none
  type (lat_struct) lat
  logical, optional :: check_sanity
end subroutine

subroutine read_digested_bmad_file (in_file_name, lat, version, err_flag)
  import
  implicit none
  type (lat_struct), target, intent(inout) :: lat
  integer version
  character(*) in_file_name
  logical, optional :: err_flag
end subroutine

subroutine reallocate_beam (beam, n_bunch, n_particle)
  import
  type (beam_struct) beam
  integer n_bunch, n_particle
end subroutine

subroutine reallocate_bunch (bunch, n_particle)
  import
  type (bunch_struct) bunch
  integer n_particle
end subroutine

function relative_mode_flip (ele1, ele2)
  import
  implicit none
  logical relative_mode_flip
  type (ele_struct) ele1
  type (ele_struct) ele2
end function

subroutine s_calc (lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine set_particle_from_rf_time (rf_time, ele, apply_hard_edge_offset, orbit)
  import
  implicit none
  type (ele_struct) ele
  type (coord_struct) orbit
  real(rp) rf_time
  logical apply_hard_edge_offset
end subroutine

subroutine set_status_flags (bookkeeping_state, stat)
  import
  implicit none
  type (bookkeeping_state_struct) bookkeeping_state
  integer stat
end subroutine

subroutine save_bunch_track (bunch, ele, s_travel)
  import
  implicit none
  type (bunch_struct) bunch
  type (ele_struct) ele
  real(rp) s_travel
end subroutine

subroutine sbend_body_with_k1_map (ele, param, n_step, orbit, mat6, make_matrix)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) orbit
  integer n_step
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine set_on (key, lat, on_switch, orb)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), optional :: orb(0:)
  integer key
  logical on_switch
end subroutine

subroutine set_ele_attribute (ele, set_string, lat, err_flag, err_print_flag)
  import
  implicit none
  type (ele_struct) ele
  type (lat_struct) lat
  character(*) set_string
  logical err_flag
  logical, optional :: err_print_flag
end subroutine

subroutine set_ele_defaults (ele, do_allocate)
  import
  implicit none
  type (ele_struct) ele
  logical, optional :: do_allocate
end subroutine

subroutine set_tune (phi_a_set, phi_b_set, dk1, lat, orb, ok)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)
  real(rp) phi_a_set
  real(rp) phi_b_set
  real(rp) dk1(:)
  logical ok
end subroutine

subroutine solenoid_track_and_mat (ele, param, start_orb, end_orb, mat6)
  import
  implicit none
  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) start_orb, end_orb
  real(rp), optional :: mat6(:,:)
end subroutine

subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_sanity, save_null_drift, err_flag)
  import
  implicit none
  type (lat_struct), target :: lat
  real(rp) s_split
  integer ix_branch
  integer ix_split
  logical split_done
  logical, optional :: add_suffix, check_sanity, save_null_drift, err_flag
end subroutine

subroutine track_a_drift_photon (orb, length, phase_relative_to_ref)
  import
  implicit none
  type (coord_struct) orb
  real(rp) length
  logical phase_relative_to_ref
end subroutine

subroutine transfer_matrix_calc (lat, xfer_mat, xfer_vec, ix1, ix2, ix_branch, one_turn)
  import
  implicit none
  type (lat_struct) lat
  real(rp) :: xfer_mat(:,:)
  real(rp), optional :: xfer_vec(:)
  integer, optional :: ix1, ix2, ix_branch
  logical, optional :: one_turn
end subroutine

subroutine transfer_map_calc (lat, t_map, err_flag, ix1, ix2, ix_branch, integrate, one_turn, unit_start)
  import
  implicit none
  type (lat_struct) lat
  type (taylor_struct) :: t_map(:)
  integer, intent(in), optional :: ix1, ix2, ix_branch
  logical err_flag
  logical, optional :: integrate, one_turn, unit_start
end subroutine

subroutine transfer_wake (wake_in, wake_out)
  import
  implicit none
  type (wake_struct), pointer :: wake_in, wake_out
end subroutine

subroutine tilt_coords (tilt_val, coord)
  import
  implicit none
  real(rp) tilt_val
  real(rp) coord(:)
end subroutine

subroutine track_all (lat, orbit, ix_branch, track_state, err_flag, orbit0)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct), allocatable, target :: orbit(:)
  type (coord_struct), optional, allocatable, target :: orbit0(:)
  integer, optional :: ix_branch, track_state
  logical, optional :: err_flag
end subroutine

subroutine track_from_s_to_s (lat, s_start, s_end, orbit_start, orbit_end, all_orb, ix_branch, track_state)
  import
  implicit none
  type (lat_struct) lat
  type (coord_struct) orbit_start, orbit_end
  type (coord_struct), optional, allocatable :: all_orb(:)
  real(rp) s_start, s_end
  integer, optional :: ix_branch, track_state
end subroutine

subroutine track_many (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  import
  implicit none
  type (lat_struct)  lat
  type (coord_struct)  orbit(0:)
  integer ix_start
  integer ix_end
  integer direction
  integer, optional :: ix_branch, track_state
end subroutine

subroutine track_backwards_time (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  import
  implicit none
  type (lat_struct)  lat
  type (coord_struct)  orbit(0:)
  integer ix_start
  integer ix_end
  integer direction
  integer, optional :: ix_branch, track_state
end subroutine

subroutine track_a_drift (orb, length, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orb
  real(rp) length
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_bend (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)
  import
  implicit none
  type (coord_struct) orbit
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp), optional :: mat6(6,6)
  logical, optional :: make_matrix
end subroutine

recursive subroutine track1 (start_orb, ele, param, end_orb, track, err_flag, ignore_radiation)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct)   :: ele
  type (lat_param_struct) :: param
  type (track_struct), optional :: track
  logical, optional :: err_flag, ignore_radiation
end subroutine

subroutine track1_backwards_time (end_orb, ele, param, start_orb, err_flag)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
end subroutine

subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
end subroutine

subroutine track1_bmad_photon (start_orb, ele, param, end_orb, err_flag)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  logical, optional :: err_flag
end subroutine

subroutine track1_linear (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine track1_symp_lie_ptc (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_symp_map (start_orb, ele, param, end_orb)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
end subroutine

subroutine track1_taylor (start_orb, ele, param, end_orb, taylor_in)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct) :: ele
  type (lat_param_struct) :: param
  type (taylor_struct), optional, target :: taylor_in(6)
end subroutine

subroutine track1_time_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)
  import
  implicit none
  type (coord_struct) :: start_orb
  type (coord_struct) :: end_orb
  type (ele_struct), target :: ele
  type (lat_param_struct), target :: param
  logical err_flag
  type (track_struct), optional :: track
end subroutine

subroutine transfer_mat_from_twiss (ele1, ele2, orb1, orb2, m)
  import
  implicit none
  type (ele_struct) ele1, ele2
  real(rp) orb1(6), orb2(6)
  real(rp) m(6,6)
end subroutine

subroutine twiss_and_track_from_s_to_s (branch, orbit_start, s_end, orbit_end, &
                                                               ele_start, ele_end, err, compute_floor_coords)
  import
  implicit none
  type (coord_struct) :: orbit_start, orbit_end
  type (ele_struct), optional :: ele_start, ele_end
  type (branch_struct) branch
  real(rp) s_end
  logical, optional, intent(inout) :: err
  logical, optional :: compute_floor_coords
end subroutine

subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, track_downstream_end, &
                                         orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords)
  import
  implicit none
  type (coord_struct), optional :: orbit_start, orbit_end
  type (ele_struct), optional :: ele_start, ele_end
  type (ele_struct) ele
  type (lat_param_struct) param
  real(rp) l_start, l_end
  logical track_upstream_end, track_downstream_end
  logical, optional :: err, compute_floor_coords
end subroutine

recursive subroutine twiss_at_element (ele, start_ele, end_ele, average)
  import
  implicit none
  type (ele_struct), target :: ele
  type (ele_struct), optional :: start_ele
  type (ele_struct), optional :: end_ele
  type (ele_struct), optional :: average
  integer ix_ele
end subroutine

subroutine twiss_at_start (lat, status, ix_branch)
  import
  implicit none
  type (lat_struct) lat
  integer, optional, intent(in) :: ix_branch
  integer, optional, intent(out) :: status
end subroutine

subroutine twiss1_propagate (twiss1, mat2, ele_type, length, twiss2, err)
  import
  implicit none
  type (twiss_struct) twiss1, twiss2
  integer ele_type
  real(rp) mat2(2,2), length
  logical err
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
  type (lat_struct), intent(inout) :: lat
  type (coord_struct), intent(in) :: ref_orb0
  real(rp), intent(in), optional :: d_orb(:)   
  real(rp), intent(out) :: symp_err
  logical err_flag
end subroutine

subroutine twiss_propagate1 (ele1, ele2, err)
  import
  implicit none
  type (ele_struct) ele1
  type (ele_struct) ele2
  logical, optional :: err
end subroutine

subroutine twiss_propagate_all (lat, ix_branch, err_flag, ie_start, ie_end, zero_uncalculated)
  import
  implicit none
  type (lat_struct) lat
  integer, optional :: ix_branch, ie_start, ie_end
  logical, optional :: err_flag, zero_uncalculated
end subroutine

subroutine type_coord (coord)
  import
  implicit none
  type (coord_struct) coord
end subroutine

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, &
      type_control, type_wake, type_floor_coords, type_field, type_wall, lines, n_lines)
  import
  implicit none
  type (ele_struct), target :: ele
  integer, optional, intent(in) :: type_mat6
  integer, optional, intent(out) :: n_lines
  integer, optional, intent(in) :: twiss_out
  logical, optional, intent(in) :: type_control, type_taylor, type_floor_coords
  logical, optional, intent(in) :: type_zero_attrib, type_wake
  logical, optional :: type_field, type_wall
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

function value_of_attribute (ele, attrib_name, err_flag, err_print_flag) result (value)
  import
  implicit none
  type (ele_struct), target :: ele
  type (all_pointer_struct) a_ptr
  real(rp) value
  character(*) attrib_name
  logical, optional :: err_print_flag, err_flag
end function

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

subroutine xsif_parser (xsif_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
  import
  implicit none
  character(*) xsif_file
  character(*), optional :: use_line
  type (lat_struct), target :: lat
  logical, optional :: make_mats6
  logical, optional :: digested_read_ok, err_flag
  end subroutine

end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
