!+
! This file defines the interfaces for the BMAD subroutines
!-

module basic_bmad_interface

interface
  subroutine compute_reference_energy (lat, compute)
    use bmad_struct, only: lat_struct
    type (lat_struct) lat
    logical, optional :: compute
  end subroutine
end interface

interface
  subroutine add_superimpose (lat, super_ele, ix_super)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct) lat
    type (ele_struct) super_ele
    integer ix_super
  end subroutine
end interface

interface
  subroutine aml_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)
    use bmad_struct, only: lat_struct
    implicit none
    character(*) lat_file
    type (lat_struct), target :: lat
    logical, optional :: make_mats6
    logical, optional :: digested_read_ok
    character(*), optional :: use_line
  end subroutine
end interface

interface
  subroutine bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)
    use bmad_struct, only: lat_struct
    implicit none
    character(*) lat_file
    type (lat_struct), target :: lat
    logical, optional :: make_mats6
    logical, optional :: digested_read_ok
    character(*), optional :: use_line
  end subroutine
end interface
  
interface
  subroutine bmad_parser2 (in_file, lat, orbit, make_mats6, &
                                      digested_file_name, digested_read_ok)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    character(*) in_file
    type (lat_struct), target :: lat
    type (coord_struct), optional :: orbit(0:)
    character(*), optional :: digested_file_name
    logical, optional :: make_mats6, digested_read_ok
  end subroutine
end interface

interface
  subroutine c_to_cbar (ele, cbar_mat)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct) ele
    real(rp) cbar_mat(2,2)
  end subroutine
end interface

interface
  subroutine cbar_to_c (cbar_mat, a, b, c_mat)
    use bmad_struct, only: twiss_struct, rp
    implicit none
    real(rp) cbar_mat(2,2), c_mat(2,2)
    type (twiss_struct) a, b
  end subroutine
end interface

interface
  subroutine calc_z_tune (lat)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
  end subroutine
end interface

interface
  subroutine check_lat_controls (lat, exit_on_error)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), target :: lat
    logical exit_on_error
  end subroutine
end interface

interface
  subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, exit_on_error)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct) lat
    real(rp) delta_e
    real(rp) chrom_x
    real(rp) chrom_y
    logical, optional :: exit_on_error
  end subroutine
end interface

interface
  subroutine chrom_tune (lat, delta_e, chrom_x, chrom_y, err_tol, err_flag)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct) lat
    real(rp) delta_e
    real(rp) chrom_x
    real(rp) chrom_y
    real(rp) err_tol
    logical err_flag
  end subroutine
end interface

interface
  subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, exit_on_error)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct) lat
    type (coord_struct), allocatable, target :: closed_orb(:)
    integer i_dim
    integer, optional :: direction
    logical, optional :: exit_on_error
  end subroutine
end interface

interface
 subroutine closed_orbit_from_tracking (lat, closed_orb, i_dim, &
                                               eps_rel, eps_abs, init_guess)
   use bmad_struct, only: lat_struct, rp, coord_struct
   type (lat_struct) lat
   type (coord_struct), allocatable :: closed_orb(:)
   type (coord_struct), optional :: init_guess
   real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)
   integer i_dim
 end subroutine
end interface

interface
  subroutine combine_consecutive_elements (lat)
    use bmad_struct, only: lat_struct
    type (lat_struct), target :: lat
  end subroutine
end interface

interface
  subroutine convert_coords (in_type_str, coord_in, ele, out_type_str, coord_out)
    use bmad_struct, only: ele_struct, coord_struct
    implicit none
    character(*) in_type_str
    character(*) out_type_str
    type (coord_struct) coord_in
    type (coord_struct) coord_out
    type (ele_struct) ele
  end subroutine
end interface

interface
  subroutine create_group (lat, ix_ele, con, err, err_print_flag)
    use bmad_struct, only: lat_struct, control_struct
    implicit none
    type (lat_struct) lat
    type (control_struct) con(:)
    integer ix_ele
    logical err
    logical, optional :: err_print_flag
  end subroutine
end interface

interface
  subroutine create_girder (lat, ix_ele, ix_slave, init_ele)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct) lat
    type (ele_struct), optional :: init_ele
    integer, intent(in) :: ix_ele
    integer, intent(in) :: ix_slave(:)
  end subroutine
end interface

interface
  subroutine create_overlay (lat, ix_overlay, attrib_value, contl, err, err_print_flag)
    use bmad_struct, only: lat_struct, control_struct
    implicit none
    type (lat_struct) lat
    integer ix_overlay
    character(*) attrib_value
    type (control_struct) contl(:)
    logical err
    logical, optional :: err_print_flag
  end subroutine
end interface

interface
  subroutine create_unique_ele_names (lat, key, suffix)
    use bmad_struct, only: lat_struct
    type (lat_struct), target :: lat
    integer key
    character(*) suffix
  end subroutine
end interface

interface
  Subroutine dispersion_to_orbit (ele, disp_orb)
    use bmad_struct, only: ele_struct, coord_struct
    implicit none
    type (ele_struct) :: ele
    type (coord_struct), intent(out) :: disp_orb
  end subroutine
end interface

interface
  subroutine do_mode_flip (ele)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele
  end subroutine
end interface

interface
  subroutine element_locator (ele_name, lat, ix_ele)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
    integer ix_ele
    character(*) ele_name
  end subroutine
end interface

interface
  subroutine elements_locator_by_key (key, lat, indx)
    use bmad_struct, only: lat_struct
    implicit none
    integer key
    type (lat_struct) lat
    integer, pointer :: indx(:)
  end subroutine
end interface

interface
  subroutine elements_locator (ele_name, lat, indx, err)
    use bmad_struct, only: lat_struct
    implicit none
    character(*) ele_name
    type (lat_struct) lat
    integer, allocatable :: indx(:)
    logical err
  end subroutine
end interface

interface
  subroutine ele_at_s (lat, s, ix_ele)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct) lat
    real(rp) s
    integer ix_ele
  end subroutine
end interface

interface
  subroutine ele_geometry (floor0, ele, floor)
    use bmad_struct, only: floor_position_struct, ele_struct
    type (ele_struct) ele
    type (floor_position_struct) floor0, floor
  end subroutine
end interface

interface
  subroutine emit_calc (lat, what, mode)
    use bmad_struct, only: lat_struct, normal_modes_struct
    implicit none
    type (lat_struct) lat
    type (normal_modes_struct) mode
    integer what
  end subroutine
end interface

interface
  subroutine find_element_ends (lat, ele, ele1, ele2)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct) lat
    type (ele_struct) ele
    type (ele_struct), pointer :: ele1, ele2
  end subroutine
end interface

interface
  subroutine get_element_lord_list (lat, ix_ele, lord_list, n_lord)
    use bmad_struct, only: lat_struct
    type (lat_struct) lat
    integer ix_ele, n_lord
    integer, allocatable :: lord_list(:)
  end subroutine
end interface

interface
  subroutine get_element_slave_list (lat, ix_ele, slave_list, n_slave)
    use bmad_struct, only: lat_struct
    type (lat_struct) lat
    integer ix_ele, n_slave
    integer, allocatable :: slave_list(:)
  end subroutine
end interface

interface
  subroutine init_LRBBI(lat, oppos_lat, LRBBI_ele, ix_LRBBI, ix_oppos)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct) lat
	type (lat_struct) :: oppos_lat
    type (ele_struct) LRBBI_ele
	integer, intent(in) :: ix_LRBBI, ix_oppos
  end subroutine
end interface

interface
  subroutine insert_element (lat, insert_ele, insert_index)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct) lat
    type (ele_struct) insert_ele
    integer insert_index
  end subroutine
end interface

interface
  subroutine insert_LRBBI (lat, oppos_lat, cross_positions, ix_LRBBI)
	use bmad_struct, only: lat_struct, rp
    type (lat_struct) lat
	type (lat_struct) oppos_lat 
   real(rp), dimension(:), intent(inout) :: cross_positions
    integer, dimension(:), intent(inout) :: ix_LRBBI
  end subroutine
end interface

interface
  subroutine make_g_mats (ele, g_mat, g_inv_mat)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct) ele
    real(rp) g_mat(4,4)
    real(rp) g_inv_mat(4,4)
  end subroutine
end interface

interface
  subroutine make_hybrid_lat (r_in, use_ele, remove_markers, &
                                           r_out, ix_out, use_taylor, orb0)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct), target :: r_in
    type (lat_struct), target :: r_out
    integer, optional :: ix_out(:)
    logical remove_markers
    logical use_ele(:)
    logical, optional :: use_taylor
    type (coord_struct), optional :: orb0(0:)
  end subroutine
end interface

interface
  subroutine make_LRBBI(master_lat_oppos, lat, ix_LRBBI, master_ix_LRBBI)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), dimension(:) :: lat
    type (lat_struct) :: master_lat_oppos
    integer, dimension(:,:) :: ix_LRBBI
    integer, dimension(:,:) :: master_ix_LRBBI
  end subroutine
end interface
 
interface
  subroutine make_mat6 (ele, param, start, end, end_in, err)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct) ele
    type (coord_struct), optional :: start, end
    type (lat_param_struct) param
    logical, optional :: end_in
    logical, optional :: err
  end subroutine
end interface

interface
  subroutine make_mat6_taylor (ele, param, start)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start
    type (lat_param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_bmad (ele, param, start, end, end_in, err)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (lat_param_struct) param
    logical, optional :: end_in
    logical, optional :: err
  end subroutine
end interface

interface
  subroutine make_mat6_runge_kutta (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (lat_param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_symp_lie_ptc (ele, param, start)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start
    type (lat_param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_tracking (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (lat_param_struct) param
  end subroutine
end interface

interface
  subroutine make_v_mats (ele, v_mat, v_inv_mat)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct) ele
    real(rp) v_mat(4,4)
    real(rp) v_inv_mat(4,4)
  end subroutine
end interface

interface
  subroutine mat6_add_offsets (ele)
    use bmad_struct, only: ele_struct
    type (ele_struct) ele
  end subroutine
end interface

interface
  subroutine mark_LRBBI(master_lat, master_lat_oppos, lat, crossings)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct), dimension(:) :: lat
    type (lat_struct) :: master_lat
    type (lat_struct) :: master_lat_oppos
    real(rp), dimension(:,:) :: crossings
  end subroutine
end interface

interface
  subroutine name_to_list (lat, ele_names, use_ele)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
    logical use_ele(:)
    character(*) ele_names(:)
  end subroutine
end interface

interface
  subroutine new_control (lat, ix_ele)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
    integer ix_ele
  end subroutine
end interface

interface
  subroutine offset_particle (ele, param, coord, set, &
           set_canonical, set_tilt, set_multipoles, set_hvkicks, s_pos)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, rp
    implicit none
    type (ele_struct) :: ele
    type (lat_param_struct), intent(in) :: param
    type (coord_struct), intent(inout) :: coord
    logical, intent(in) :: set
    logical, optional, intent(in) :: set_canonical, set_multipoles, &
                                                      set_tilt, set_hvkicks
    real(rp), optional, intent(in) :: s_pos
  end subroutine
end interface

interface
  subroutine transfer_matrix_calc (lat, rf_on, xfer_mat, xfer_vec, ix1, ix2)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct) lat
    logical, intent(in) :: rf_on
    real(rp), intent(out) :: xfer_mat(:,:)
    real(rp), intent(out), optional :: xfer_vec(:)
    integer, intent(in), optional :: ix1, ix2
  end subroutine
end interface

interface
  subroutine transfer_map_calc (lat, t_map, ix1, ix2, &
                                           integrate, one_turn, unit_start)
    use bmad_struct, only: lat_struct, rp, taylor_struct
    implicit none
    type (lat_struct) lat
    type (taylor_struct) :: t_map(:)
    integer, intent(in), optional :: ix1, ix2
    logical, optional :: integrate, one_turn, unit_start
  end subroutine
end interface

interface
  subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)
    use bmad_struct, only: ele_struct, rp
    type (ele_struct) ele
    real(rp) phi_a
    real(rp) phi_b
    real(rp) mat4(4,4)
  end subroutine
end interface

interface
  subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele, &
                                                    stable, growth_rate, chi)
    use bmad_struct, only: ele_struct, rp, coord_struct
    implicit none
    type (coord_struct), intent(in) :: track(:)
    type (coord_struct), intent(out) :: track0
    type (ele_struct) :: ele
    real(rp), intent(out) :: growth_rate, chi
    integer, intent(in) :: i_dim
    logical, intent(out) :: stable
  end subroutine
end interface

interface
  subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, map0, track0, chi)
    use bmad_struct, only: rp, coord_struct
    implicit none
    type (coord_struct), intent(in), target :: track(:)
    type (coord_struct), intent(out) :: track0
    real(rp), intent(out) :: mat1(:,:), map0(:)
    real(rp), intent(out) :: chi
    integer, intent(in) :: i_dim
  end subroutine
end interface

interface
  Subroutine orbit_to_dispersion (orb_diff, ele)
    use bmad_struct, only: ele_struct, coord_struct
    implicit none
    type (coord_struct), intent(in) :: orb_diff
    type (ele_struct) :: ele
  end subroutine
end interface

interface
  subroutine order_super_lord_slaves (lat, ix_lord)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), target :: lat
    integer ix_lord
  end subroutine
end interface

interface
  subroutine phase_space_fit (x, xp, twiss, tune, emit, x_0, xp_0, chi, tol)
    use bmad_struct, only: rp, twiss_struct
    implicit none
    type (twiss_struct) twiss
    real(rp), optional :: tol
    real(rp) x(:), xp(:)
    real(rp) tune, emit
    real(rp) x_0, xp_0, chi
  end subroutine
end interface

interface
  Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, &
                    ptr_attrib, err_flag, err_print_flag, ix_attrib)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct), target :: ele
    real(rp), pointer :: ptr_attrib
    character(*) attrib_name
    logical err_flag
    logical do_allocation
    logical, optional :: err_print_flag
    integer, optional :: ix_attrib
  end subroutine
end interface

interface
  Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                    ptr_array, err_flag, err_print_flag, ix_eles, ix_attrib)
    use bmad_struct, only: lat_struct, real_pointer_struct
    implicit none
    type (lat_struct) lat
    type (real_pointer_struct), allocatable :: ptr_array(:)
    character(*) ele_name, attrib_name
    logical err_flag
    logical do_allocation
    logical, optional :: err_print_flag
    integer, optional, allocatable :: ix_eles(:)
    integer, optional :: ix_attrib
  end subroutine
end interface

interface
  function attribute_free (ix_ele, ix_attrib, lat, err_print_flag, except_overlay) result (free)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct), target :: lat
    integer ix_ele, ix_attrib
    logical free
    logical, optional :: err_print_flag, except_overlay
  end function
end interface

interface
  subroutine quad_beta_ave (lat, ix_ele, beta_a_ave, beta_b_ave)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct) lat
    integer ix_ele
    real(rp) beta_a_ave
    real(rp) beta_b_ave
  end subroutine
end interface

interface
  subroutine radiation_integrals (lat, orb, mode, ix_cache)
    use bmad_struct, only: lat_struct, coord_struct, normal_modes_struct
    implicit none
    type (lat_struct), target :: lat
    type (coord_struct), target :: orb(0:)
    type (normal_modes_struct) mode
    integer, optional :: ix_cache
  end subroutine
end interface

interface
  subroutine remove_eles_from_lat (lat, check_controls)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
    logical, optional :: check_controls
  end subroutine
end interface

interface
  subroutine read_digested_bmad_file (in_file_name, lat, version)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), target, intent(inout) :: lat
    integer version
    character(*) in_file_name
  end subroutine
end interface

interface
  function relative_mode_flip (ele1, ele2)
    use bmad_struct, only: ele_struct
    implicit none
    logical relative_mode_flip
    type (ele_struct) ele1
    type (ele_struct) ele2
  end function
end interface

interface
  subroutine lat_geometry (lat)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
  end subroutine
end interface

interface
  recursive subroutine lat_make_mat6 (lat, ix_ele, coord)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct), target :: lat
    type (coord_struct), optional :: coord(0:)
    integer, optional :: ix_ele
  end subroutine
end interface

interface
  subroutine s_calc (lat)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
  end subroutine
end interface

interface
  subroutine set_on (key, lat, on_switch, orb)
    use bmad_struct, only: lat_struct, coord_struct
    type (lat_struct) lat
    type (coord_struct), optional :: orb(0:)
    integer key
    logical on_switch
  end subroutine
end interface

interface
  subroutine set_tune (phi_a_set, phi_b_set, dk1, lat, orb, ok)
    use bmad_struct, only: lat_struct, rp, coord_struct
    implicit none
    type (lat_struct) lat
    type (coord_struct), allocatable :: orb(:)
    real(rp) phi_a_set
    real(rp) phi_b_set
    real(rp) dk1(:)
    logical ok
  end subroutine
end interface

interface
  subroutine set_z_tune (lat, z_tune)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct), target :: lat
    real(rp), optional :: z_tune
  end subroutine
end interface

interface
  subroutine slice_ele_calc (ele, param, i_slice, n_slice_tot, sliced_ele)
    use bmad_struct, only: ele_struct, lat_param_struct
    implicit none
    type (ele_struct) ele, sliced_ele
    type (lat_param_struct) param
    integer i_slice, n_slice_tot
  end subroutine
end interface

interface
  subroutine split_lat (lat, s_split, ix_split, split_done, add_suffix, check_controls)
    use bmad_struct, only: lat_struct, rp
    implicit none
    type (lat_struct), target :: lat
    real(rp) s_split
    integer ix_split
    logical split_done
    logical, optional :: add_suffix, check_controls
  end subroutine
end interface

interface
  subroutine tilt_coords (tilt_val, coord, set)
    use precision_def
    implicit none
    real(rp) tilt_val
    real(rp) coord(:)
    logical set
  end subroutine
end interface

interface
  subroutine track_all (lat, orbit)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct) lat
    type (coord_struct), allocatable :: orbit(:)
  end subroutine
end interface

interface
  subroutine track_many (lat, orbit, ix_start, ix_end, direction)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct)  lat
    type (coord_struct)  orbit(0:)
    integer ix_start
    integer ix_end
    integer direction
  end subroutine
end interface

interface
  subroutine track1 (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct)   :: ele
    type (lat_param_struct) :: param
  end subroutine
end interface

interface
  subroutine track1_runge_kutta (start, ele, param, end, track)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, track_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct), target :: ele
    type (lat_param_struct), target :: param
    type (track_struct) track
  end subroutine
end interface

interface
  subroutine track1_linear (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
  end subroutine
end interface

interface
  subroutine track1_taylor (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
  end subroutine
end interface

interface
  subroutine track1_bmad (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
  end subroutine
end interface

interface
  subroutine track1_symp_lie_ptc (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
  end subroutine
end interface

interface
  subroutine track1_symp_map (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
  end subroutine
end interface

interface
  subroutine twiss_and_track (lat, orb, ok)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct) lat
    type (coord_struct), allocatable :: orb(:)
    logical, optional :: ok
  end subroutine
end interface

interface
  subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
                                                  start, end, body_only, err)
    use bmad_struct, only: coord_struct, ele_struct, lat_param_struct, rp
    implicit none
    type (ele_struct), optional :: ele3
    type (ele_struct) ele1
    type (ele_struct) ele2
    type (coord_struct), optional :: start
    type (coord_struct), optional :: end
    type (lat_param_struct) param
    logical, optional :: body_only, err
    real(rp) del_s
  end subroutine
end interface

interface
  subroutine twiss_and_track_body (ele1, ele2, param, del_s, ele3, &
                                                                 start, end)
    use bmad_struct, only: ele_struct, rp, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), optional :: ele3
    type (ele_struct) ele1
    type (ele_struct) ele2
    type (coord_struct), optional :: start
    type (coord_struct), optional :: end
    type (lat_param_struct) param
    real(rp) del_s
  end subroutine
end interface

interface
  recursive subroutine twiss_at_element (lat, ix_ele, start, end, average)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (lat_struct), target :: lat
    type (ele_struct), optional :: start
    type (ele_struct), optional :: end
    type (ele_struct), optional :: average
    integer ix_ele
  end subroutine
end interface

interface
  subroutine twiss_and_track_at_s (lat, s, ele, orb, here, err)
    use bmad_struct, only: lat_struct, ele_struct, rp, coord_struct
    implicit none
    type (lat_struct) :: lat
    type (ele_struct), optional :: ele
    real(rp) s
    type (coord_struct), optional :: orb(0:)
    type (coord_struct), optional :: here
    logical, optional :: err
  end subroutine
end interface

interface
  subroutine twiss_at_start (lat)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
  end subroutine
end interface

interface
  subroutine twiss1_propagate (twiss1, mat2, length, twiss2, err)
    use bmad_struct, only: twiss_struct, rp
    implicit none
    type (twiss_struct) twiss1, twiss2
    real(rp) mat2(2,2), length
    logical err
  end subroutine
end interface

interface
  subroutine twiss_from_mat6 (mat6, map0, ele, stable, growth_rate)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct) :: ele
    real(rp), intent(in) :: mat6(:,:), map0(:)
    real(rp), intent(out) :: growth_rate
    logical, intent(out) :: stable
  end subroutine
end interface

interface
  subroutine twiss_from_tracking (lat, ref_orb0, error, d_orb) 
    use bmad_struct, only: lat_struct, rp, coord_struct
    type (lat_struct), intent(inout) :: lat
    type (coord_struct), intent(in) :: ref_orb0
    real(rp), intent(in), optional :: d_orb(:)   
    real(rp), intent(out) :: error
  end subroutine
end interface

interface
  subroutine twiss_propagate1 (ele1, ele2, err)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele1
    type (ele_struct) ele2
    logical, optional :: err
  end subroutine
end interface

interface
  subroutine twiss_propagate_all (lat)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
  end subroutine
end interface

interface
  subroutine twiss_propagate_many (lat, ix_start, ix_end, direction)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) :: lat
    integer, intent(in) :: ix_start
    integer, intent(in) :: ix_end
    integer, intent(in) :: direction
  end subroutine
end interface

interface
  subroutine type_coord (coord)
    use bmad_struct, only: coord_struct
    implicit none
    type (coord_struct) coord
  end subroutine
end interface

interface
  subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, &
        twiss_type, type_control, lattice, type_wake, type_floor_coords, type_wig_terms, nunit)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (ele_struct) ele
    type (lat_struct), optional :: lattice
    integer, optional :: type_mat6
    integer, optional :: twiss_type
    logical, optional :: type_zero_attrib, type_floor_coords
    logical, optional :: type_control, type_taylor, type_wake
    logical, optional :: type_wig_terms
    integer, optional :: nunit
  end subroutine
end interface

interface
  subroutine type_twiss (ele, frequency_units)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele
    integer, optional :: frequency_units
  end subroutine
end interface
 
interface
  subroutine type2_ele (ele, lines, n_lines, type_zero_attrib, type_mat6, &
        type_taylor, twiss_type, type_control, lattice, type_wake, &
        type_floor_coords, type_wig_terms)
    use bmad_struct, only: lat_struct, ele_struct
    implicit none
    type (ele_struct), target :: ele
    type (lat_struct), optional :: lattice
    integer, optional, intent(in) :: type_mat6
    integer, intent(out) :: n_lines
    integer, optional, intent(in) :: twiss_type
    logical, optional, intent(in) :: type_control, type_taylor, type_floor_coords
    logical, optional, intent(in) :: type_zero_attrib, type_wake
    logical, optional :: type_wig_terms
    character(100), pointer :: lines(:)
  end subroutine
end interface

interface
  subroutine type2_twiss (ele, lines, n_lines, frequency_units, compact_format)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele
    integer, optional :: frequency_units
    integer n_lines
    character(*) lines(:)
    logical, optional :: compact_format
  end subroutine
end interface
 
interface
  subroutine write_digested_bmad_file (digested_name, lat,  &
                                                    n_files, file_names)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), target, intent(in) :: lat
    integer, intent(in), optional :: n_files
    character(*) digested_name
    character(*), optional :: file_names(:)
  end subroutine
end interface

interface
  recursive subroutine update_hybrid_list (lat, n_in, use_ele)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
    logical use_ele(:)
    integer n_in
  end subroutine
end interface

interface
  subroutine add_lattice_control_structs (lat, ix_ele)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), target :: lat
    integer ix_ele
  end subroutine
end interface

interface
  subroutine delete_lattice_control_struct (lat, ix_con)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct), target :: lat
    integer ix_con
  end subroutine
end interface

interface
  subroutine orbit_amplitude_calc (ele, orb, amp_a, amp_b, &
                                            amp_na, amp_nb, particle)
    use bmad_struct, only: ele_struct, coord_struct, rp
    implicit none
    type (ele_struct) ele
    type (coord_struct) orb
    integer, optional :: particle
    real(rp), optional :: amp_a, amp_b, amp_na, amp_nb
  end subroutine
end interface

interface
  subroutine set_design_linear (lat)
    use bmad_struct, only: lat_struct
    implicit none
    type (lat_struct) lat
  end subroutine
end interface

interface
  subroutine xsif_parser (xsif_file, lat, make_mats6, digested_read_ok, use_line)
    use bmad_struct, only: lat_struct
    implicit none
    character(*) xsif_file
    character(*), optional :: use_line
    type (lat_struct), target :: lat
    logical, optional :: make_mats6
    logical, optional :: digested_read_ok
  end subroutine
end interface

interface
  subroutine bbi_kick_matrix (ele, param, orb, s_pos, mat6)
    use bmad_struct, only: ele_struct, rp, coord_struct, lat_param_struct
    implicit none
    type (ele_struct) ele
    type (lat_param_struct) param
    type (coord_struct) orb
    real(rp) s_pos
    real(rp) mat6(6,6)
  end subroutine
end interface

end module
