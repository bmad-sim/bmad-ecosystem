!+
! This file defines the interfaces for the BMAD subroutines
!-

#include "CESR_platform.inc"

module bmad_interface

use matrix_mod
use bmad_basic_mod
use equal_mod
use nrutil, only: reallocate

!---------------------------------------------------------

interface
  subroutine lost_particle_info (lattice, orbit, ix_lost, plane_lost)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct) lattice
    type (coord_struct) orbit(0:)
    integer ix_lost, plane_lost
  end subroutine
end interface

interface
  subroutine compute_element_energy (ring)
    use bmad_struct, only: ring_struct
    type (ring_struct) ring
  end subroutine
end interface

interface
  subroutine add_superimpose (ring, super_ele, ix_super)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ring_struct) ring
    type (ele_struct) super_ele
    integer ix_super
  end subroutine
end interface

interface
  function attribute_index (ele, name)
    use bmad_struct, only: ele_struct
    implicit none
    integer attribute_index
    type (ele_struct) ele
    character(*) name
  end function
end interface

interface
  function attribute_name(ele, index)
    use bmad_struct, only: ele_struct
    implicit none
    character(16) attribute_name
    type (ele_struct) ele
    integer index
  end function
end interface

interface
  subroutine bmad_parser (in_file, ring, make_mats6, digested_read_ok)
    use bmad_struct, only: ring_struct
    implicit none
    character(*) in_file
    type (ring_struct), target :: ring
    logical, optional :: make_mats6
    logical, optional :: digested_read_ok
  end subroutine
end interface
  
interface
  subroutine bmad_parser2 (in_file, ring, orbit_, make_mats6)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    character(*) in_file
    type (ring_struct), target :: ring
    type (coord_struct), optional :: orbit_(0:)
    logical, optional :: make_mats6
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
  subroutine cbar_to_c (cbar_mat, ele)
    use bmad_struct, only: ele_struct, rp
    implicit none
    real(rp) cbar_mat(2,2)
    type (ele_struct) ele
  end subroutine
end interface

interface
  subroutine calc_z_tune (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
  end subroutine
end interface

interface
  subroutine check_ring_controls (ring, exit_on_error)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target :: ring
    logical exit_on_error
  end subroutine
end interface

interface
  subroutine chrom_calc (ring, delta_e, chrom_x, chrom_y)
    use bmad_struct, only: ring_struct, rp
    implicit none
    type (ring_struct) ring
    real(rp) delta_e
    real(rp) chrom_x
    real(rp) chrom_y
  end subroutine
end interface

interface
  subroutine chrom_tune (ring, delta_e, chrom_x, chrom_y, err_tol, err_flag)
    use bmad_struct, only: ring_struct, rp
    implicit none
    type (ring_struct) ring
    real(rp) delta_e
    real(rp) chrom_x
    real(rp) chrom_y
    real(rp) err_tol
    logical err_flag
  end subroutine
end interface

interface
  subroutine closed_orbit_at_start (ring, co, i_dim, iterate)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct) ring
    type (coord_struct) co
    integer i_dim
    logical iterate
  end subroutine
end interface

interface
  subroutine closed_orbit_calc (ring, closed_orb, i_dim, direction)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct) ring
    type (coord_struct), allocatable, target :: closed_orb(:)
    integer i_dim
    integer, optional :: direction
  end subroutine
end interface

interface
 subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, &
                                               eps_rel, eps_abs, init_guess)
   use bmad_struct, only: ring_struct, rp, coord_struct
   type (ring_struct) ring
   type (coord_struct), allocatable :: closed_orb_(:)
   type (coord_struct), optional :: init_guess
   real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)
   integer i_dim
 end subroutine
end interface

interface
  subroutine compress_ring (ring, ok)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target :: ring
    logical ok
  end subroutine
end interface

interface
  subroutine convert_coords (in_type_str, coord_in, ele, out_type_str, &
                                                                 coord_out)
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
  subroutine create_group (ring, ix_ele, con_)
    use bmad_struct, only: ring_struct, control_struct
    implicit none
    type (ring_struct) ring
    type (control_struct) con_(:)
    integer ix_ele
  end subroutine
end interface

interface
  subroutine create_i_beam (ring, ix_ele, ix_slave, init_ele)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ring_struct) ring
    type (ele_struct), optional :: init_ele
    integer, intent(in) :: ix_ele
    integer, intent(in) :: ix_slave(:)
  end subroutine
end interface

interface
  subroutine create_overlay (ring, ix_overlay, attrib_value, contl)
    use bmad_struct, only: ring_struct, control_struct
    implicit none
    type (ring_struct) ring
    integer ix_overlay
    character(*) attrib_value
    type (control_struct) contl(:)
  end subroutine
end interface

interface
  subroutine custom_emitt_calc (ring, ir, i2, i3, i5a, i5b)
    use bmad_struct, only: ring_struct, rp
    type (ring_struct) ring
    integer ir
    real(rp) i2, i3, i5a, i5b
  end subroutine
end interface

interface
  Subroutine dispersion_to_orbit (ele, disp_orb)
    use bmad_struct, only: ele_struct, coord_struct
    implicit none
    type (ele_struct), intent(in) :: ele
    type (coord_struct), intent(out) :: disp_orb
  end subroutine
end interface

interface
  subroutine do_mode_flip (ele, ele_flip)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele
    type (ele_struct) ele_flip
  end subroutine
end interface

interface
  subroutine element_locator (ele_name, ring, ix_ele)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
    integer ix_ele
    character(*) ele_name
  end subroutine
end interface

interface
  subroutine elements_locator (key, ring, indx)
    use bmad_struct, only: ring_struct
    implicit none
    integer key
    type (ring_struct) ring
    integer, pointer :: indx(:)
  end subroutine
end interface

interface
  subroutine emitt_calc (ring, what, mode)
    use bmad_struct, only: ring_struct, modes_struct
    implicit none
    type (ring_struct) ring
    type (modes_struct) mode
    integer what
  end subroutine
end interface

interface
  subroutine find_element_ends (ring, ix_ele, ix_start, ix_end)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
    integer ix_ele
    integer ix_start
    integer ix_end
  end subroutine
end interface

interface
  subroutine get_lattice_list (lat_list, num_lats, directory)
    use precision_def
    implicit none
    integer num_lats
    character(*) directory
    character(40) lat_list(:)
  end subroutine
end interface

interface
  subroutine init_LRBBI(ring, oppos_ring, LRBBI_ele, ix_LRBBI, ix_oppos)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ring_struct) ring
	type (ring_struct) :: oppos_ring
    type (ele_struct) LRBBI_ele
	integer, intent(in) :: ix_LRBBI, ix_oppos
  end subroutine
end interface

interface
  subroutine insert_element (ring, insert_ele, insert_index)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ring_struct) ring
    type (ele_struct) insert_ele
    integer insert_index
  end subroutine
end interface

interface
  subroutine insert_LRBBI (ring, oppos_ring, cross_positions, ix_LRBBI)
	use bmad_struct, only: ring_struct, rp
    type (ring_struct) ring
	type (ring_struct) oppos_ring 
   real(rp), dimension(:), intent(inout) :: cross_positions
    integer, dimension(:), intent(inout) :: ix_LRBBI
  end subroutine
end interface

interface
  subroutine lattice_to_bmad_file_name (lattice, bmad_file_name)
    implicit none
    character(*) lattice
    character(*) bmad_file_name
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
  subroutine make_hybrid_ring (r_in, use_ele, remove_markers, &
                                           r_out, ix_out, use_taylor, orb0_)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct), target :: r_in
    type (ring_struct), target :: r_out
    integer ix_out(:)
    logical remove_markers
    logical use_ele(:)
    logical, optional :: use_taylor
    type (coord_struct), optional :: orb0_(0:)
  end subroutine
end interface

interface
  subroutine make_LRBBI(master_ring_oppos, ring, ix_LRBBI, master_ix_LRBBI)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), dimension(:) :: ring
    type (ring_struct) :: master_ring_oppos
    integer, dimension(:,:) :: ix_LRBBI
    integer, dimension(:,:) :: master_ix_LRBBI
  end subroutine
end interface
 
interface
  subroutine make_mat6 (ele, param, start, end, end_in)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct) ele
    type (coord_struct), optional :: start, end
    type (param_struct) param
    logical, optional :: end_in
  end subroutine
end interface

interface
  subroutine make_mat6_custom (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_taylor (ele, param, start)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start
    type (param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_bmad (ele, param, start, end, end_in)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (param_struct) param
    logical, optional :: end_in
  end subroutine
end interface

interface
  subroutine make_mat6_runge_kutta (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_symp_lie_ptc (ele, param, start)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start
    type (param_struct) param
  end subroutine
end interface

interface
  subroutine make_mat6_tracking (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (param_struct) param
  end subroutine
end interface

interface
  subroutine custom_radiation_integrals (ring, ir, orb_)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct) ring
    type (coord_struct) orb_(0:)
    integer ir
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
  subroutine mark_LRBBI(master_ring, master_ring_oppos, ring, crossings)
    use bmad_struct, only: ring_struct, rp
    implicit none
    type (ring_struct), dimension(:) :: ring
    type (ring_struct) :: master_ring
    type (ring_struct) :: master_ring_oppos
    real(rp), dimension(:,:) :: crossings
  end subroutine
end interface

interface
  subroutine name_to_list (ring, ele_names, use_ele)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
    logical use_ele(:)
    character(*) ele_names(:)
  end subroutine
end interface

interface
  subroutine new_control (ring, ix_ele)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
    integer ix_ele
  end subroutine
end interface

interface
  subroutine offset_particle (ele, param, coord, set, &
           set_canonical, set_tilt, set_multipoles, set_hvkicks, s_pos)
    use bmad_struct, only: ele_struct, coord_struct, param_struct, rp
    implicit none
    type (ele_struct), intent(in) :: ele
    type (param_struct), intent(in) :: param
    type (coord_struct), intent(inout) :: coord
    logical, intent(in) :: set
    logical, optional, intent(in) :: set_canonical, set_multipoles, &
                                                      set_tilt, set_hvkicks
    real(rp), optional, intent(in) :: s_pos
  end subroutine
end interface

interface
  subroutine transfer_matrix_calc (lat, rf_on, t_mat, ix1, ix2)
    use bmad_struct, only: ring_struct, rp
    implicit none
    type (ring_struct) lat
    logical, intent(in) :: rf_on
    real(rp), intent(out) :: t_mat(:,:)
    integer, intent(in), optional :: ix1, ix2
  end subroutine
end interface

interface
  subroutine transfer_map_calc (lat, t_map, ix1, ix2)
    use bmad_struct, only: ring_struct, rp, taylor_struct
    implicit none
    type (ring_struct) lat
    type (taylor_struct) :: t_map(:)
    integer, intent(in), optional :: ix1, ix2
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
  subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)
    use bmad_struct, only: rp, coord_struct
    implicit none
    type (coord_struct), intent(in), target :: track(:)
    type (coord_struct), intent(out) :: track0
    real(rp), intent(out) :: mat1(:,:)
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
  subroutine order_super_lord_slaves (ring, ix_lord)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target :: ring
    integer ix_lord
  end subroutine
end interface

interface
  subroutine phase_space_fit (x, xp, twiss, tune, emitt, x_0, xp_0, chi, tol)
    use bmad_struct, only: rp, twiss_struct
    implicit none
    type (twiss_struct) twiss
    real(rp), optional :: tol
    real(rp) x(:), xp(:)
    real(rp) tune, emitt
    real(rp) x_0, xp_0, chi
  end subroutine
end interface

interface
  Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, &
                    ptr_attrib, ix_attrib, cannot_vary_flag, err_print_flag)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct), target :: ele
    real(rp), pointer :: ptr_attrib
    integer ix_attrib
    character(*) attrib_name
    logical cannot_vary_flag
    logical do_allocation
    logical, optional :: err_print_flag
  end subroutine
end interface

interface
  subroutine check_attrib_free (ele, ix_attrib, ring, &
                                                err_flag, err_print_flag)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ring_struct) :: ring
    type (ele_struct) :: ele
    integer ix_attrib
    logical err_flag
    logical, optional :: err_print_flag
  end subroutine
end interface

interface
  subroutine quad_beta_ave (ring, ix_ele, beta_x_ave, beta_y_ave)
    use bmad_struct, only: ring_struct, rp
    implicit none
    type (ring_struct) ring
    integer ix_ele
    real(rp) beta_x_ave
    real(rp) beta_y_ave
  end subroutine
end interface

interface
  subroutine radiation_integrals (ring, orb_, mode, ix_cache)
    use bmad_struct, only: ring_struct, coord_struct, modes_struct
    implicit none
    type (ring_struct), target :: ring
    type (coord_struct), target :: orb_(0:)
    type (modes_struct) mode
    integer, optional :: ix_cache
  end subroutine
end interface

interface
  subroutine read_digested_bmad_file (in_file_name, ring, version)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target, intent(inout) :: ring
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
  subroutine ring_geometry (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
  end subroutine
end interface

interface
  recursive subroutine ring_make_mat6 (ring, ix_ele, coord_)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct), target :: ring
    type (coord_struct), optional :: coord_(0:)
    integer ix_ele
  end subroutine
end interface

interface
  subroutine set_ele_attribute (ring, i_ele, attrib_name, &
                              attrib_value, err_flag, make_mat6_flag, orbit_)
    use bmad_struct, only: ring_struct, rp, coord_struct
    implicit none
    type (ring_struct) :: ring
    type (coord_struct), optional :: orbit_(0:)
    real(rp) attrib_value
    integer i_ele
    character(*) attrib_name
    logical make_mat6_flag
    logical err_flag
  end subroutine
end interface

interface
  subroutine s_calc (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
  end subroutine
end interface

interface
  subroutine set_on (key, ring, on_switch, orb_)
    use bmad_struct, only: ring_struct, coord_struct
    type (ring_struct) ring
    type (coord_struct), optional :: orb_(0:)
    integer key
    logical on_switch
  end subroutine
end interface

interface
  subroutine set_tune (phi_x_set, phi_y_set, dk1, ring, orb_, ok)
    use bmad_struct, only: ring_struct, rp, coord_struct
    implicit none
    type (ring_struct) ring
    type (coord_struct), allocatable :: orb_(:)
    real(rp) phi_x_set
    real(rp) phi_y_set
    real(rp) dk1(:)
    logical ok
  end subroutine
end interface

interface
  subroutine set_z_tune (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target :: ring
  end subroutine
end interface

interface
  subroutine split_ring (ring, s_split, ix_split, split_done)
    use bmad_struct, only: ring_struct, rp
    implicit none
    type (ring_struct), target :: ring
    real(rp) s_split
    integer ix_split
    logical split_done
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
  subroutine track_all (ring, orbit_)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct) ring
    type (coord_struct), allocatable :: orbit_(:)
  end subroutine
end interface

interface
  subroutine track_many (ring, orbit_, ix_start, ix_end, direction)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct)  ring
    type (coord_struct)  orbit_(0:)
    integer ix_start
    integer ix_end
    integer direction
  end subroutine
end interface

interface
  subroutine track1 (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_runge_kutta (start, ele, param, end, track)
    use bmad_struct, only: ele_struct, coord_struct, param_struct, track_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), target, intent(inout) :: ele
    type (param_struct), target, intent(inout) :: param
    type (track_struct) track
  end subroutine
end interface

interface
  subroutine track1_linear (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_taylor (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_custom (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_bmad (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_symp_lie_ptc (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_symp_map (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct), intent(in) :: start
    type (coord_struct), intent(out) :: end
    type (ele_struct), intent(inout) :: ele
    type (param_struct), intent(inout) :: param
  end subroutine
end interface

interface
  subroutine track1_wiedemann_wiggler (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, param_struct
    implicit none
    type (coord_struct) start
    type (coord_struct) end
    type (param_struct) param
    type (ele_struct) ele
  end subroutine
end interface

interface
  subroutine twiss_and_track (ring, orb)
    use bmad_struct, only: ring_struct, coord_struct
    implicit none
    type (ring_struct) ring
    type (coord_struct), allocatable :: orb(:)
  end subroutine
end interface

interface
  subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
                                                       start, end, body_only)
    use bmad_struct, only: coord_struct, ele_struct, param_struct, rp
    implicit none
    type (ele_struct), optional :: ele3
    type (ele_struct) ele1
    type (ele_struct) ele2
    type (coord_struct), optional :: start
    type (coord_struct), optional :: end
    type (param_struct) param
    logical, optional :: body_only
    real(rp) del_s
  end subroutine
end interface

interface
  subroutine twiss_and_track_body (ele1, ele2, param, del_s, ele3, &
                                                                 start, end)
    use bmad_struct, only: ele_struct, rp, coord_struct, param_struct
    implicit none
    type (ele_struct), optional :: ele3
    type (ele_struct) ele1
    type (ele_struct) ele2
    type (coord_struct), optional :: start
    type (coord_struct), optional :: end
    type (param_struct) param
    real(rp) del_s
  end subroutine
end interface

interface
  subroutine twiss_at_element (ring, ix_ele, start, end, average)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ring_struct), target :: ring
    type (ele_struct), optional :: start
    type (ele_struct), optional :: end
    type (ele_struct), optional :: average
    integer ix_ele
  end subroutine
end interface

interface
  subroutine twiss_at_s (ring, s, ele)
    use bmad_struct, only: ring_struct, ele_struct, rp
    implicit none
    type (ring_struct) :: ring
    type (ele_struct) :: ele
    real(rp) s
  end subroutine
end interface

interface
  subroutine twiss_and_track_at_s (ring, s, ele, orb_, here)
    use bmad_struct, only: ring_struct, ele_struct, rp, coord_struct
    implicit none
    type (ring_struct) :: ring
    type (ele_struct) :: ele
    real(rp) s
    type (coord_struct), optional :: orb_(0:)
    type (coord_struct), optional :: here
  end subroutine
end interface

interface
  subroutine twiss_at_start (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
  end subroutine
end interface

interface
  subroutine twiss_from_mat6 (mat6, ele, stable, growth_rate)
    use bmad_struct, only: ele_struct, rp
    implicit none
    type (ele_struct) :: ele
    real(rp), intent(in) :: mat6(6,6)
    real(rp), intent(out) :: growth_rate
    logical, intent(out) :: stable
  end subroutine
end interface

interface
  subroutine twiss_from_tracking (ring, ref_orb0, error, d_orb) 
    use bmad_struct, only: ring_struct, rp, coord_struct
    type (ring_struct), intent(inout) :: ring
    type (coord_struct), intent(in) :: ref_orb0
    real(rp), intent(in), optional :: d_orb(:)   
    real(rp), intent(out) :: error
  end subroutine
end interface

interface
  subroutine twiss_propagate1 (ele1, ele2)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele1
    type (ele_struct) ele2
  end subroutine
end interface

interface
  subroutine twiss_propagate_all (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
  end subroutine
end interface

interface
  subroutine twiss_propagate_many (ring, ix_start, ix_end, direction)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) :: ring
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
                                    twiss_type, type_control, ring)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ele_struct) ele
    type (ring_struct), optional :: ring
    integer type_mat6
    integer twiss_type
    logical type_control, type_taylor
    logical type_zero_attrib
  end subroutine
end interface

interface
  subroutine type_twiss (ele, frequency_units)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele
    integer frequency_units
  end subroutine
end interface
 
interface
  subroutine type2_ele (ele, type_zero_attrib, type_mat6, type_taylor, &
                          twiss_type, type_control, lines, n_lines, ring)
    use bmad_struct, only: ring_struct, ele_struct
    implicit none
    type (ele_struct), target, intent(in) :: ele
    type (ring_struct), optional, intent(in) :: ring
    integer, intent(in) :: type_mat6
    integer, intent(out) :: n_lines
    integer, intent(in) :: twiss_type
    logical, intent(in) :: type_control, type_taylor
    logical, intent(in) :: type_zero_attrib
    character(80), pointer :: lines(:)
  end subroutine
end interface

interface
  subroutine type2_twiss (ele, frequency_units, lines, n_lines)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct) ele
    integer frequency_units
    integer n_lines
    character(*) lines(:)
  end subroutine
end interface
 
interface
  subroutine write_digested_bmad_file (digested_name, ring,  &
                                                    n_files, file_names)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target, intent(in) :: ring
    integer, intent(in), optional :: n_files
    character(*) digested_name
    character(*), optional :: file_names(:)
  end subroutine
end interface

interface
  recursive subroutine update_hybrid_list (ring, n_in, use_ele)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
    logical use_ele(:)
    integer n_in
  end subroutine
end interface

interface
  subroutine adjust_control_struct (ring, ix_ele)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct), target :: ring
    integer ix_ele
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
  subroutine set_design_linear (ring)
    use bmad_struct, only: ring_struct
    implicit none
    type (ring_struct) ring
  end subroutine
end interface

interface
  subroutine xsif_parser (xsif_file, ring, make_mats6)
    use bmad_struct, only: ring_struct
    implicit none
    character(*) xsif_file
    type (ring_struct), target :: ring
    logical, optional :: make_mats6
  end subroutine
end interface

interface
  subroutine bbi_kick_matrix (ele, orb, s_pos, mat6)
    use bmad_struct, only: ele_struct, rp, coord_struct
    implicit none
    type (ele_struct) ele
    type (coord_struct) orb
    real(rp) s_pos
    real(rp) mat6(6,6)
  end subroutine
end interface

end module
