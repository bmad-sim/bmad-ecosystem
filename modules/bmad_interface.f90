!+
! This file defines the interfaces for the BMAD subroutines
!-
!$Id$
!$Log$
!Revision 1.8  2002/07/16 20:44:19  dcs
!*** empty log message ***
!
!Revision 1.7  2002/06/13 14:54:59  dcs
!Interfaced with FPP/PTC
!
!Revision 1.6  2002/02/23 20:32:31  dcs
!Double/Single Real toggle added
!
!Revision 1.5  2002/01/08 21:45:22  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.4  2001/11/29 19:40:11  helms
!Updates from DCS including (*) -> (:)
!
!Revision 1.3  2001/10/12 20:53:50  rwh24
!DCS changes
!
!Revision 1.2  2001/09/27 18:32:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

module bmad_interface

  interface assignment (=)

    subroutine equal_ele_ele (ele1, ele2)
      use bmad_struct
      implicit none
      type (ele_struct), intent(out) :: ele1
      type (ele_struct), intent(in) :: ele2
    end subroutine

    subroutine equal_ring_ring (ring1, ring2)
      use bmad_struct
      implicit none
      type (ring_struct), intent(out) :: ring1
      type (ring_struct), intent(in) :: ring2
    end subroutine

    subroutine equal_coord_coord (coord1, coord2)
      use bmad_struct
      implicit none
      type (coord_struct), intent(out) :: coord1
      type (coord_struct), intent(in) :: coord2
    end subroutine

  end interface


!-------------------

  interface
    subroutine init_taylor (bmad_taylor)
      use bmad_struct
      type (taylor_struct) bmad_taylor(:)
    end subroutine
  end interface

  interface
    subroutine kill_taylor (bmad_taylor)
      use bmad_struct
      type (taylor_struct) bmad_taylor(:)
    end subroutine
  end interface

  interface
    subroutine kill_gen_field (gen_field)
      use bmad_struct
      type (genfield), pointer :: gen_field
    end subroutine
  end interface

  interface
    subroutine concat_taylor (taylor1, taylor2, taylor_out)
      use bmad_struct
      type (taylor_struct) taylor1(:), taylor2(:), taylor_out(:)
    end subroutine
  end interface

  interface
    subroutine accel_sol_mat_calc (ls, c_m, c_e, gamma_old, gamma_new, b_x,  &
        b_y, coord, mat4, vec_st)
      use bmad_struct
      implicit none
      type (coord_struct) coord
      real(rdef) ls
      real(rdef) c_m
      real(rdef) c_e
      real(rdef) gamma_old
      real(rdef) gamma_new
      real(rdef) b_x
      real(rdef) b_y
      real(rdef) mat4(4,4)
      real(rdef) vec_st(4)
    end subroutine
  end interface

  interface
    subroutine add_superimpose (ring, super_ele, ix_super)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (ele_struct) super_ele
      integer ix_super
    end subroutine
  end interface

  interface
    subroutine attribute_bookkeeper (ele, param)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (param_struct) param
    end subroutine
  end interface
 
  interface
    function attribute_index (ele, name)
      use bmad_struct
      implicit none
      integer attribute_index
      type (ele_struct) ele
      character*(*) name
    end function
  end interface

  interface
    function attribute_name(ele, index)
      use bmad_struct
      implicit none
      character*16 attribute_name
      type (ele_struct) ele
      integer index
    end function
  end interface

  interface
    subroutine bmad_parser (in_file, ring, make_mats6)
      use bmad_struct
      implicit none
      character*(*) in_file
      type (ring_struct) ring
      logical, optional :: make_mats6
    end subroutine
  end interface

  interface
    subroutine bmad_parser2 (in_file, ring, orbit_, make_mats6)
      use bmad_struct
      implicit none
      character*(*) in_file
      type (ring_struct) ring
      type (coord_struct), optional :: orbit_(0:n_ele_maxx)
      logical, optional :: make_mats6
    end subroutine
  end interface

  interface
    subroutine bmad_to_cesr (ring, cesr)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (cesr_struct) cesr
    end subroutine
  end interface

  interface
    subroutine bmad_to_db (ring, db)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (db_struct) db
    end subroutine
  end interface

  interface
    subroutine c_to_cbar (ele, cbar_mat)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      real(rdef) cbar_mat(2,2)
    end subroutine
  end interface

  interface
    subroutine calc_z_tune (ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine cesr_crossings(i_train, j_car, species, n_trains_tot, n_cars, &
                                 cross_positions, n_car_spacing, train_spacing)
			use bmad_struct
			implicit none
      integer, intent(in) :: i_train
      integer, intent(in) :: j_car
      integer, intent(in) :: species
      integer, intent(in) :: n_trains_tot
      integer, intent(in) :: n_cars
      integer, optional, intent(in) :: train_spacing(1:10)
      integer, optional, intent(in) :: n_car_spacing(1:10)
      real(rdef), dimension(:), intent(out) :: cross_positions
    end subroutine
  end interface

  interface
    subroutine cesr_elements_get (name, n_found, ele)
      use bmad_struct
      use precision_def
      implicit none
      type (ring_master_struct) :: ele(:)
      integer n_found
      character name*(*)
    end subroutine
  end interface
 
  interface
    subroutine change_basis (coord, ref_energy, ref_z, to_cart, time_disp)
      use bmad_struct
      implicit none
      type (coord_struct) coord
      real(rdef) ref_energy
      real(rdef) ref_z
      real(rdef) time_disp
      logical to_cart
    end subroutine
  end interface

  interface
    Subroutine check_ele_attribute_set (ring, i_ele, attrib_name, &
                                          ix_attrib, err_flag, err_print_flag)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
      integer i_ele
      integer ix_attrib
      character*(*) attrib_name
      logical err_print_flag
      logical err_flag
    end subroutine
  end interface

  interface
    subroutine check_ring_controls (ring, exit_on_error)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
      logical exit_on_error
    end subroutine
  end interface

  interface
    subroutine choose_cesr_lattice (lattice, lat_file, current_lat, &
                                                                ring, choice)
      use bmad_struct
      implicit none
      type (ring_struct), optional :: ring
      character(len=*), optional :: choice
      character*(*) lat_file
      character*40 lattice
      character*40 current_lat
    end subroutine
  end interface

  interface
    subroutine chrom_calc (ring, delta_e, chrom_x, chrom_y)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      real(rdef) delta_e
      real(rdef) chrom_x
      real(rdef) chrom_y
    end subroutine
  end interface

  interface
    subroutine closed_orbit_at_start (ring, co, i_dim, iterate)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct) co
      integer i_dim
      logical iterate
    end subroutine
  end interface

 interface
   subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, &
                                                 eps_rel, eps_abs, init_guess)
     use bmad_struct
     type (ring_struct) ring
     type (coord_struct) closed_orb_(0:n_ele_maxx)
     type (coord_struct), optional :: init_guess
     real(rdef), intent(in) :: eps_rel(:), eps_abs(:)
     integer i_dim
   end subroutine
 end interface

  interface
    subroutine compress_ring (ring, ok)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
      logical ok
    end subroutine
  end interface

  interface
    recursive subroutine control_bookkeeper (ring, ix_ele)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer ix_ele
    end subroutine
  end interface

  interface
    subroutine convert_coords (in_type_str, coord_in, ele, out_type_str, &
                                                                   coord_out)
      use bmad_struct
      implicit none
      character*(*) in_type_str
      character*(*) out_type_str
      type (coord_struct) coord_in
      type (coord_struct) coord_out
      type (ele_struct) ele
    end subroutine
  end interface

  interface
    subroutine create_group (ring, ix_ele, n_control, control_)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (control_struct) control_(:)
      integer ix_ele
      integer n_control
    end subroutine
  end interface

  interface
    subroutine create_overlay (ring, ix_overlay, ix_value, n_slave, con_)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (control_struct) con_(:)
      integer ix_overlay
      integer n_slave
      integer ix_value
    end subroutine
  end interface

  interface
    subroutine create_vsp_volt_elements (ring, ele_type)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer ele_type
    end subroutine
  end interface

  interface
    subroutine custom_emitt_calc (ring, ir, i2, i3, i5a, i5b)
      use bmad_struct
      type (ring_struct) ring
      integer ir
      real(rdef) i2, i3, i5a, i5b
    end subroutine
  end interface

  interface
    subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, ring, db, &
                                                con_, n_con, ok, type_err)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (db_struct) db
      type (control_struct) con_(:)
      integer n_con
      integer ing_num
      integer biggrp_set
      character*12 ing_name
      logical ok, type_err
    end subroutine
  end interface

  interface
    subroutine db_group_to_bmad_group (group_name, group_num, i_biggrp, &
                                           ring, db, ix_ele, ok, type_err)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (db_struct) db
      integer group_num
      integer ix_ele
      integer i_biggrp
      character*12 group_name
      logical ok
      logical type_err
    end subroutine
  end interface

  interface
    subroutine deallocate_ele_pointers (ele)
      use bmad_struct
      implicit none
      type (ele_struct) ele
    end subroutine
  end interface

  interface
    subroutine do_mode_flip (ele, ele_flip)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (ele_struct) ele_flip
    end subroutine
  end interface

  interface
    subroutine ele_to_taylor (ele, orb0, param)
      use bmad_struct
      implicit none
      type (ele_struct), intent(inout) :: ele
      type (coord_struct), optional, intent(in) :: orb0
      type (param_struct), optional, intent(in) :: param
    end subroutine
  end interface

  interface
    subroutine element_locator (ele_name, ring, ix_ele)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer ix_ele
      character*(*) ele_name
    end subroutine
  end interface

  interface
    subroutine emitt_calc (ring, what, mode)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (modes_struct) mode
      integer what
    end subroutine
  end interface

  interface
    real(rdef) function field_interpolate_3d (position, field_mesh, &
                                                         deltas, position0)
      use precision_def
      implicit none
      real(rdef), intent(in) :: position(3), deltas(3)
      real(rdef), intent(in) :: field_mesh(:,:,:)
      real(rdef), intent(in), optional :: position0(3)
    end function
  end interface

  interface
    subroutine find_element_ends (ring, ix_ele, ix_start, ix_end)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer ix_ele
      integer ix_start
      integer ix_end
    end subroutine
  end interface

  interface
    subroutine fitpoly(coe, x, y, order, samples)
      use precision_def
      implicit none
      integer order
      integer samples
      real(rdef) coe(0:)
      real(rdef) x(:)
      real(rdef) y(:)
    end subroutine
  end interface

  interface
    subroutine get_lattice_list (lat_list, num_lats, directory)
      use precision_def
      implicit none
      integer num_lats
      character*(*) directory
      character*40 lat_list(*)
    end subroutine
  end interface

  interface
      function hypergeom(hgcx, arg)
      use precision_def
      implicit none
      integer hgcx
      real(rdef) arg
      real(rdef) hypergeom
    end function
  end interface

  interface
    subroutine identify_db_node (db_name, db, db_ptr, ok, type_err)
      use bmad_struct
      implicit none
      type (db_struct), target :: db
      type (db_element_struct), pointer :: db_ptr(:)
      character*(*) db_name
      logical ok
      logical type_err
    end subroutine
  end interface

  interface
    subroutine init_ele (ele)
      use bmad_struct
      implicit none
      type (ele_struct) ele
    end subroutine
  end interface

  interface
    subroutine init_LRBBI(ring, oppos_ring, LRBBI_ele, ix_LRBBI, ix_oppos)
      use bmad_struct
      implicit none
      type (ring_struct) ring
			type (ring_struct) oppos_ring
      type (ele_struct) LRBBI_ele
			integer, intent(in) :: ix_LRBBI, ix_oppos
    end subroutine
  end interface

  interface
    subroutine insert_element (ring, insert_ele, insert_index)
      use bmad_struct
      implicit none
      type (ring_struct) ring
			type (ring_struct) oppos_ring
      type (ele_struct) insert_ele
      integer insert_index
    end subroutine
  end interface

  interface
    subroutine insert_LRBBI (ring, oppos_ring, cross_positions, ix_LRBBI)
			use bmad_struct
      type (ring_struct) ring
			type (ring_struct) oppos_ring 
     real(rdef), dimension(:), intent(inout) :: cross_positions
      integer, dimension(:), intent(inout) :: ix_LRBBI
    end subroutine
  end interface

  interface
    subroutine k_to_quad_calib(k_theory, energy, cu_theory, k_base,  &
                                                     dk_gev_dcu, cu_per_k_gev)
      use precision_def
      implicit none
      real(rdef) energy
      real(rdef) k_theory(0:*)
      real(rdef) k_base(0:120)
      real(rdef) cu_per_k_gev(0:120)
      real(rdef) dk_gev_dcu(0:*)
      integer cu_theory(0:*)
    end subroutine
  end interface

  interface
    subroutine lattice_to_bmad_file_name (lattice, bmad_file_name)
      implicit none
      character*(*) lattice
      character*(*) bmad_file_name
    end subroutine
  end interface
 
  interface
    subroutine LRBBI_crossings (n_bucket, oppos_buckets, cross_positions)
      use precision_def
      implicit none
			real(rdef), intent(in) :: n_bucket
      real(rdef), dimension(:), intent(in) :: oppos_buckets
      real(rdef), dimension(:), intent(inout) :: cross_positions
    end subroutine LRBBI_crossings 
  end interface
 
  interface
    subroutine make_g_mats (ele, g_mat, g_inv_mat)
      use bmad_struct       
      implicit none
      type (ele_struct) ele
      real(rdef) g_mat(4,4)
      real(rdef) g_inv_mat(4,4)
    end subroutine
  end interface

  interface
    subroutine make_g2_mats (twiss, g2_mat, g2_inv_mat)
      use bmad_struct
      implicit none
      type (twiss_struct) twiss
      real(rdef) g2_mat(2,2)
      real(rdef) g2_inv_mat(2,2)
    end subroutine
  end interface

  interface
    subroutine make_hybrid_ring (r_in, use_ele, remove_markers, &
                                             r_out, ix_out, use_taylor, orb0_)
      use bmad_struct
      implicit none
      type (ring_struct) r_in
      type (ring_struct) r_out
      integer ix_out(:)
      logical remove_markers
      logical use_ele(:)
      logical, optional :: use_taylor
      type (coord_struct), optional :: orb0_(0:)
    end subroutine
  end interface

  interface
    subroutine make_LRBBI(master_ring, master_ring_oppos, ring, &
    													ix_LRBBI, master_ix_LRBBI)
      use bmad_struct
      implicit none
      type (ring_struct), dimension(:) :: ring
      type (ring_struct) :: master_ring
      type (ring_struct) :: master_ring_oppos
      integer, dimension(:,:) :: ix_LRBBI
      integer, dimension(:,:) :: master_ix_LRBBI
    end subroutine
  end interface
 
  interface
    subroutine make_mat6 (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct), optional :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine make_mat6_custom (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine make_mat6_taylor (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine make_mat6_bmad (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine make_mat6_runge_kutta (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine make_mat6_symp_lie (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine make_mat6_tracking (ele, param, c0, c1)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) :: c0, c1
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine custom_radiation_integrals (ring, ir, orb_)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct) orb_(0:n_ele_maxx)
      integer ir
    end subroutine
  end interface

  interface
    subroutine make_v_mats (ele, v_mat, v_inv_mat)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      real(rdef) v_mat(4,4)
      real(rdef) v_inv_mat(4,4)
    end subroutine
  end interface

  interface
    subroutine mark_LRBBI(master_ring, master_ring_oppos, ring, crossings)
      use bmad_struct
      implicit none
      type (ring_struct), dimension(:) :: ring
      type (ring_struct) :: master_ring
      type (ring_struct) :: master_ring_oppos
      real(rdef), dimension(:,:) :: crossings
    end subroutine
  end interface

  interface
    subroutine mat6_dispersion (mat6, e_vec)
      use bmad_struct
      implicit none
      real(rdef), intent(inout) :: mat6(6,6)
      real(rdef), intent(in) :: e_vec(:)
    end subroutine
  end interface    

  interface
    subroutine mat_inverse (mat, mat_inv)
      use precision_def
      implicit none
      real(rdef), intent(in)  :: mat(:,:)
      real(rdef), intent(out) :: mat_inv(:,:)
    end subroutine
  end interface

  interface
    subroutine mat_symp_check (mat, error)
      use precision_def
      implicit none
      real(rdef), intent(in) :: mat(:,:)
      real(rdef) error
    end subroutine
  end interface

  interface
    subroutine mat_symp_decouple(t0, tol, stat, U, V, Ubar, Vbar, G,  &
                                                 twiss1, twiss2, type_out)
      use bmad_struct
      implicit none
      type (twiss_struct) twiss1, twiss2
      real(rdef) t0(4,4), U(4,4), V(4,4), tol
      real(rdef) Ubar(4,4), Vbar(4,4), G(4,4)
      integer stat
      logical type_out
    end subroutine
  end interface

  interface
    subroutine mat_symplectify (mat_in, mat_symp)
      use precision_def
      real(rdef), intent(in)  :: mat_in(:,:)
      real(rdef), intent(out) :: mat_symp(:,:)
    end subroutine
  end interface

  interface
    subroutine mobius_twiss_calc (ele, v_mat)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      real(rdef) v_mat(4,4)
    end subroutine
  end interface

  interface
    subroutine multipole_ab_to_kt (an, bn, knl, tn)
      use bmad_struct
      implicit none
      real(rdef) an(0:), bn(0:)
      real(rdef) knl(0:), tn(0:)
    end subroutine
  end interface

  interface
    function c_multi (n, m)
      use precision_def
      implicit none
      real(rdef) c_multi
      integer, intent(in) :: n, m
    end function
  end interface

  interface
    subroutine multipole_ele_to_ab (ele, particle, a, b, use_tilt)
      use bmad_struct
      type (ele_struct) ele
      integer particle
      real(rdef) a(0:), b(0:)
      real(rdef) value(n_attrib_maxx)
      logical use_tilt
    end subroutine
  end interface

  interface  
    subroutine multipole_kick (knl, tilt, n, coord)
      use bmad_struct
      implicit none
      type (coord_struct) coord
      real(rdef) knl
      real(rdef) tilt
      integer n
    end subroutine
  end interface

  interface
    subroutine multipole_kt_to_ab (knl, tn, an, bn)
      use bmad_struct
      implicit none
      real(rdef) an(0:), bn(0:)
      real(rdef) knl(0:), tn(0:)
    end subroutine
  end interface

  interface
    subroutine multipole_ele_to_kt (ele, particle, knl, tilt, use_ele_tilt)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      real(rdef) knl(0:)
      real(rdef) tilt(0:)
      logical use_ele_tilt
      integer particle
    end subroutine
  end interface

  interface
    subroutine name_to_list (ring, ele_names, use_ele)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      logical use_ele(:)
      character*(*) ele_names(:)
    end subroutine
  end interface

  interface
    subroutine new_control (ring, ix_ele)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer ix_ele
    end subroutine
  end interface

  interface
    subroutine offset_particle (ele, param, coord, set, &
                    set_canonical, set_tilt, set_multipoles, set_hvkicks)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) coord
      type (param_struct) param
      logical set
      logical, optional :: set_canonical, set_multipoles, set_tilt, set_hvkicks
    end subroutine
  end interface

  interface
    subroutine one_turn_matrix (ring, mat6)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      real(rdef) mat6(6,6)
    end subroutine
  end interface

  interface
    subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)
      use bmad_struct
      type (ele_struct) ele
      real(rdef) phi_a
      real(rdef) phi_b
      real(rdef) mat4(4,4)
    end subroutine
  end interface

  interface
    subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele, &
                                                      stable, growth_rate, chi)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: track(:)
      type (coord_struct), intent(out) :: track0
      type (ele_struct), intent(out) :: ele
      real(rdef), intent(out) :: growth_rate, chi
      integer, intent(in) :: i_dim
      logical, intent(out) :: stable
    end subroutine
  end interface

  interface
    subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in), target :: track(:)
      type (coord_struct), intent(out) :: track0
      real(rdef), intent(out) :: mat1(:,:)
      real(rdef), intent(out) :: chi
      integer, intent(in) :: i_dim
    end subroutine
  end interface

  interface
    subroutine order_super_lord_slaves (ring, ix_lord)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
      integer ix_lord
    end subroutine
  end interface

  interface
    subroutine phase_space_fit (x, xp, twiss, tune, emitt, x_0, xp_0, chi, tol)
      use bmad_struct
      implicit none
      type (twiss_struct) twiss
      real(rdef), optional :: tol
      real(rdef) x(:), xp(:)
      real(rdef) tune, emitt
      real(rdef) x_0, xp_0, chi
    end subroutine
  end interface

  interface
    subroutine quad_beta_ave (ring, ix_ele, beta_x_ave, beta_y_ave)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer ix_ele
      real(rdef) beta_x_ave
      real(rdef) beta_y_ave
    end subroutine
  end interface

  interface
    subroutine quad_calib (lattice, k_theory, k_base,  &
                     len_quad, cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)
      use precision_def
      implicit none
      character lattice*(*)
      real(rdef) k_theory(0:*)
      real(rdef) k_base(0:*)
      real(rdef) len_quad(0:*)
      real(rdef) cu_per_k_gev(0:*)
      real(rdef) dk_gev_dcu(0:*)
      real(rdef) quad_rot(0:*)
      integer cu_theory(0:*)
    end subroutine
  end interface

  interface
    subroutine radiation_integrals (ring, orb_, mode)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
      type (coord_struct), target :: orb_(0:)
      type (modes_struct) mode
    end subroutine
  end interface

  interface
    subroutine read_butns_file (butns_num, butns, db, ok, type_err)
      use bmad_struct
      implicit none
      type (db_struct) db
      type (butns_struct) butns
      integer butns_num
      logical ok, type_err
    end subroutine
  end interface

  interface
    subroutine read_digested_bmad_file (in_file_name, ring, version)
      use bmad_struct
      implicit none
      type (ring_struct), intent(out) :: ring
      integer version
      character*(*) in_file_name
    end subroutine
  end interface

  interface
    function relative_mode_flip (ele1, ele2)
      use bmad_struct
      implicit none
      logical relative_mode_flip
      type (ele_struct) ele1
      type (ele_struct) ele2
    end function
  end interface

  interface
    subroutine ring_geometry (ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    recursive subroutine ring_make_mat6 (ring, ix_ele, coord_)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct), optional :: coord_(0:n_ele_maxx)
      integer ix_ele
    end subroutine
  end interface

  interface
    subroutine set_ele_attribute (ring, i_ele, attrib_name, &
                                attrib_value, err_flag, make_mat6_flag, orbit_)
      use bmad_struct
      implicit none
      type (ring_struct) :: ring
      type (coord_struct), optional :: orbit_(0:n_ele_maxx)
      real(rdef) attrib_value
      integer i_ele
      character*(*) attrib_name
      logical make_mat6_flag
      logical err_flag
    end subroutine
  end interface

  interface
    subroutine ring_to_quad_calib (ring, cesr, k_theory, k_base,  &
                     len_quad, cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)
      use bmad_struct
      implicit none
      type (cesr_struct)  cesr
      type (ring_struct)  ring
      real(rdef) k_theory(0:*)
      real(rdef) k_base(0:*)
      real(rdef) len_quad(0:*)
      real(rdef) cu_per_k_gev(0:*)
      real(rdef) dk_gev_dcu(0:*)
      real(rdef) quad_rot(0:*)
      integer cu_theory(0:*)
    end subroutine
  end interface

  interface
    subroutine s_calc (ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine set_on (key, ring, on_switch, orb_)
      use bmad_struct
      type (ring_struct) ring
      type (coord_struct), optional :: orb_(0:)
      integer key
      logical on_switch
    end subroutine
  end interface

 interface
    subroutine set_ptc (param, taylor_order, integ_order, &
                                      num_steps, no_cavity, exact_calc)
      use bmad_struct
      implicit none
      type (param_struct), optional :: param
      integer, optional :: taylor_order
      integer, optional :: integ_order
      integer, optional :: num_steps
      logical, optional :: no_cavity
      logical, optional :: exact_calc
    end subroutine
  end interface

  interface
    subroutine set_symmetry (symmetry, ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      integer symmetry
    end subroutine
  end interface

  interface
    subroutine set_taylor_order (order, override_flag)
      use bmad_struct
      implicit none
        integer, intent(in) :: order
        logical, optional, intent(in) :: override_flag
    end subroutine
  end interface

  interface
    subroutine set_tune (phi_x_set, phi_y_set, dk1, ring, orb_, ok)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct) orb_(0:)
      real(rdef) phi_x_set
      real(rdef) phi_y_set
      real(rdef) dk1(:)
      logical ok
    end subroutine
  end interface

  interface
    subroutine set_z_tune (ring)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
    end subroutine
  end interface

  interface
    subroutine sort_taylor_terms (taylor_in, taylor_sorted)
      use bmad_struct
      implicit none
      type (taylor_struct), intent(in)  :: taylor_in
      type (taylor_struct) :: taylor_sorted
    end subroutine
  end interface

  interface
    subroutine split_ring (ring, s_split, ix_split, split_done)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      real(rdef) s_split
      integer ix_split
      logical split_done
    end subroutine
  end interface

  interface
    subroutine taylor_to_mat6 (bmad_taylor, c0, mat6, c1)
      use bmad_struct
      implicit none
      type (taylor_struct), target, intent(in) :: bmad_taylor(6)
      type (coord_struct), intent(in) :: c0
      type (coord_struct), intent(out) :: c1
      real(rdef), intent(out) :: mat6(6,6)
    end subroutine
  end interface

  interface 
    subroutine taylor_propagate1 (bmad_taylor, ele, param)
      use bmad_struct
      implicit none
      type (taylor_struct) bmad_taylor(:)
      type (ele_struct) ele
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine tilt_coords (tilt_val, coord, set)
      use precision_def
      implicit none
      real(rdef) tilt_val
      real(rdef) coord(:)
      logical set
    end subroutine
  end interface

  interface
    subroutine track_all (ring, orbit_)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct) orbit_(0:)
    end subroutine
  end interface

  interface
    subroutine track_bend (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct) start
      type (coord_struct) end
      type (ele_struct) ele
      type (param_struct) param
      logical is_lost
    end subroutine
  end interface

  interface
    subroutine track_long (ring, orbit_, ix_start, direction, mats627)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct) orbit_(0:*)
      type (mat627_struct) mats627(*)
      integer ix_start
      integer direction
    end subroutine
  end interface
 
  interface
    subroutine track_many (ring, orbit_, ix_start, ix_end, direction)
      use bmad_struct
      implicit none
      type (ring_struct)  ring
      type (coord_struct)  orbit_(0:)
      integer ix_start
      integer ix_end
      integer direction
    end subroutine
  end interface

  interface
    subroutine track_taylor (start, bmad_taylor, end)
      use bmad_struct
      implicit none
      type (taylor_struct), intent(in) :: bmad_taylor(6)
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
    end subroutine
  end interface

  interface
    subroutine transfer_ele_pointers (ele1, ele2)
      use bmad_struct
      implicit none
      type (ele_struct), intent(in)  :: ele2
      type (ele_struct), intent(out) :: ele1
    end subroutine
  end interface

  interface
    subroutine transfer_mat_from_tracking (ele, param, start, d_orb, end, error)
      use bmad_struct
      implicit none
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
      type (coord_struct), intent(in) :: start
      type (coord_struct), optional, intent(out) :: end
      type (coord_struct), optional :: d_orb
      real(rdef), optional, intent(out) :: error
    end subroutine
  end interface

  interface
    subroutine transfer_mat_from_twiss (twiss1, twiss2, mat)
      use bmad_struct
      implicit none
      type (twiss_struct) twiss1
      type (twiss_struct) twiss2
      real(rdef) mat(2,2)
    end subroutine
  end interface

  interface
    subroutine track1 (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_runge_kutta (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_linear (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_taylor (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_custom (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_bmad (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_symp_lie (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_symp_map (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct), intent(in) :: start
      type (coord_struct), intent(out) :: end
      type (ele_struct), intent(inout) :: ele
      type (param_struct), intent(inout) :: param
    end subroutine
  end interface

  interface
    subroutine track1_wiedemann_wiggler (start, ele, param, end)
      use bmad_struct
      implicit none
      type (coord_struct) start
      type (coord_struct) end
      type (param_struct) param
      type (ele_struct) ele
      logical is_lost
    end subroutine
  end interface

  interface
    subroutine twiss_and_track (ring, orb)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (coord_struct) orb(0:)
    end subroutine
  end interface

  interface
    subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
                                                                   start, end)
      use bmad_struct
      implicit none
      type (ele_struct), optional :: ele3
      type (ele_struct) ele1
      type (ele_struct) ele2
      type (coord_struct), optional :: start
      type (coord_struct), optional :: end
      type (param_struct) param
      real(rdef) del_s
    end subroutine
  end interface

  interface
    subroutine twiss_and_track_body (ele1, ele2, param, del_s, ele3, &
                                                                   start, end)
      use bmad_struct
      implicit none
      type (ele_struct), optional :: ele3
      type (ele_struct) ele1
      type (ele_struct) ele2
      type (coord_struct), optional :: start
      type (coord_struct), optional :: end
      type (param_struct) param
      real(rdef) del_s
    end subroutine
  end interface

  interface
    subroutine twiss_at_element (ring, ix_ele, start, end, average)
      use bmad_struct
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
      use bmad_struct
      implicit none
      type (ring_struct) :: ring
      type (ele_struct) :: ele
      real(rdef) s
    end subroutine
  end interface

  interface
    subroutine twiss_at_start (ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine twiss_from_mat6 (mat6, ele, stable, growth_rate)
      use bmad_struct
      implicit none
      type (ele_struct), intent(out) :: ele
      real(rdef), intent(in) :: mat6(6,6)
      real(rdef), intent(out) :: growth_rate
      logical, intent(out) :: stable
    end subroutine
  end interface

  interface
    subroutine twiss_from_tracking (ring, closed_orb_, d_orb, error)
      use bmad_struct
      type (ring_struct), intent(inout) :: ring
      type (coord_struct), intent(in) :: closed_orb_(0:n_ele_maxx)
      type (coord_struct), intent(in) :: d_orb
      real(rdef), intent(out) :: error
    end subroutine
  end interface

  interface
    subroutine twiss_propagate1 (ele1, ele2)
      use bmad_struct
      implicit none
      type (ele_struct) ele1
      type (ele_struct) ele2
    end subroutine
  end interface

  interface
    subroutine twiss_propagate_all (ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine type_coord (coord)
      use bmad_struct
      implicit none
      type (coord_struct) coord
    end subroutine
  end interface

  interface
    subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, &
                                      twiss_type, type_control, ring)
      use bmad_struct
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
      use bmad_struct
      implicit none
      type (ele_struct) ele
      integer frequency_units
    end subroutine
  end interface
 
  interface
    subroutine type2_ele (ele, type_zero_attrib, type_mat6, type_taylor, &
                            twiss_type, type_control, lines, n_lines, ring)
      use bmad_struct
      implicit none
      type (ele_struct), intent(in) :: ele
      type (ring_struct), optional, intent(in) :: ring
      integer, intent(in) :: type_mat6
      integer, intent(out) :: n_lines
      integer, intent(in) :: twiss_type
      logical, intent(in) :: type_control, type_taylor
      logical, intent(in) :: type_zero_attrib
      character*80, pointer :: lines(:)
    end subroutine
  end interface

  interface
    subroutine type_taylors (bmad_taylor)
      use bmad_struct
      implicit none
      type (taylor_struct) bmad_taylor(:)
    end subroutine
  end interface

  interface
    subroutine type2_taylors (bmad_taylor, lines, n_lines)
      use bmad_struct
      implicit none
      type (taylor_struct), intent(in) :: bmad_taylor(6)
      integer, intent(out) :: n_lines
      character*80, pointer :: lines(:)
    end subroutine
  end interface

  interface
    subroutine type2_twiss (ele, frequency_units, lines, n_lines)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      integer frequency_units
      integer n_lines
      character*(*) lines(:)
    end subroutine
  end interface
 
  interface
    subroutine write_digested_bmad_file (digested_name, ring,  &
                                                      n_files, file_names)
      use bmad_struct
      implicit none
      type (ring_struct), intent(in) :: ring
      integer n_files
      character*(*) digested_name
      character*(*) file_names(*)
    end subroutine
  end interface

  interface
    subroutine update_hybrid_list (ring, n_in, use_ele)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      logical use_ele(:)
      integer n_in
    end subroutine
  end interface

  interface
    subroutine adjust_control_struct (ring, ix_ele)
      use bmad_struct
      implicit none
      type (ring_struct), target :: ring
      integer ix_ele
    end subroutine
  end interface

  interface
    subroutine make_mat627 (ele, param, direction, mat627)
      use bmad_struct
      implicit none
      type (ele_struct), target :: ele
      type (param_struct) param
      real(rdef) mat627(6,27)
      integer direction
    end subroutine
  end interface

  interface
    recursive subroutine ring_make_mat627 (ring, ix_ele, direction, mats627)
      use bmad_struct
      implicit none
      type (ring_struct) ring
      type (mat627_struct) mats627(:)
      integer direction
      integer ix_ele
    end subroutine
  end interface

  interface
    subroutine set_design_linear (ring)
      use bmad_struct
      implicit none
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine track1_627 (start, ele, param, mat627, end)
      use bmad_struct
      implicit none
      type (coord_struct) start
      type (coord_struct) end
      type (ele_struct) ele
      type (param_struct) param
      real(rdef) mat627(6,27)
    end subroutine
  end interface

  interface
    subroutine twiss_from_mat2 (mat, det, twiss, stat, tol, type_out)
      use bmad_struct
      implicit none
      type (twiss_struct) twiss
      integer psize
      integer stat
      real(rdef) mat(:, :)
      real(rdef) det
      real(rdef) tol
      logical type_out
    end subroutine
  end interface

  interface
    subroutine twiss_to_1_turn_mat (twiss, phi, mat2)
      use bmad_struct
      type (twiss_struct) twiss
      real(rdef) phi
      real(rdef) mat2(2,2)
    end subroutine
  end interface

  interface
    subroutine bbi_kick_matrix (ele, orb, s_pos, mat6)
      use bmad_struct
      implicit none
      type (ele_struct) ele
      type (coord_struct) orb
      real(rdef) s_pos
      real(rdef) mat6(6,6)
    end subroutine
  end interface

end module
