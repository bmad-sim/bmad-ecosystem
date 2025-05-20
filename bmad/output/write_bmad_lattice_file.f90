!+ 
! Subroutine write_bmad_lattice_file (bmad_file, lat, err, output_form, orbit0)
!
! Subroutine to write a Bmad lattice file using the information in a lat_struct.
! Optionally only part of the lattice can be generated.
! Also see: write_lattice_in_foreign_format
!
! Note: bmad_com parameters that are changed from their default value are
! saved in the lattice file.
!
! Input:
!   bmad_file     -- Character(*): Name of the output lattice file.
!   lat           -- lat_struct: Holds the lattice information.
!   output_form   -- integer, optional: 
!                       binary$   -> Write grid_field info in binary hdf5 form in separate files. Default.
!                                      All other fields are writen in separate files in ASCII
!                       ascii$    -> Fields will be put in separate ASCII files.
!                       one_file$ -> Everything in one file.
!   orbit0        -- coord_struct, optional: Initial orbit. Used to write the inital orbit if the 
!                       lattice geometry is closed.
!
! Output:
!   err    -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_bmad_lattice_file (bmad_file, lat, err, output_form, orbit0)

use write_lattice_file_mod, dummy => write_bmad_lattice_file
use expression_mod, only: end_stack$, variable$, split_expression_string, expression_stack_to_string

implicit none

type (multipass_region_lat_struct), target :: mult_lat
type (multipass_region_ele_struct), pointer :: mult_ele(:), m_ele

type (lat_struct), target :: lat
type (coord_struct), optional :: orbit0
type (ele_attribute_struct) attrib
type (branch_struct), pointer :: branch, branch2
type (ele_struct), pointer :: ele, super, slave, lord, lord2, s1, s2, multi_lord
type (ele_struct), pointer :: slave1, slave2, ele2, ele_dflt, ele0, girder
type (ele_struct), target :: ele_default(n_key$), this_ele
type (ele_pointer_struct), allocatable :: named_eles(:)  ! List of unique element names 
type (ele_attribute_struct) info
type (control_struct), pointer :: ctl, ctl2
type (control_ramp1_struct), pointer :: rmp
type (taylor_term_struct) tm
type (multipass_all_info_struct), target :: m_info
type (multipass_ele_info_struct), pointer :: e_info
type (wake_lr_struct), pointer :: lr
type (wake_sr_struct), pointer :: sr
type (wake_lr_mode_struct), pointer :: lrm
type (wake_sr_mode_struct), pointer :: srm
type (wake_sr_z_long_struct), pointer :: srz
type (ele_pointer_struct), pointer :: ss1(:), ss2(:)
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (gen_grad_map_struct), pointer :: gg_map
type (grid_field_struct), pointer :: g_field
type (em_taylor_term_struct), pointer :: t_term
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v
type (bmad_common_struct), parameter :: bmad_com_default = bmad_common_struct()
type (space_charge_common_struct), parameter :: space_charge_com_default = space_charge_common_struct()
type (ac_kicker_struct), pointer :: ac
type (expression_atom_struct), pointer :: stack(:)
type (str_index_struct) str_index
type (lat_ele_order_struct) order
type (material_struct), pointer :: material

real(rp) s0, x_lim, y_lim, val, x, y, z, fid, f

character(*) bmad_file
character(4000) line
character(2000) line2
character(200) file_name, path, basename, fname
character(120), allocatable :: list(:)
character(60) alias
character(40) name, look_for, attrib_name
character(40), allocatable :: names(:)
character(16) polar, dependence
character(40) angle
character(4) end_str, last
character(2), parameter :: spin_quat_name(0:3) = ['S1', 'Sx', 'Sy', 'Sz']
character(*), parameter :: r_name = 'write_bmad_lattice_file'

integer, optional :: output_form
integer i, j, k, n, ii, ix, iu, im, ix_ptr, iu2, iuw, ios, ixs, ie1, ie, ib, ib1, ic
integer unit(6), n_names, ix_match, ie2, id1, id2, id3, j1, j2, ip, it
integer ix_slave, ix_ss, ix_l, ix_r, ix_pass
integer ix_lord, ix_super, default_val, imax, ibr
integer, allocatable :: an_indexx(:), index_list(:)

logical, optional :: err
logical unit_found, write_term, found, in_multi_region, have_expand_lattice_line
logical x_lim_good, y_lim_good, is_default, need_new_region, err_flag, has_been_added

! Init...
! Init default parameters.

do i = 1, size(ele_default)
  call init_ele (ele_default(i), i)
  call deallocate_ele_pointers(ele_default(i))   ! don't need.
enddo

if (present(err)) err = .true.
call ele_order_calc(lat, order)

! Open the file

iu = lunget()
call fullfilename (bmad_file, file_name)
open (iu, file = file_name, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(bmad_file))
  return
endif

ix = splitfilename(file_name, path, basename)
if (path == '') path = '.'

! Custom attribute names

call custom_ele_attrib_name_list(index_list, list)
do i = 1, size(index_list)
  write (iu, '(a, i0, 2a)') 'parameter[custom_attribute', index_list(i), '] = ', trim(list(i))
enddo

! Global custom attribute values

if (allocated(lat%custom)) then
  write (iu, '(a)') 
  do i = 1, size(lat%custom)
    name = attribute_name(def_parameter$, i+custom_attribute0$)
    if (name(1:1) == '!') cycle
    write (iu, '(4a)') 'parameter[', trim(name), '] = ', re_str(lat%custom(i))
  enddo
endif

! Non-elemental stuff

if (lat%title /= ' ')            write (iu, '(4a)')    'title, "', trim(lat%title), '"'
if (lat%lattice /= ' ')          write (iu, '(4a)')    'parameter[lattice]     = "', trim(lat%lattice), '"'

write (iu, '(4a)') 'parameter[geometry] = ', geometry_name(lat%param%geometry)
if (.not. lat%param%live_branch) write (iu, '(a)') 'parameter[live_branch] = F'
if (lat%input_taylor_order /= 0) write (iu, '(a, i0)') 'parameter[taylor_order] = ', lat%input_taylor_order

write (iu, '(a)')
write (iu, '(4a)')    'parameter[p0c]                      = ', re_str(lat%ele(0)%value(p0c_start$))
write (iu, '(4a)')    'parameter[particle]                 = ', trim(species_name(lat%param%particle))
if (lat%param%default_tracking_species /= ref_particle$) then
  write (iu, '(4a)')    'parameter[default_tracking_species] = ', species_name(lat%param%default_tracking_species)
endif

if (lat%param%n_part /= 0)             write (iu, '(2a)') 'parameter[n_part]                 = ', re_str(lat%param%n_part)

ele => lat%ele(lat%n_ele_track)
if (ele%name /= 'END' .or. ele%key /= marker$) then
  write (iu, '(a)') 'parameter[no_end_marker]          =  T'
endif

if (lat%photon_type /= incoherent$) then
  write (iu, '(3a)') 'parameter[photon_type] = ', photon_type_name(lat%photon_type)
endif

! Bmad_com

call write_if_real_param_changed  (bmad_com%max_aperture_limit,              bmad_com_default%max_aperture_limit,               'bmad_com[max_aperture_limit]')
call write_if_real_param_changed  (bmad_com%default_ds_step,                 bmad_com_default%default_ds_step,                  'bmad_com[default_ds_step]')
call write_if_real_param_changed  (bmad_com%significant_length,              bmad_com_default%significant_length,               'bmad_com[significant_length]')
call write_if_real_param_changed  (bmad_com%rel_tol_tracking,                bmad_com_default%rel_tol_tracking,                 'bmad_com[rel_tol_tracking]')
call write_if_real_param_changed  (bmad_com%abs_tol_tracking,                bmad_com_default%abs_tol_tracking,                 'bmad_com[abs_tol_tracking]')
call write_if_real_param_changed  (bmad_com%rel_tol_adaptive_tracking,       bmad_com_default%rel_tol_adaptive_tracking,        'bmad_com[rel_tol_adaptive_tracking]')
call write_if_real_param_changed  (bmad_com%abs_tol_adaptive_tracking,       bmad_com_default%abs_tol_adaptive_tracking,        'bmad_com[abs_tol_adaptive_tracking]')
call write_if_real_param_changed  (bmad_com%init_ds_adaptive_tracking,       bmad_com_default%init_ds_adaptive_tracking,        'bmad_com[init_ds_adaptive_tracking]')
call write_if_real_param_changed  (bmad_com%min_ds_adaptive_tracking,        bmad_com_default%min_ds_adaptive_tracking,         'bmad_com[min_ds_adaptive_tracking]')
call write_if_real_param_changed  (bmad_com%fatal_ds_adaptive_tracking,      bmad_com_default%fatal_ds_adaptive_tracking,       'bmad_com[fatal_ds_adaptive_tracking]')
call write_if_real_param_changed  (bmad_com%autoscale_amp_abs_tol,           bmad_com_default%autoscale_amp_abs_tol,            'bmad_com[autoscale_amp_abs_tol]')
call write_if_real_param_changed  (bmad_com%autoscale_amp_rel_tol,           bmad_com_default%autoscale_amp_rel_tol,            'bmad_com[autoscale_amp_rel_tol]')
call write_if_real_param_changed  (bmad_com%autoscale_phase_tol,             bmad_com_default%autoscale_phase_tol,              'bmad_com[autoscale_phase_tol]')
call write_if_real_param_changed  (bmad_com%electric_dipole_moment,          bmad_com_default%electric_dipole_moment,           'bmad_com[electric_dipole_moment]')
call write_if_real_param_changed  (bmad_com%synch_rad_scale,                 bmad_com_default%synch_rad_scale,                  'bmad_com[synch_rad_scale]')
call write_if_real_param_changed  (bmad_com%sad_eps_scale,                   bmad_com_default%sad_eps_scale,                    'bmad_com[sad_eps_scale]')
call write_if_real_param_changed  (bmad_com%sad_amp_max,                     bmad_com_default%sad_amp_max,                      'bmad_com[sad_amp_max]')
call write_if_int_param_changed   (bmad_com%sad_n_div_max,                   bmad_com_default%sad_n_div_max,                    'bmad_com[sad_n_div_max]')
call write_if_int_param_changed   (bmad_com%taylor_order,                    bmad_com_default%taylor_order,                     'bmad_com[taylor_order]')
call write_if_int_param_changed   (bmad_com%runge_kutta_order,               bmad_com_default%runge_kutta_order,                'bmad_com[runge_kutta_order]')
call write_if_int_param_changed   (bmad_com%default_integ_order,             bmad_com_default%default_integ_order,              'bmad_com[default_integ_order]')
call write_if_int_param_changed   (bmad_com%max_num_runge_kutta_step,        bmad_com_default%max_num_runge_kutta_step,         'bmad_com[max_num_runge_kutta_step]')
call write_if_logic_param_changed (bmad_com%rf_phase_below_transition_ref,   bmad_com_default%rf_phase_below_transition_ref,    'bmad_com[rf_phase_below_transition_ref]')
call write_if_logic_param_changed (bmad_com%sr_wakes_on,                     bmad_com_default%sr_wakes_on,                      'bmad_com[sr_wakes_on]')
call write_if_logic_param_changed (bmad_com%lr_wakes_on,                     bmad_com_default%lr_wakes_on,                      'bmad_com[lr_wakes_on]')
call write_if_logic_param_changed (bmad_com%high_energy_space_charge_on,     bmad_com_default%high_energy_space_charge_on,      'bmad_com[high_energy_space_charge_on]')
call write_if_logic_param_changed (bmad_com%csr_and_space_charge_on,         bmad_com_default%csr_and_space_charge_on,          'bmad_com[csr_and_space_charge_on]')
call write_if_logic_param_changed (bmad_com%spin_tracking_on,                bmad_com_default%spin_tracking_on,                 'bmad_com[spin_tracking_on]')
call write_if_logic_param_changed (bmad_com%spin_sokolov_ternov_flipping_on, bmad_com_default%spin_sokolov_ternov_flipping_on,  'bmad_com[spin_sokolov_ternov_flipping_on]')
call write_if_logic_param_changed (bmad_com%radiation_damping_on,            bmad_com_default%radiation_damping_on,             'bmad_com[radiation_damping_on]')
call write_if_logic_param_changed (bmad_com%radiation_zero_average,          bmad_com_default%radiation_zero_average,           'bmad_com[radiation_zero_average]')
call write_if_logic_param_changed (bmad_com%radiation_fluctuations_on,       bmad_com_default%radiation_fluctuations_on,        'bmad_com[radiation_fluctuations_on]')
call write_if_logic_param_changed (bmad_com%conserve_taylor_maps,            bmad_com_default%conserve_taylor_maps,             'bmad_com[conserve_taylor_maps]')
call write_if_logic_param_changed (bmad_com%absolute_time_ref_shift,         bmad_com_default%absolute_time_ref_shift,          'bmad_com[absolute_time_ref_shift]')
call write_if_logic_param_changed (bmad_com%absolute_time_tracking,          bmad_com_default%absolute_time_tracking,           'bmad_com[absolute_time_tracking]')
call write_if_logic_param_changed (bmad_com%convert_to_kinetic_momentum,     bmad_com_default%convert_to_kinetic_momentum,      'bmad_com[convert_to_kinetic_momentum]')
call write_if_logic_param_changed (bmad_com%aperture_limit_on,               bmad_com_default%aperture_limit_on,                'bmad_com[aperture_limit_on]')

call write_if_logic_param_changed (ptc_com%use_orientation_patches,     ptc_com_default%use_orientation_patches,      'ptc_com[use_orientation_patches]')
call write_if_logic_param_changed (ptc_com%print_info_messages,         ptc_com_default%print_info_messages,          'ptc_com[print_info_messages]')
call write_if_real_param_changed  (ptc_com%cut_factor,                  ptc_com_default%cut_factor,                   'ptc_com[cut_factor]')
call write_if_int_param_changed   (ptc_com%max_fringe_order,            ptc_com_default%max_fringe_order,             'ptc_com[max_fringe_order]')
call write_if_real_param_changed  (ptc_com%vertical_kick,               ptc_com_default%vertical_kick,                'ptc_com[vertical_kick]')
call write_if_int_param_changed   (ptc_com%old_integrator,              ptc_com_default%old_integrator,               'ptc_com[old_integrator]')
call write_if_logic_param_changed (ptc_com%exact_model,                 ptc_com_default%exact_model,                  'ptc_com[exact_model]')
call write_if_logic_param_changed (ptc_com%exact_misalign,              ptc_com_default%exact_misalign,               'ptc_com[exact_misalign]')
call write_if_logic_param_changed (ptc_com%translate_patch_drift_time,  ptc_com_default%translate_patch_drift_time,   'ptc_com[translate_patch_drift_time]')

call write_if_real_param_changed (space_charge_com%ds_track_step,              space_charge_com_default%ds_track_step,              'space_charge_com[ds_track_step]')
call write_if_real_param_changed (space_charge_com%dt_track_step,              space_charge_com_default%dt_track_step,              'space_charge_com[dt_track_step]')
call write_if_real_param_changed (space_charge_com%cathode_strength_cutoff,    space_charge_com_default%cathode_strength_cutoff,    'space_charge_com[cathode_strength_cutoff]')
call write_if_real_param_changed (space_charge_com%rel_tol_tracking,           space_charge_com_default%rel_tol_tracking,           'space_charge_com[rel_tol_tracking]')
call write_if_real_param_changed (space_charge_com%abs_tol_tracking,           space_charge_com_default%abs_tol_tracking,           'space_charge_com[abs_tol_tracking]')
call write_if_real_param_changed (space_charge_com%beam_chamber_height,        space_charge_com_default%beam_chamber_height,        'space_charge_com[beam_chamber_height]')
call write_if_real_param_changed (space_charge_com%lsc_sigma_cutoff,           space_charge_com_default%lsc_sigma_cutoff,           'space_charge_com[lsc_sigma_cutoff]')
call write_if_real_param_changed (space_charge_com%particle_sigma_cutoff,      space_charge_com_default%particle_sigma_cutoff,      'space_charge_com[particle_sigma_cutoff]')
call write_if_int_param_changed  (space_charge_com%n_bin,                      space_charge_com_default%n_bin,                      'space_charge_com[n_bin]')
call write_if_int_param_changed  (space_charge_com%particle_bin_span,          space_charge_com_default%particle_bin_span,          'space_charge_com[particle_bin_span]')
call write_if_int_param_changed  (space_charge_com%n_shield_images,            space_charge_com_default%n_shield_images,            'space_charge_com[n_shield_images]')
call write_if_int_param_changed  (space_charge_com%sc_min_in_bin,              space_charge_com_default%sc_min_in_bin,              'space_charge_com[sc_min_in_bin]')
call write_if_logic_param_changed (space_charge_com%lsc_kick_transverse_dependence, space_charge_com_default%lsc_kick_transverse_dependence, 'space_charge_com[lsc_kick_transverse_dependence]')
call write_if_logic_param_changed (space_charge_com%debug,                          space_charge_com_default%debug,                          'space_charge_com[debug]')

if (space_charge_com%diagnostic_output_file /= '') write (iu, '(2a)') 'space_charge_com[diagnostic_output_file] = ', quote(space_charge_com%diagnostic_output_file)

ele => lat%ele(0) 

if (ele%floor%r(1) /= 0)   write (iu, '(2a)') 'beginning[x_position]     = ', re_str(ele%floor%r(1))
if (ele%floor%r(2) /= 0)   write (iu, '(2a)') 'beginning[y_position]     = ', re_str(ele%floor%r(2))
if (ele%floor%r(3) /= 0)   write (iu, '(2a)') 'beginning[z_position]     = ', re_str(ele%floor%r(3))
if (ele%floor%theta /= 0)  write (iu, '(2a)') 'beginning[theta_position] = ', re_str(ele%floor%theta)
if (ele%floor%phi /= 0)    write (iu, '(2a)') 'beginning[phi_position]   = ', re_str(ele%floor%phi)
if (ele%floor%psi /= 0)    write (iu, '(2a)') 'beginning[psi_position]   = ', re_str(ele%floor%psi)

if (ele%s /= 0)            write (iu, '(2a)') 'beginning[s]        = ', re_str(ele%s)
if (ele%ref_time /= 0)     write (iu, '(2a)') 'beginning[ref_time] = ', re_str(ele%ref_time)

! Write beginning Twiss even for closed lattices as that is useful info.

write (iu, '(2a)')
if (ele%a%beta /= 0)        write (iu, '(2a)') 'beginning[beta_a]       = ', re_str(ele%a%beta)
if (ele%a%alpha /= 0)       write (iu, '(2a)') 'beginning[alpha_a]      = ', re_str(ele%a%alpha)
if (ele%a%phi /= 0)         write (iu, '(2a)') 'beginning[phi_a]        = ', re_str(ele%a%phi)
if (ele%x%eta /= 0)         write (iu, '(2a)') 'beginning[eta_x]        = ', re_str(ele%x%eta)
if (ele%x%etap /= 0)        write (iu, '(2a)') 'beginning[etap_x]       = ', re_str(ele%x%etap)
if (ele%b%beta /= 0)        write (iu, '(2a)') 'beginning[beta_b]       = ', re_str(ele%b%beta)
if (ele%b%alpha /= 0)       write (iu, '(2a)') 'beginning[alpha_b]      = ', re_str(ele%b%alpha)
if (ele%b%phi /= 0)         write (iu, '(2a)') 'beginning[phi_b]        = ', re_str(ele%b%phi)
if (ele%y%eta /= 0)         write (iu, '(2a)') 'beginning[eta_y]        = ', re_str(ele%y%eta)
if (ele%y%etap /= 0)        write (iu, '(2a)') 'beginning[etap_y]       = ', re_str(ele%y%etap)
if (ele%c_mat(1,1) /= 0)    write (iu, '(2a)') 'beginning[cmat_11]      = ', re_str(ele%c_mat(1,1))
if (ele%c_mat(1,2) /= 0)    write (iu, '(2a)') 'beginning[cmat_12]      = ', re_str(ele%c_mat(1,2))
if (ele%c_mat(2,1) /= 0)    write (iu, '(2a)') 'beginning[cmat_21]      = ', re_str(ele%c_mat(2,1))
if (ele%c_mat(2,2) /= 0)    write (iu, '(2a)') 'beginning[cmat_22]      = ', re_str(ele%c_mat(2,2))
if (ele%mode_flip)          write (iu, '(a)')  'beginning[mode_flip]    = T'
if (ele%a%dbeta_dpz /= 0)   write (iu, '(2a)') 'beginning[dbeta_dpz_a]  = ', re_str(ele%a%dbeta_dpz)
if (ele%b%dbeta_dpz /= 0)   write (iu, '(2a)') 'beginning[dbeta_dpz_b]  = ', re_str(ele%b%dbeta_dpz)
if (ele%a%dalpha_dpz /= 0)  write (iu, '(2a)') 'beginning[dalpha_dpz_a] = ', re_str(ele%a%dalpha_dpz)
if (ele%b%dalpha_dpz /= 0)  write (iu, '(2a)') 'beginning[dalpha_dpz_b] = ', re_str(ele%b%dalpha_dpz)
if (ele%x%deta_dpz /= 0)    write (iu, '(2a)') 'beginning[deta_dpz_x]   = ', re_str(ele%x%deta_dpz)
if (ele%y%deta_dpz /= 0)    write (iu, '(2a)') 'beginning[deta_dpz_y]   = ', re_str(ele%y%deta_dpz)
if (ele%x%detap_dpz /= 0)   write (iu, '(2a)') 'beginning[detap_dpz_x]  = ', re_str(ele%x%detap_dpz)
if (ele%y%detap_dpz /= 0)   write (iu, '(2a)') 'beginning[detap_dpz_y]  = ', re_str(ele%y%detap_dpz)

! particle_start. Note: For an open geometry, orbit0 should be the same as lat%particle_start

if (lat%param%geometry == closed$ .and. present(orbit0)) then
  if (orbit0%vec(1) /= 0) write (iu, '(2a)') 'particle_start[x]  = ', re_str(orbit0%vec(1))
  if (orbit0%vec(2) /= 0) write (iu, '(2a)') 'particle_start[px] = ', re_str(orbit0%vec(2))
  if (orbit0%vec(3) /= 0) write (iu, '(2a)') 'particle_start[y]  = ', re_str(orbit0%vec(3))
  if (orbit0%vec(4) /= 0) write (iu, '(2a)') 'particle_start[py] = ', re_str(orbit0%vec(4))
  if (orbit0%vec(5) /= 0) write (iu, '(2a)') 'particle_start[z]  = ', re_str(orbit0%vec(5))
  if (orbit0%vec(6) /= 0) write (iu, '(2a)') 'particle_start[pz] = ', re_str(orbit0%vec(6))

  if (orbit0%spin(1) /= 0) write (iu, '(2a)') 'particle_start[spin_x] = ', re_str(orbit0%spin(1))
  if (orbit0%spin(2) /= 0) write (iu, '(2a)') 'particle_start[spin_y] = ', re_str(orbit0%spin(2))
  if (orbit0%spin(3) /= 0) write (iu, '(2a)') 'particle_start[spin_z] = ', re_str(orbit0%spin(3))

else
  if (lat%particle_start%vec(1) /= 0) write (iu, '(2a)') 'particle_start[x]  = ', re_str(lat%particle_start%vec(1))
  if (lat%particle_start%vec(2) /= 0) write (iu, '(2a)') 'particle_start[px] = ', re_str(lat%particle_start%vec(2))
  if (lat%particle_start%vec(3) /= 0) write (iu, '(2a)') 'particle_start[y]  = ', re_str(lat%particle_start%vec(3))
  if (lat%particle_start%vec(4) /= 0) write (iu, '(2a)') 'particle_start[py] = ', re_str(lat%particle_start%vec(4))
  if (lat%particle_start%vec(5) /= 0) write (iu, '(2a)') 'particle_start[z]  = ', re_str(lat%particle_start%vec(5))
  if (lat%particle_start%vec(6) /= 0) write (iu, '(2a)') 'particle_start[pz] = ', re_str(lat%particle_start%vec(6))
  if (lat%particle_start%t /= 0)      write (iu, '(2a)') 'particle_start[t]  = ', re_str(lat%particle_start%t)

  if (lat%particle_start%spin(1) /= 0) write (iu, '(2a)') 'particle_start[spin_x] = ', re_str(lat%particle_start%spin(1))
  if (lat%particle_start%spin(2) /= 0) write (iu, '(2a)') 'particle_start[spin_y] = ', re_str(lat%particle_start%spin(2))
  if (lat%particle_start%spin(3) /= 0) write (iu, '(2a)') 'particle_start[spin_z] = ', re_str(lat%particle_start%spin(3))
endif

! Named constants

write (iu, '(a)')

do i = 1, lat%n_control_max
  if (.not. allocated(lat%control(i)%stack)) cycle
  stack => lat%control(i)%stack
  do j = 1, size(stack)
    if (stack(j)%type == end_stack$) exit
    if (stack(j)%type /= variable$) cycle
    if (stack(j)%name == '') cycle
    if (any(stack(j)%name == physical_const_list%name)) cycle
    call find_index(stack(j)%name, str_index, ix, add_to_list = .true., has_been_added = has_been_added)
    if (.not. (has_been_added)) cycle  ! Avoid duplicates
    write (iu, '(3a)') trim(stack(j)%name), ' = ', re_str(stack(j)%value)
  enddo
enddo

! Element stuff

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)')

n_names = 0
n = lat%n_ele_max
allocate (names(n), an_indexx(n), named_eles(n))

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (ib > 0) then
    write (iu, '(a)')
    write (iu, '(a)')  '!-------------------------------------------------------'
    write (iu, '(2a)') '! Branch: ', trim(branch%name)
    write (iu, '(a)')
  endif

  ele_loop: do ie = 1, branch%n_ele_max

    ele => branch%ele(ie)
    if (ie == ele%branch%n_ele_track .and. ele%name == 'END' .and. ele%key == marker$) cycle
    if (ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$) cycle   ! Handled in next section.
    if (ele%key == null_ele$) cycle

    ele_dflt => ele_default(ele%key) ! Element with default attributes.

    ! Superposition stragegy: Swap drifts for super_slaves

    multi_lord => pointer_to_multipass_lord (ele, ix_pass) 
    if (ele%lord_status == super_lord$ .and. ix_pass > 0) cycle
    if (ele%slave_status == super_slave$ .and. ix_pass > 1) cycle

    if (ele%slave_status == super_slave$) then
      lord => pointer_to_lord(ele, 1)
      slave => pointer_to_slave(lord, 1)
      slave2 => pointer_to_slave(lord, lord%n_slave)
      write (iu, '(2(a, i0), 2a)') 'slave_drift_', ib, '_', ele%ix_ele, ': drift, l = ', re_str(ele%value(l$))
      cycle
    endif

    if (ix_pass > 0) cycle

    ! Do not write anything for elements that have a duplicate name.

    call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles)
    if (.not. has_been_added) cycle

    ! Girder

    if (ele%key == girder$) then
      write (line, '(2a)') trim(ele%name), ': girder = {'
      do j = 1, ele%n_slave
        slave => pointer_to_slave(ele, j)
        if (j == ele%n_slave) then
          write (line, '(3a)') trim(line), trim(slave%name), '}'
        else
          write (line, '(3a)') trim(line), trim(slave%name), ', '
        endif
      enddo

    ! Everything but a girder

    else
      line = trim(ele%name) // ': ' // key_name(ele%key)
    endif

    ! Branch

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      n = nint(ele%value(ix_to_branch$))
      line = trim(line) // ', to_line = ' // trim(lat%branch(n)%name)
      if (ele%value(ix_to_element$) > 0) then
        i = nint(ele%value(ix_to_element$))
        line = trim(line) // ', to_element = ' // trim(lat%branch(n)%ele(i)%name)
      endif
    endif

    ! Other elements

    if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
    if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
    if (associated(ele%descrip)) line = trim(line) // ', descrip = "' // trim(ele%descrip) // '"'

    ! AC_Kicker

    if (associated(ele%ac_kick)) then
      ac => ele%ac_kick
      if (allocated(ac%amp_vs_time)) then
        line = trim(line) // ', amp_vs_time = {(' // re_str(ac%amp_vs_time(1)%time) // &
                             ', ' // re_str(ac%amp_vs_time(1)%amp) // ')'
        do i = 2, size(ac%amp_vs_time)
          line = trim(line) // ', (' // re_str(ac%amp_vs_time(i)%time) // &
                             ', ' // re_str(ac%amp_vs_time(i)%amp) // ')'
        enddo
        line = trim(line) // '}'

      else
        line = trim(line) // ', frequencies = {(' // re_str(ac%frequency(1)%f) // &
                  ', ' // re_str(ac%frequency(1)%amp) // ', ' // re_str(ac%frequency(1)%phi)  // ')'
        do i = 2, size(ac%frequency)
          line = trim(line) // ', (' // re_str(ac%frequency(i)%f) // &
                  ', ' // re_str(ac%frequency(i)%amp) // ', ' // re_str(ac%frequency(i)%phi) // ')'
        enddo
        line = trim(line) // '}'
      endif
    endif

    ! Wall3d

    if (associated(ele%wall3d)) then

      ! First find out out if a wall file has been written
      found = .false.
      do ibr = 0, ubound(lat%branch, 1)
        branch2 => lat%branch(ibr)
        imax = branch2%n_ele_max
        if (ibr == branch%ix_branch) imax = ie-1
        do ie2 = 1, imax
          ele2 => branch2%ele(ie2)
          if (ele2%slave_status == multipass_slave$) cycle
          if (.not. associated(ele2%wall3d)) cycle
          if (.not. associated(ele2%wall3d, ele%wall3d)) cycle
          found = .true.
          exit
        enddo
        if (found) exit
      enddo

      if (found) then
        call str_downcase(name, ele2%name)
        line = trim(line) // ', wall = call::wall_' // trim(name)
      else
        call str_downcase(name, ele%name)
        line = trim(line) // ', wall = call::wall_' // trim(name)
        iu2 = lunget()
 
        open (iu2, file = trim(path) // '/wall_' // trim(name))
        write (iu2, '(a)') '{'

        if (ele%key /= diffraction_plate$ .and. ele%key /= mask$) then
          write (iu2, '(2x, 3a)') 'ele_anchor_pt = ', trim(anchor_pt_name(ele%wall3d(1)%ele_anchor_pt)), ','
        endif

        do i = 1, size(ele%wall3d(1)%section)
          section => ele%wall3d(1)%section(i)
          write (iu2, '(2x, a)')   'section = {'

          if (ele%key == diffraction_plate$ .or. ele%key == mask$) then
            write (iu2, '(4x, 3a)') 'type = ', trim(wall3d_section_type_name(section%type)), ','
          else
            write (iu2, '(4x, 3a)')  's     = ', re_str(section%s), ','
            if (section%dr_ds /= real_garbage$) write (iu2, '(4x, 3a)')  'dr_ds = ', re_str(section%s), ','
          endif

          if (any(section%r0 /= 0)) then
            write (iu2, '(4x, 3a)')  'r0    = ', trim(array_re_str(section%r0, '()')), ','
          endif

          if (section%vertices_state == shifted_to_relative$) then
            write (iu2, '(4x, 3a)')  'absolute_vertices = True,'
          endif

          end_str = ','
          do j = 1, size(section%v)
            if (j == size(section%v)) then
              end_str = '},'
              if (i == size(ele%wall3d(1)%section)) end_str = '}}'
            endif
            v => section%v(j)

            if (section%vertices_state == shifted_to_relative$) then
              x = v%x + section%r0(1)
              y = v%y + section%r0(2)
            else
              x = v%x
              y = v%y
            endif

            if (v%tilt /= 0) then
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([x, y, v%radius_x, v%radius_y, v%tilt], '{}')), end_str
            elseif (v%radius_y /= 0) then
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([x, y, v%radius_x, v%radius_y], '{}')), end_str
            elseif (v%radius_x /= 0) then
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([x, y, v%radius_x], '{}')), end_str
            else
              write (iu2, '(4x, a, i0, 3a)') 'v(', j, ') = ', &
                    trim(array_re_str([x, y], '{}')), end_str
            endif
          enddo
          if (len_trim(line) > 1000) call write_lat_line(line, iu, .false.)
        enddo
        close (iu2)

      endif
    endif

    ! field overlap

    if (ele%n_slave_field /= 0) then
      slave => pointer_to_slave (ele, 1, lord_type = field_lord$)
      line = trim(line) // ', field_overlaps = {' // slave%name
      do n = 2, ele%n_slave_field
        slave => pointer_to_slave (ele, n, lord_type = field_lord$)
        line = trim(line) // ', ' // slave%name
      enddo
      line = trim(line) // '}'
    endif

    ! Foil

    if (associated(ele%foil)) then
      if (size(ele%foil%material) > 1) then
        if (any(ele%foil%material%density /= real_garbage$)) then
          line = trim(line) // ', density = ('
          do n = 1, size(ele%foil%material)
            if (n == 1) then; line = trim(line) // re_str(ele%foil%material(n)%density)
            else;             line = trim(line) // ', ' // re_str(ele%foil%material(n)%density)
            endif
          enddo
          line = trim(line) // ')'
        endif

        if (any(ele%foil%material%area_density /= real_garbage$)) then
          line = trim(line) // ', area_density = ('
          do n = 1, size(ele%foil%material)
            if (n == 1) then; line = trim(line) // re_str(ele%foil%material(n)%area_density)
            else;             line = trim(line) // ', ' // re_str(ele%foil%material(n)%area_density)
            endif
          enddo
          line = trim(line) // ')'
        endif

        if (any(ele%foil%material%radiation_length /= real_garbage$)) then
          line = trim(line) // ', radiation_length = ('
          do n = 1, size(ele%foil%material)
            if (n == 1) then; line = trim(line) // re_str(ele%foil%material(n)%radiation_length)
            else;             line = trim(line) // ', ' // re_str(ele%foil%material(n)%radiation_length)
            endif
          enddo
          line = trim(line) // ')'
        endif

      else
        material => ele%foil%material(1)
        if (material%density /= real_garbage$)          line = trim(line) // ', density = ' // re_str(material%density)
        if (material%area_density /= real_garbage$)     line = trim(line) // ', area_density = ' // re_str(material%area_density)
        if (material%radiation_length /= real_garbage$) line = trim(line) // ', radiation_length = ' // re_str(material%radiation_length)
      endif
    endif

    ! Cartesian_map.

    if (associated(ele%cartesian_map)) then
      do im = 1, size(ele%cartesian_map)
        ct_map => ele%cartesian_map(im)

        call find_matching_fieldmap (ct_map%ptr%file, ele, cartesian_map$, ele2, ix_ptr, ignore_slaves = .true.) 

        if (integer_option(binary$, output_form) == one_file$) then
          line = trim(line) // ', cartesian_map ='
          call write_lat_line (line, iu, .true.)
          call write_this_cartesian_map (ct_map, ele, iu, line)

        elseif (ix_ptr > 0) then  ! A file has been created so refer to that

          call form_this_fieldmap_name(fname, '.cartesian_map', ele2, ix_ptr, ascii$)
          write (line, '(3a)')  trim(line), ', cartesian_map = call::', trim(fname)

        else
          call form_this_fieldmap_name(fname, '.cartesian_map', ele, im, ascii$)
          line = trim(line) // ', cartesian_map = call::' // trim(fname)
          fname = trim(path) // '/' // trim(fname)
          iu2 = lunget()
          open (iu2, file = fname)
          call write_this_cartesian_map (ct_map, ele, iu2)
          close (iu2)
        endif
      enddo
    endif

    ! cylindrical_map

    if (associated(ele%cylindrical_map)) then
      do im = 1, size(ele%cylindrical_map)
        cl_map => ele%cylindrical_map(im)

        call find_matching_fieldmap (cl_map%ptr%file, ele, cylindrical_map$, ele2, ix_ptr, ignore_slaves = .true.) 

        if (integer_option(binary$, output_form) == one_file$) then
          line = trim(line) // ', cylindrical_map ='
          call write_lat_line (line, iu, .true.)
          call write_this_cylindrical_map (cl_map, ele, iu, line)

        elseif (ix_ptr > 0) then
          call form_this_fieldmap_name(fname, '.cylindrical_map', ele2, ix_ptr, ascii$)
          write (line, '(3a)')  trim(line), ', cylindrical_map = call::', trim(fname)

        else
          call form_this_fieldmap_name(fname, '.cylindrical_map', ele, im, ascii$)
          line = trim(line) // ', cylindrical_map = call::' // trim(fname)
          fname = trim(path) // '/' // trim(fname)
          iu2 = lunget()
          open (iu2, file = fname)
          call write_this_cylindrical_map (cl_map, ele, iu2)
          close (iu2)
        endif
      enddo
    endif

    ! gen_grad_map

    if (associated(ele%gen_grad_map)) then
      do im = 1, size(ele%gen_grad_map)
        gg_map => ele%gen_grad_map(im)

        ! First find out out if an file has been written
        call find_matching_fieldmap (gg_map%file, ele, gen_grad_map$, ele2, ix_ptr, ignore_slaves = .true.) 

        if (integer_option(binary$, output_form) == one_file$) then
          line = trim(line) // ', gen_grad_map_map ='
          call write_lat_line (line, iu, .true.)
          call write_this_gen_grad_map_map (gg_map, ele, iu, line)

        elseif (ix_ptr > 0) then
          call form_this_fieldmap_name(fname, '.gen_grad_map', ele2, ix_ptr, ascii$)
          write (line, '(3a)')  trim(line), ', gen_grad_map = call::', trim(fname)

        else
          call form_this_fieldmap_name(fname, '.gen_grad_map', ele, im, ascii$)
          line = trim(line) // ', gen_grad_map = call::' // trim(fname)
          fname = trim(path) // '/' // trim(fname)
          iu2 = lunget()
          open (iu2, file = fname, recl = 500)
          call write_this_gen_grad_map_map (gg_map, ele, iu2)
          close (iu2)
        endif
      enddo
    endif

    ! grid_field

    if (associated(ele%grid_field)) then
      do im = 1, size(ele%grid_field)
        g_field => ele%grid_field(im)

        ! First find out out if an file has been written
        call find_matching_fieldmap (g_field%ptr%file, ele, grid_field$, ele2, ix_ptr, ignore_slaves = .true.) 

        if (integer_option(binary$, output_form) == one_file$) then
          line = trim(line) // ', grid_fieldmap ='
          call write_lat_line (line, iu, .true.)
          call write_this_grid_fieldmap (g_field, ele, iu, line)

        elseif (ix_ptr > 0) then
          call form_this_fieldmap_name(fname, '.grid_field', ele2, ix_ptr, output_form)
          write (line, '(3a)')  trim(line), ', grid_field = call::', trim(fname)

        else
          call form_this_fieldmap_name(fname, '.grid_field', ele, im, output_form)
          line = trim(line) // ', grid_field = call::' // trim(fname)
          fname = trim(path) // '/' // trim(fname)

          if (integer_option(binary$, output_form) == binary$) then
            call hdf5_write_grid_field (fname, ele, ele%grid_field(im:im), err_flag)
          else
            iu2 = lunget()
            open (iu2, file = fname)
            call write_this_grid_fieldmap (g_field, ele, iu2)
            close (iu2)
          endif
        endif
      enddo
    endif

    ! Wake

    if (associated(ele%wake)) then
      lr => ele%wake%lr
      if (size(lr%mode) /= 0) then
        line = trim(line) // ', lr_wake = @'
        if (lr%freq_spread /= 0) line = trim(line) // ', freq_spread = ' // re_str(lr%freq_spread)
        if (.not. lr%self_wake_on) line = trim(line) // ', self_wake_on = ' // logic_str(lr%self_wake_on)
        if (lr%amp_scale /= 1) line = trim(line) // ', amp_scale = ' // re_str(lr%amp_scale)
        if (lr%time_scale /= 1) line = trim(line) // ', time_scale = ' // re_str(lr%time_scale)
        if (lr%t_ref /= 0) line = trim(line) // ', t_ref = ' // re_str(lr%t_ref)
        do i = 1, size(lr%mode)
          lrm => lr%mode(i)
          line = trim(line) // ', mode = {' // re_str(lrm%freq_in) // ', ' // re_str(lrm%R_over_Q) // &
                     ', ' // re_str(lrm%damp) // ', ' // re_str(lrm%phi) // ', ' // int_str(lrm%m) 

          if (lrm%polarized) then
            line = trim(line) // ', ' // re_str(lrm%angle)
          else
            line = trim(line) // ', unpolarized'
          endif

          if (lrm%b_sin == 0 .and. lrm%b_cos == 0 .and. lrm%a_sin == 0 .and. lrm%a_cos == 0) then
            line = trim(line) // '}'
          else
            line = trim(line) // ', ' // re_str(lrm%b_sin) // ', ' // re_str(lrm%b_cos) // ', ' // &
                                         re_str(lrm%a_sin) // ', ' // re_str(lrm%a_cos) // '}'
          endif

          if (i == 1) then
            ix = index(line, '@,')
            line = line(1:ix-1) // '{' //line(ix+2:)
          endif

          if (len_trim(line) > 1000) call write_lat_line(line, iu, .false.)
        enddo
        line = trim(line) // '}'
      endif

      sr => ele%wake%sr
      if (size(sr%long) /= 0 .or. size(sr%trans) /= 0 .or. sr%z_long%dz /= 0) then
        line = trim(line) // ', sr_wake = @'
        if (sr%z_max /= 0) line = trim(line) // ', z_max = ' // re_str(sr%z_max)
        if (sr%amp_scale /= 1) line = trim(line) // ', amp_scale = ' // re_str(sr%amp_scale)
        if (sr%z_scale /= 1) line = trim(line) // ', z_scale = ' // re_str(sr%z_scale)

        do i = 1, size(sr%long)
          ix = index(line, '@,')
          if (ix /= 0) line = line(1:ix-1) // '{' //line(ix+2:)

          srm => sr%long(i)
          line = trim(line) // ', longitudinal = {' // re_str(srm%amp) // ', ' // re_str(srm%damp) // ', ' // re_str(srm%k) // &
                               ', ' // re_str(srm%phi) // ', ' // trim(sr_longitudinal_position_dep_name(srm%position_dependence)) // '}'
          if (len_trim(line) > 1000) call write_lat_line(line, iu, .false.)
        enddo 

        do i = 1, size(sr%trans)
          ix = index(line, '@,')
          if (ix /= 0) line = line(1:ix-1) // '{' //line(ix+2:)
          srm => sr%trans(i)
          line = trim(line) // ', transverse = {' // re_str(srm%amp) // ', ' // re_str(srm%damp) // ', ' // re_str(srm%k) // &
                               ', ' // re_str(srm%phi) // ', ' // trim(sr_transverse_polarization_name(srm%polarization)) // &
                               ', ' // trim(sr_transverse_position_dep_name(srm%position_dependence)) // '}'
          if (len_trim(line) > 1000) call write_lat_line(line, iu, .false.)
        enddo

        srz => sr%z_long
        if (srz%dz /= 0) then
          ix = index(line, '@,')
          if (ix /= 0) line = line(1:ix-1) // '{' //line(ix+2:)
          name = trim(ele%name) // '.sr_z_long'
          f = 1
          if (srz%time_based) f = 1.0_rp / c_light
          line = trim(line) // ', z_long = {time_based = ' // logic_str(srz%time_based) // ', position_dependence = ' // &
                  trim(sr_transverse_position_dep_name(srz%position_dependence)) // ', smoothing_sigma = ' // re_str(f*srz%smoothing_sigma) // &
                  ', w = {call::' // trim(name) // '}'

          iu2 = lunget()
          open (iu2, file = trim(path) // '/' // trim(name))
          n = size(srz%w)
          do i = 1, n
            z = -srz%z0 + (i-1) * srz%dz
            if (srz%time_based) then
              write (iu2, '(es16.8, es20.12, a)') f*z, srz%w(n+1-i), ','
            else
              write (iu2, '(es16.8, es20.12, a)') z, srz%w(i), ','
            endif
          enddo
          close(iu2)
          line = trim(line) // '}'
        endif

        line = trim(line) // '}'

        ix = index(line, '@,')
        if (ix /= 0) line = line(1:ix-1) // '{' //line(ix+2:)
      endif

      call write_lat_line (line, iu, .false.)
    endif

    ! Decide if x1_limit, etc. are to be output directly or combined. 

    x_lim = ele%value(x1_limit$) 
    x_lim_good = .false.
    if (x_lim /=0 .and. ele%value(x2_limit$) == x_lim) x_lim_good = .true.

    y_lim = ele%value(y1_limit$) 
    y_lim_good = .false.
    if (y_lim /=0 .and. ele%value(y2_limit$) == y_lim) y_lim_good = .true.

    !----------------------------------------------------------------------------
    ! Write the element attributes.

    fid = nint(ele%value(fiducial_pt$))
    attribute_loop: do j = 1, num_ele_attrib$
      attrib = attribute_info(ele, j)
      val = ele%value(j)
      if (val == ele_dflt%value(j)) cycle
      if (ele%key == sbend$) then
        if (j == fintx$ .and. ele%value(fintx$) == ele%value(fint$)) cycle
        if (j == hgapx$ .and. ele%value(hgapx$) == ele%value(hgap$)) cycle

        if (j == l$ .and. (fid == entrance_end$ .or. fid == entrance_end$)) cycle
        if (j == l_rectangle$ .and. (fid == none_pt$ .or. fid == center_pt$)) cycle
      endif
      if (j == check_sum$) cycle
      if (x_lim_good .and. (j == x1_limit$ .or. j == x2_limit$)) cycle
      if (y_lim_good .and. (j == y1_limit$ .or. j == y2_limit$)) cycle
      if (.not. attribute_free (ele, attrib%name, .false., .true.)) cycle
      if ((attrib%name == 'P0C' .or. attrib%name == 'P0C_START') .and. &
                          (ele%lord_status /= multipass_lord$ .or. nint(ele%value(multipass_ref_energy$)) == first_pass$)) cycle

      if (ele%key == mask$ .or. ele%key == diffraction_plate$) then
        if (j == x1_limit$ .or. j == x2_limit$ .or. j == y1_limit$ .or. j == y2_limit$)  cycle
      endif

      ! Default for ds_step and integrator_order is determined by attribute_bookkeeper based upon the
      ! settings of other parameters like the element's strength.
      if (attrib%name == 'DS_STEP' .or. attrib%name == 'INTEGRATOR_ORDER') then
        call transfer_ele (ele, this_ele, .true.) 
        this_ele%value(ds_step$) = 0
        this_ele%value(num_steps$) = 0
        this_ele%value(integrator_order$) = 0
        call attribute_bookkeeper (this_ele, .true.)
        if (attrib%name == 'DS_STEP' .and. abs(val - this_ele%value(ds_step$)) < 1e-6*val) cycle
        if (attrib%name == 'INTEGRATOR_ORDER' .and. val == this_ele%value(integrator_order$)) cycle        
      endif

      if (attrib%name == 'E_TOT') cycle        ! Will use p0c instead.
      if (attrib%name == 'E_TOT_START') cycle  ! Will use p0c_start instead.
      if (attrib%name == 'E_TOT_STRONG') cycle  ! Will use pc_strong instead.
      if (attrib%name == null_name$) then
        call out_io (s_error$, r_name, 'ELEMENT: ' // ele%name, 'HAS AN UNKNOWN ATTRIBUTE INDEX: \i0\ ', i_array = [j])
        if (global_com%exit_on_error) call err_exit
        return
      endif

      if (attrib%name == 'COUPLER_AT') then
        if (nint(val) /= downstream_end$) then
          line = trim(line) // ', coupler_at = ' // end_at_name(nint(val))
        endif
        cycle
      endif

      select case (attribute_type(attrib%name))
      case (is_logical$)
        write (line, '(4a, l1)') trim(line), ', ', trim(attrib%name), ' = ', (val /= 0)
      case (is_integer$)
        write (line, '(4a, i0)') trim(line), ', ', trim(attrib%name), ' = ', int(val)
      case (is_real$)
        line = trim(line) // ', ' // trim(attrib%name) // ' = ' // re_str(val)
      case (is_switch$)
        name = switch_attrib_value_name (attrib%name, val, ele, is_default)
          if (.not. is_default) then
            line = trim(line) // ', ' // trim(attrib%name) // ' = ' // name
          endif
      end select

    enddo attribute_loop ! attribute loop

    ! Custom attributes

    do j = 1, custom_attribute_num$
      attrib = attribute_info(ele, j+custom_attribute0$)
      if (attrib%name(1:1) == '!') cycle
      val = value_of_attribute(ele, attrib%name, err_flag)
      if (val == 0) cycle
      line = trim(line) // ', ' // trim(attrib%name) // ' = ' // re_str(val)
    enddo

    !----------------------------------------------------------------------------
    ! Print the combined limits if needed.

    if (x_lim_good .and. y_lim_good .and. x_lim == y_lim) then
      line = trim(line) // ', aperture = ' // re_str(x_lim)
    else
      if (x_lim_good) line = trim(line) // ', x_limit = ' // re_str(x_lim)
      if (y_lim_good) line = trim(line) // ', y_limit = ' // re_str(y_lim)
    endif

    ! Encode methods, etc.

    if (has_attribute (ele, 'MAT6_CALC_METHOD') .and. (ele%mat6_calc_method /= ele_dflt%mat6_calc_method)) &
                                      line = trim(line) // ', mat6_calc_method = ' // mat6_calc_method_name(ele%mat6_calc_method)
    if (has_attribute (ele, 'TRACKING_METHOD') .and. (ele%tracking_method /= ele_dflt%tracking_method)) &
                                      line = trim(line) // ', tracking_method = ' // tracking_method_name(ele%tracking_method)
    if (has_attribute (ele, 'SPIN_TRACKING_METHOD') .and. (ele%spin_tracking_method /= ele_dflt%spin_tracking_method)) &
                                      line = trim(line) // ', spin_tracking_method = ' // spin_tracking_method_name(ele%spin_tracking_method)
    if (has_attribute (ele, 'CSR_METHOD') .and. (ele%csr_method /= ele_dflt%csr_method)) &
                                      line = trim(line) // ', csr_method = ' // csr_method_name(ele%csr_method)
    if (has_attribute (ele, 'SPACE_CHARGE_METHOD') .and. (ele%space_charge_method /= ele_dflt%space_charge_method)) &
                                      line = trim(line) // ', space_charge_method = ' // space_charge_method_name(ele%space_charge_method)
    if (has_attribute (ele, 'PTC_INTEGRATION_TYPE') .and. (ele%ptc_integration_type /= ele_dflt%ptc_integration_type)) &
                                      line = trim(line) // ', ptc_integration_type = ' // ptc_integration_type_name(ele%ptc_integration_type)
    if (has_attribute (ele, 'FIELD_CALC') .and. (ele%field_calc /= ele_dflt%field_calc)) &
                                      line = trim(line) // ', field_calc = ' // field_calc_name(ele%field_calc)

    if (has_attribute (ele, 'APERTURE_AT') .and. (ele%aperture_at /= ele_dflt%aperture_at)) &
                                      line = trim(line) // ', aperture_at = ' // aperture_at_name(ele%aperture_at)
    if (has_attribute (ele, 'APERTURE_TYPE') .and. (ele%aperture_type /= ele_dflt%aperture_type)) &
                                      line = trim(line) // ', aperture_type = ' // aperture_type_name(ele%aperture_type)

    if (has_attribute (ele, 'SYMPLECTIFY') .and. ele%symplectify) line = trim(line) // ', symplectify = T'

    if (has_attribute (ele, 'FIELD_MASTER') .and. (ele%field_master .neqv. ele_dflt%field_master)) &
                                      write (line, '(2a, l1)') trim(line), ', field_master = ', ele%field_master
    if (has_attribute (ele, 'IS_ON') .and. (ele%is_on .neqv. ele_dflt%is_on)) &
                                      write (line, '(2a, l1)') trim(line), ', is_on = ', ele%is_on
    if (has_attribute (ele, 'SCALE_MULTIPOLES') .and. (ele%scale_multipoles .neqv. ele_dflt%scale_multipoles)) &
                                      write (line, '(2a, l1)') trim(line), ', scale_multipoles = ', ele%scale_multipoles
    if (has_attribute (ele, 'MULTIPOLES_ON') .and. (ele%multipoles_on .neqv. ele_dflt%multipoles_on)) &
                                      write (line, '(2a, l1)') trim(line), ', multipoles_on = ', ele%multipoles_on
    if (has_attribute (ele, 'TAYLOR_MAP_INCLUDES_OFFSETS') .and. (ele%taylor_map_includes_offsets .neqv. ele_dflt%taylor_map_includes_offsets)) &
                                      write (line, '(2a, l1)') trim(line), ', taylor_map_includes_offsets = ', ele%taylor_map_includes_offsets
    if (has_attribute (ele, 'OFFSET_MOVES_APERTURE') .and. (ele%offset_moves_aperture .neqv. ele_dflt%offset_moves_aperture)) &
                                      write (line, '(2a, l1)') trim(line), ', offset_moves_aperture = ', ele%offset_moves_aperture

    if (has_attribute (ele, 'ORIGIN_ELE') .and. ele%component_name /= '')      line = trim(line) // ', origin_ele = ' // ele%component_name 
    if (has_attribute (ele, 'CRYSTAL_TYPE') .and. ele%component_name /= '')    line = trim(line) // ', crystal_type = ' // ele%component_name 
    if (has_attribute (ele, 'MATERIAL_TYPE') .and. ele%component_name /= '')   line = trim(line) // ', material_type = ' // ele%component_name 
    if (has_attribute (ele, 'PHYSICAL_SOURCE') .and. ele%component_name /= '') line = trim(line) // ', physical_source = ' // ele%component_name 


    call write_lat_line (line, iu, .false.)

    ! Encode taylor map. Hybrid elements do not have default terms.

    if (ele%key == taylor$ .or. (ele%key == hybrid$ .and. associated(ele%taylor(1)%term))) then
      do j = 1, 6
        unit_found = .false.
        unit = 0
        unit(j:j) = 1

        do k = 1, size(ele%taylor(j)%term)
          tm = ele%taylor(j)%term(k)
          write_term = .false.
          if (all(tm%expn == unit)) then
            unit_found = .true.
            if (tm%coef /= 1) write_term = .true.
          else
            write_term = .true.
          endif

          if (write_term .or. ele%key == hybrid$) then
            if (sum(tm%expn) < 6) then
              name = ''
              do ix = 1, 6
                do n = 1, tm%expn(ix)
                  write (name, '(a, i1)') trim(name), ix
                enddo
              enddo
              write (line, '(2a, i0, 5a)') trim(line), ', {', j, ': ', re_str(tm%coef), ' | ', trim(name), '}'
            else
              write (line, '(2a, i0, 3a, 6(1x, i0), a)') trim(line), ', {', j, ': ', re_str(tm%coef), ',', tm%expn, '}'
            endif
          endif

          call write_lat_line (line, iu, .false.)
        enddo

        if (ele%key == taylor$ .and. .not. unit_found) write (line, '(2a, i0, a, 6i2, a)') trim(line), ', {', j, ': 0,', tm%expn, '}'
      enddo

      do j1 = 0, 3
        if (.not. associated(ele%spin_taylor(j1)%term)) cycle
        do k = 1, size(ele%spin_taylor(j1)%term)
          tm = ele%spin_taylor(j1)%term(k)
          write (line, '(6a, 6i2, a)') trim(line), ', {', spin_quat_name(j1), ': ', re_str(tm%coef), ',', tm%expn, '}'
        enddo
      enddo

      if (any(ele%taylor%ref /= 0)) then
        write (line, '(16a)') trim(line), ', ref_orbit = (', &
                re_str(ele%taylor(1)%ref), ', ', re_str(ele%taylor(2)%ref), ', ', &
                re_str(ele%taylor(3)%ref), ', ', re_str(ele%taylor(4)%ref), ', ', &
                re_str(ele%taylor(5)%ref), ', ', re_str(ele%taylor(6)%ref), ')'
      endif
    endif

    ! Encode multipoles

    if (associated(ele%a_pole)) then
      do j = 0, ubound(ele%a_pole, 1)
        if (ele%a_pole(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+a0$)) // ' = ' // re_str(ele%a_pole(j))
        if (ele%b_pole(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+b0$)) // ' = ' // re_str(ele%b_pole(j))
      enddo
    endif
    
    if (associated(ele%a_pole_elec)) then
      do j = 0, ubound(ele%a_pole_elec, 1)
        if (ele%a_pole_elec(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+a0_elec$)) // ' = ' // re_str(ele%a_pole_elec(j))
        if (ele%b_pole_elec(j) /= 0) line = trim(line) // ', ' // &
                trim(attribute_name(ele, j+b0_elec$)) // ' = ' // re_str(ele%b_pole_elec(j))
      enddo
    endif
    
    call write_lat_line (line, iu, .true.)  

  enddo ele_loop
enddo  ! branch loop

!----------------------------------------------------------
! Overlays, groups, and superimpose

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)') '! Overlays, groups, rampers, and superimpose'
write (iu, '(a)')

do ie = lat%n_ele_track+1, lat%n_ele_max
  ele => lat%ele(ie)

  ! Superimpose. Only first pass multipass_slaves are superimpsed and with first pass multipass_slaves
  ! use the lord name in the superimpose statement.

  if (ele%lord_status == super_lord$) then
    if (ele%slave_status == multipass_slave$) then
      ele2 => pointer_to_multipass_lord (ele, ix_pass) 
      if (ix_pass /= 1) cycle
    else
      ele2 => ele
    endif  
    slave => pointer_to_slave(ele, 1)
    s0 = (ele%s_start + 0.5_rp * ele%value(l$)) - (slave%s_start + 0.5_rp * slave%value(l$))   ! Center to center distance
    name = 'slave_drift_' // int_str(slave%ix_branch) // '_' // int_str(slave%ix_ele)
    line = 'superimpose, element = ' // trim(ele2%name) // ', ref = ' // trim(name) // ', offset = ' // re_str(s0)
    call write_lat_line (line, iu, .true.)
    cycle
  endif

  !

  call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles)
  if (.not. has_been_added) cycle
  
  ! Overlays, rampers, and groups

  if (ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$) then
    select case (ele%key)
    case (overlay$);  write (line, '(2a)') trim(ele%name), ': overlay = {'
    case (group$);    write (line, '(2a)') trim(ele%name), ': group = {'
    case (ramper$);   write (line, '(2a)') trim(ele%name), ': ramper = {'
    end select


    if (ele%key == ramper$) then
      do j = 1, size(ele%control%ramp)
        rmp => ele%control%ramp(j)
        if (j /= 1) line = trim(line) // ','
        line = trim(line) // trim(rmp%slave_name) // ' [' // trim(rmp%attribute) // ']'

        if (allocated(rmp%stack)) then
          call split_expression_string(expression_stack_to_string(rmp%stack), 120, -min(len_trim(line), 60), list)
          write (line, '(3a)') trim(line), ': ', trim(list(1))
          if (size(list) > 1) call write_lat_line(line, iu, .false., .false.)
          do ixs = 2, size(list)
            line = trim(line) // list(ixs)
            call write_lat_line(line, iu, .false., (ixs == size(list)))
          enddo
        else
          if (j > 1) then
            if (all(rmp%y_knot == ele%control%ramp(j-1)%y_knot)) cycle
          endif
          write (line, '(1000a)') trim(line), ':{', (re_str(rmp%y_knot(ix)), ', ', ix = 1, size(rmp%y_knot))
          n = len_trim(line)
          line(n:) = '}'
          call write_lat_line(line, iu, .false.)
        endif
      enddo

    else
      j_loop: do j = 1, ele%n_slave
        slave => pointer_to_slave(ele, j, ctl)
        ! do not use slaves w/ duplicate name & attribute
        do k = 1, j-1 
          slave2 => pointer_to_slave(ele, k, ctl2)
          if (slave2%name == slave%name .and. ctl2%attribute == ctl%attribute) cycle j_loop
        enddo
        ! Now write the slave info
        if (j == 1) then
          write (line, '(3a)') trim(line), trim(slave%name)
        else
          write (line, '(3a)') trim(line), ', ', trim(slave%name)
        endif
        name = ctl%attribute  
        if (name /= ele%control%var(1)%name) line = trim(line) // '[' // trim(name) // ']'

        if (allocated(ctl%stack)) then
          call split_expression_string(expression_stack_to_string(ctl%stack), 120, -min(len_trim(line), 60), list)
          if (size(list) /= 1 .or. list(1) /= ele%control%var(1)%name) then
            write (line, '(3a)') trim(line), ': ', trim(list(1))
            if (size(list) > 1) call write_lat_line(line, iu, .false., .false.)
            do ixs = 2, size(list)
              line = trim(line) // ' ' // list(ixs)
              call write_lat_line(line, iu, .false., (ixs == size(list)))
            enddo
          endif
        else
          if (j > 1) then
            if (all(ctl%y_knot == ctl2%y_knot)) cycle
          endif
          write (line, '(1000a)') trim(line), ':{', (re_str(ctl%y_knot(ix)), ', ', ix = 1, size(ctl%y_knot))
          n = len_trim(line)
          line(n:) = '}'
        endif

        if (len_trim(line) > 1000) call write_lat_line(line, iu, .false.)
      enddo j_loop
    endif

    line = trim(line) // '}, var = {' // ele%control%var(1)%name

    do j = 2, size(ele%control%var)
      line = trim(line) // ', ' // ele%control%var(j)%name
    enddo

    line = trim(line) // '}'

    if (allocated(ele%control%x_knot)) then
      write (line, '(1000a)') trim(line), ', x_knot = {', (re_str(ele%control%x_knot(ix)), ', ', &
                                                                      ix = 1, size(ele%control%x_knot))
      n = len_trim(line)
      line(n:) = '}'
    endif

    do j = 1, size(ele%control%var)
      if (ele%control%var(j)%value /= 0) then
        line = trim(line) // ', ' // trim(ele%control%var(j)%name) // ' = ' // re_str(ele%control%var(j)%value)
      endif
      if (ele%control%var(j)%old_value /= 0) then
        line = trim(line) // ', old_' // trim(ele%control%var(j)%name) // ' = ' // re_str(ele%control%var(j)%value)
      endif
    enddo

    if (nint(ele%value(interpolation$)) == linear$) line = trim(line) // ', interpolation = linear'
    if (ele%key /= ramper$ .and. is_false(ele%value(gang$))) line = trim(line) // ', gang = False'
    if (ele%type /= ' ') line = trim(line) // ', type = "' // trim(ele%type) // '"'
    if (ele%alias /= ' ') line = trim(line) // ', alias = "' // trim(ele%alias) // '"'
    if (associated(ele%descrip)) line = trim(line) // ', descrip = "' // trim(ele%descrip) // '"'
    if (has_attribute (ele, 'IS_ON') .and. .not. ele%is_on) write (line, '(2a)') trim(line), ', is_on = F'
    call write_lat_line (line, iu, .true.)
    cycle
  endif

enddo

!----------------------------------------------------------
! Lattice Layout...

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)') '! Lattice lines'
write (iu, '(a)')

call multipass_region_info(lat, mult_lat, m_info)

! Each 1st pass region is now a valid multipass line.
! Write out this info.

if (size(m_info%lord) /= 0) then
  write (iu, '(a)')
  write (iu, '(a)') '!-------------------------------------------------------'

  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    mult_ele => mult_lat%branch(ib)%ele

    in_multi_region = .false.

    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      ix_pass = m_info%branch(ib)%ele(ie)%ix_pass
      if (ix_pass /= 1) cycle 

      if (mult_ele(ie)%region_start_pt) then
        if (in_multi_region) then
          call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #1! PLEASE REPORT THIS!')
        endif
        in_multi_region = .true.
        ix_r = mult_ele(ie)%ix_region
        write (iu, '(a)')
        write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ': line[multipass] = ('
      endif

      if (mult_ele(ie)%ix_region /= ix_r) then
        call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
      endif

      call write_line_element (line, iu, ele, lat)

      if (mult_ele(ie)%region_stop_pt) then
        line = line(:len_trim(line)-1) // ')'
        call write_lat_line (line, iu, .true.)
        in_multi_region = .false.
      endif
    enddo

    if (in_multi_region) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
    endif

  enddo  ! ib branch loop
endif

! Lines for all the branches.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  write (iu, '(a)')
  name = branch%name
  if (name == '') name = 'lat_line'
  line = trim(name) // ': line = ('

  in_multi_region = .false.
  do ie = 1, branch%n_ele_track
    e_info => m_info%branch(ib)%ele(ie)
    ele => branch%ele(ie)
    if (ie == ele%branch%n_ele_track .and. ele%name == 'END' .and. ele%key == marker$) cycle

    if (.not. e_info%multipass) then
      call write_line_element (line, iu, ele, lat)
      cycle
    endif

    ix_lord = e_info%ix_lord(1)
    ix_super = e_info%ix_super(1)
    ie1 = m_info%lord(ix_lord)%slave(1,ix_super)%ele%ix_ele
    ib1 = m_info%lord(ix_lord)%slave(1,ix_super)%ele%ix_branch
    m_ele => mult_lat%branch(ib1)%ele(ie1)
    ix_r = m_ele%ix_region

    ! If entering new multipass region
    if (.not. in_multi_region) then
      in_multi_region = .true.
      if (m_ele%region_start_pt) then
        write (line, '(2a, i2.2, a)') trim(line), ' multi_line_', ix_r, ','
        look_for = 'stop'
      else
        write (line, '(2a, i2.2, a)') trim(line), ' -multi_line_', ix_r, ','
        look_for = 'start'
      endif
    endif

    if (look_for == 'start' .and. m_ele%region_start_pt .or. &
        look_for == 'stop' .and. m_ele%region_stop_pt) then 
      in_multi_region = .false.
    endif

  enddo

  line = line(:len_trim(line)-1) // ')'
  call write_lat_line (line, iu, .true.)

  ! Branch line info

  if (ib == 0) cycle

  write (iu, '(a)')
  write (iu, '(3a)') trim(branch%name), '[geometry] = ', trim(geometry_name(branch%param%geometry))
  if (branch%param%default_tracking_species /= ref_particle$) write (iu, '(3a)') trim(branch%name), &
                        '[default_tracking_species] = ', trim(species_name(branch%param%default_tracking_species))
  if (.not. branch%param%live_branch) write (iu, '(2a)') trim(branch%name), '[live_branch] = F'

  ele0 => branch%ele(0)

  if (ele0%floor%r(1) /= 0)     write (iu, '(3a)') trim(branch%name), '[x_position]     = ', re_str(ele0%floor%r(1))
  if (ele0%floor%r(2) /= 0)     write (iu, '(3a)') trim(branch%name), '[y_position]     = ', re_str(ele0%floor%r(2))
  if (ele0%floor%r(3) /= 0)     write (iu, '(3a)') trim(branch%name), '[z_position]     = ', re_str(ele0%floor%r(3))
  if (ele0%floor%theta /= 0)    write (iu, '(3a)') trim(branch%name), '[theta_position] = ', re_str(ele0%floor%theta)
  if (ele0%floor%phi /= 0)      write (iu, '(3a)') trim(branch%name), '[phi_position]   = ', re_str(ele0%floor%phi)
  if (ele0%floor%psi /= 0)      write (iu, '(3a)') trim(branch%name), '[psi_position]   = ', re_str(ele0%floor%psi)

  if (ele0%s /= 0)              write (iu, '(3a)') trim(branch%name), '[s]              = ', re_str(ele0%s)
  if (ele0%ref_time /= 0)       write (iu, '(3a)') trim(branch%name), '[ref_time]       = ', re_str(ele0%ref_time)
  if (branch%param%n_part /= 0) write (iu, '(3a)') trim(branch%name), '[n_part]         = ', re_str(lat%param%n_part)

  if (ele0%a%beta /= 0)         write (iu, '(3a)') trim(branch%name), '[beta_a]       = ', re_str(ele0%a%beta)
  if (ele0%a%alpha /= 0)        write (iu, '(3a)') trim(branch%name), '[alpha_a]      = ', re_str(ele0%a%alpha)
  if (ele0%a%phi /= 0)          write (iu, '(3a)') trim(branch%name), '[phi_a]        = ', re_str(ele0%a%phi)
  if (ele0%x%eta /= 0)          write (iu, '(3a)') trim(branch%name), '[eta_x]        = ', re_str(ele0%x%eta)
  if (ele0%x%etap /= 0)         write (iu, '(3a)') trim(branch%name), '[etap_x]       = ', re_str(ele0%x%etap)
  if (ele0%b%beta /= 0)         write (iu, '(3a)') trim(branch%name), '[beta_b]       = ', re_str(ele0%b%beta)
  if (ele0%b%alpha /= 0)        write (iu, '(3a)') trim(branch%name), '[alpha_b]      = ', re_str(ele0%b%alpha)
  if (ele0%b%phi /= 0)          write (iu, '(3a)') trim(branch%name), '[phi_b]        = ', re_str(ele0%b%phi)
  if (ele0%y%eta /= 0)          write (iu, '(3a)') trim(branch%name), '[eta_y]        = ', re_str(ele0%y%eta)
  if (ele0%y%etap /= 0)         write (iu, '(3a)') trim(branch%name), '[etap_y]       = ', re_str(ele0%y%etap)
  if (ele0%c_mat(1,1) /= 0)     write (iu, '(3a)') trim(branch%name), '[cmat_11]      = ', re_str(ele0%c_mat(1,1))
  if (ele0%c_mat(1,2) /= 0)     write (iu, '(3a)') trim(branch%name), '[cmat_12]      = ', re_str(ele0%c_mat(1,2))
  if (ele0%c_mat(2,1) /= 0)     write (iu, '(3a)') trim(branch%name), '[cmat_21]      = ', re_str(ele0%c_mat(2,1))
  if (ele0%c_mat(2,2) /= 0)     write (iu, '(3a)') trim(branch%name), '[cmat_22]      = ', re_str(ele0%c_mat(2,2))
  if (ele0%a%dbeta_dpz /= 0)    write (iu, '(3a)') trim(branch%name), '[dbeta_dpz_a]  = ', re_str(ele0%a%dbeta_dpz)
  if (ele0%b%dbeta_dpz /= 0)    write (iu, '(3a)') trim(branch%name), '[dbeta_dpz_b]  = ', re_str(ele0%b%dbeta_dpz)
  if (ele0%a%dalpha_dpz /= 0)   write (iu, '(3a)') trim(branch%name), '[dalpha_dpz_a] = ', re_str(ele0%a%dalpha_dpz)
  if (ele0%b%dalpha_dpz /= 0)   write (iu, '(3a)') trim(branch%name), '[dalpha_dpz_b] = ', re_str(ele0%b%dalpha_dpz)

  if (is_false(ele0%value(inherit_from_fork$))) write (iu, '(3a)') trim(branch%name), '[particle] = ', trim(species_name(branch%param%particle))
  write (iu, '(3a)') trim(branch%name), '[p0c]      = ', re_str(ele0%value(p0c_start$))
  if (branch%ix_from_branch >= 0) write (iu, '(2a, l1)') trim(branch%name), '[inherit_from_fork]      = ', is_true(ele0%value(inherit_from_fork$))

enddo

! Use line

line = 'use'
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_from_branch > -1) cycle
  name = branch%name
  if (name == '') name = 'lat_line'
  line = trim(line) // ', ' // name
enddo

write (iu, '(a)')
write (iu, '(a)') trim(line)

! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

have_expand_lattice_line = .false.
do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  if (ele%slave_status == super_slave$) cycle

  if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
    if (ele%value(phi0_multipass$) /= 0) then
      if (.not. have_expand_lattice_line) call write_expand_lat_header
      write (iu, '(3a)') trim(ele%name), '[phi0_multipass] = ', re_str(ele%value(phi0_multipass$))
    endif
  endif

  if (ele%slave_status == multipass_slave$) then
    multi_lord => pointer_to_multipass_lord (ele, ix_pass) 
    slave1 => pointer_to_slave(multi_lord, 1)
    if (ele%space_charge_method /= slave1%space_charge_method) then
      if (.not. have_expand_lattice_line) call write_expand_lat_header
      write (iu, '(3a)') trim(ele%name), '[space_charge_method] = ', space_charge_method_name(ele%space_charge_method)
    endif
  endif
enddo

! If there are lattice elements with duplicate names but differing parameters then
! Write the differences.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 1, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%key == marker$ .and. ele%name == 'END') cycle
    if (ele%slave_status == super_slave$) cycle
    if (ele%slave_status == multipass_slave$) cycle
    call eles_with_same_name_handler(ele, named_eles, an_indexx, names, n_names, order)
  enddo
enddo

! cleanup

close(iu)
deallocate (names, an_indexx)
deallocate (mult_lat%branch)

if (present(err)) err = .false.

!--------------------------------------------------------------------------------
contains

subroutine eles_with_same_name_handler(ele, named_eles, an_indexx, names, n_names, order)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_struct), pointer :: lat
type (ele_pointer_struct) :: named_eles(:)
type (lat_ele_order_struct) order

real(rp), pointer :: a0(:), b0(:), ksl0(:), a(:), b(:), ksl(:)
real(rp), target :: az(0:n_pole_maxx) = 0, bz(0:n_pole_maxx) = 0
character(40), allocatable :: names(:)
integer, allocatable :: an_indexx(:)
integer n_names, ix_match
integer i, iv

!

lat => ele%branch%lat
if (ele%slave_status == multipass_slave$) return
call find_index (ele%name, names, an_indexx, n_names, ix_match)
ele0 => named_eles(ix_match)%ele   ! Element with this name whose attributes were written to the lattice file.
if (ele%ix_ele == ele0%ix_ele .and. ele%ix_branch == ele0%ix_branch) return

do iv = 1, num_ele_attrib$
  if (ele%value(iv) == ele0%value(iv)) cycle
  info = attribute_info(ele, iv)
  if (info%state /= is_free$ .and. info%state /= quasi_free$) cycle
  if (info%state == quasi_free$) then
    if (.not. attribute_free(ele, info%name, .false.)) cycle
  endif
  ! Have a differing attribute
  call write_this_differing_attrib(iu, ele, attribute_name(ele, iv), ele%value(iv), order)
enddo

if (associated(ele%a_pole) .or. associated(ele0%a_pole)) then
  call pointer_to_ele_multipole(ele0, a0, b0, ksl0, magnetic$)
  if (.not. associated(a0)) a0 => az
  if (.not. associated(b0)) b0 => bz
  call pointer_to_ele_multipole(ele, a, b, ksl, magnetic$)
  if (.not. associated(a)) a => az
  if (.not. associated(b)) b => bz

  if (ele%key == multipole$) then
    do i = 0, n_pole_maxx
      if (a(i) /= a0(i)) call write_this_differing_attrib(iu, ele, 'k' // int_str(i) // 'l', a(i), order)
      if (b(i) /= b0(i)) call write_this_differing_attrib(iu, ele, 't' // int_str(i), b(i), order)
      if (ksl(i) /= ksl0(i)) call write_this_differing_attrib(iu, ele, 'k' // int_str(i) // 'sl', ksl(i), order)
    enddo
  else
    do i = 0, n_pole_maxx
      if (a(i) /= a0(i)) call write_this_differing_attrib(iu, ele, 'a' // int_str(i), a(i), order)
      if (b(i) /= b0(i)) call write_this_differing_attrib(iu, ele, 'b' // int_str(i), b(i), order)
    enddo
  endif
endif

if (associated(ele%b_pole_elec) .or. associated(ele0%b_pole_elec)) then
  call pointer_to_ele_multipole(ele0, a0, b0, ksl0, electric$)
  if (.not. associated(a0)) a0 => az
  if (.not. associated(b0)) b0 => bz
  call pointer_to_ele_multipole(ele, a, b, ksl, electric$)
  if (.not. associated(a)) a => az
  if (.not. associated(b)) b => bz

  do i = 0, n_pole_maxx
    if (a(i) /= a0(i)) call write_this_differing_attrib(iu, ele, 'a' // int_str(i) // '_elec', a(i), order)
    if (b(i) /= b0(i)) call write_this_differing_attrib(iu, ele, 'b' // int_str(i) // '_elec', b(i), order)
  enddo
endif

end subroutine eles_with_same_name_handler

!--------------------------------------------------------------------------------
! contains

subroutine write_this_differing_attrib(iu, ele, attrib_name, value, order)

type (ele_struct) ele
type (lat_ele_order_struct) order

integer iu
real(rp) value
character(*) attrib_name

!

if (.not. have_expand_lattice_line) call write_expand_lat_header

write (iu, '(5a)') trim(ele_unique_name(ele, order)), '[', trim(attrib_name), '] = ', re_str(value)

end subroutine write_this_differing_attrib

!--------------------------------------------------------------------------------
! contains

subroutine write_if_real_param_changed (param_now, param_default, param_name)

real(rp) param_now, param_default
character(*) param_name

!

if (abs(param_now - param_default) <= 1e-12 * (abs(param_now) + abs(param_default))) return
write (iu, '(3a)') param_name, ' = ', re_str(param_now)

end subroutine write_if_real_param_changed

!--------------------------------------------------------------------------------
! contains

subroutine write_if_int_param_changed (param_now, param_default, param_name)

integer param_now, param_default
character(*) param_name

!

if (param_now == param_default) return
write (iu, '(2a, i0)') param_name, ' = ', param_now

end subroutine write_if_int_param_changed

!--------------------------------------------------------------------------------
! contains

subroutine write_if_logic_param_changed (param_now, param_default, param_name)

logical param_now, param_default
character(*) param_name

!

if (param_now .eqv. param_default) return
write (iu, '(2a, l1)') param_name, ' = ', param_now

end subroutine write_if_logic_param_changed

!--------------------------------------------------------------------------------
! contains

subroutine write_expand_lat_header ()

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)')
write (iu, '(a)') 'expand_lattice'
write (iu, '(a)')
have_expand_lattice_line = .true.

end subroutine write_expand_lat_header

!--------------------------------------------------------------------------------
! contains

subroutine form_this_fieldmap_name(string, field_type, ele, ix_map, output_form)

type (ele_struct) ele
integer ix_map
integer, optional :: output_form
character(*) string, field_type
character(40) name

!

call str_downcase(name, ele%name)
if (ix_map == 1) then
  write (string, '(2a)') trim(name), field_type
else
  write (string, '(2a, i0, a)') trim(name), '_', ix_map, field_type
endif

if (integer_option(binary$, output_form) == binary$) string = trim(string) // '.h5'

end subroutine form_this_fieldmap_name

!--------------------------------------------------------------------------------
! contains

subroutine write_map_coef (tag, array)

real(rp) :: array(:)
character(*) tag
integer j, k0, k, n, n5

!

write (iu2, '(a)') ','
write (iu2, '(2x, 2a)', advance = 'NO') tag, ' = ('

n = size(array)
n5 = (n - 1) / 5

do j = 1, n5
  k0 = 5 * (j - 1)
  write (iu2, '(5(es14.6, a))') (array(k0+k), ',', k = 1, 5) 
  write (iu2, '(15x)', advance = 'NO')
enddo

if (n5 == 0) write (iu2, '(15x)', advance = 'NO')
write (iu2, '(5(es14.6, a))', advance = 'NO') (array(k), ',', k = 5*n5+1, n-1), array(n), ')' 

end subroutine write_map_coef

!--------------------------------------------------------------------------------
! contains

subroutine write_this_cartesian_map (ct_map, ele, iu9, line)

type (cartesian_map_struct) :: ct_map
type (ele_struct) ele
integer j, iu9
character(*), optional :: line

!

write (iu9, '(a)') '{'
if (ct_map%master_parameter > 0) write (iu9, '(2x, 3a)') &
                        'master_parameter  = ', trim(attribute_name(ele, ct_map%master_parameter)), ','
write (iu9, '(2x, 3a)') 'field_scale       = ', re_str(ct_map%field_scale), ','
write (iu9, '(2x, 4a)') 'r0                = ', trim(array_re_str(ct_map%r0)), ','
write (iu9, '(2x, 3a)') 'ele_anchor_pt     = ', trim(anchor_pt_name(ct_map%ele_anchor_pt)), ','
write (iu9, '(2x, 3a)') 'field_type        = ', trim(em_field_type_name(ct_map%field_type)), ','

do j = 1, size(ct_map%ptr%term)
  ct_term => ct_map%ptr%term(j)
  last = '}, &'
  if (j == size(ct_map%ptr%term)) last = '} &'
  select case (ct_term%family)
  case (family_x$)
    name = 'X'
  case (family_y$)
    name = 'Y'
  case (family_qu$)
    name = 'QU'
  case (family_sq$)
    name = 'SQ'
  end select
  write (iu9, '(17a)') '  term = {', re_str(ct_term%coef), ', ', &
    re_str(ct_term%kx), ', ', re_str(ct_term%ky), ', ', re_str(ct_term%kz), &
    ', ', re_str(ct_term%x0), ', ', re_str(ct_term%y0), ', ', re_str(ct_term%phi_z), ', ', trim(name), trim(last)
enddo

! present(line) = T when single file is being constructed

if (present(line)) then
  line = '}'
else
  write (iu9, '(a)') '}'
endif

end subroutine write_this_cartesian_map

!--------------------------------------------------------------------------------
! contains

subroutine write_this_cylindrical_map (cl_map, ele, iu9, line)

type (cylindrical_map_struct) :: cl_map
type (ele_struct) ele
integer iu9
character(*), optional :: line

!

write (iu9, '(a)') '{'
if (cl_map%master_parameter > 0) write (iu9, '(2x, 3a)') &
                                          'master_parameter  = ', trim(attribute_name(ele, cl_map%master_parameter)), ','
write (iu9, '(2x, 3a)')       'field_scale       = ', re_str(cl_map%field_scale), ','
write (iu9, '(2x, 3a)')       'ele_anchor_pt     = ', trim(anchor_pt_name(cl_map%ele_anchor_pt)), ','
write (iu9, '(2x, a, i0, a)') 'm                 = ', cl_map%m, ','
write (iu9, '(2x, a, i0, a)') 'harmonic          = ', cl_map%harmonic, ','
write (iu9, '(2x, 3a)')       'dz                = ', re_str(cl_map%dz), ','
write (iu9, '(2x, 4a)')       'r0                = ', trim(array_re_str(cl_map%r0)), ','
write (iu9, '(2x, 3a)')       'phi0_fieldmap     = ', re_str(cl_map%phi0_fieldmap), ','
write (iu9, '(2x, 3a)', advance = 'NO') 'theta0_azimuth      = ', re_str(cl_map%theta0_azimuth)

if (any(real(cl_map%ptr%term%e_coef) /= 0)) call write_map_coef ('E_coef_re', real(cl_map%ptr%term%e_coef))
if (any(aimag(cl_map%ptr%term%e_coef) /= 0)) call write_map_coef ('E_coef_im', aimag(cl_map%ptr%term%e_coef))
if (any(real(cl_map%ptr%term%b_coef) /= 0)) call write_map_coef ('B_coef_re', real(cl_map%ptr%term%b_coef))
if (any(aimag(cl_map%ptr%term%b_coef) /= 0)) call write_map_coef ('B_coef_im', aimag(cl_map%ptr%term%b_coef))

! present(line) = T when single file is being constructed

if (present(line)) then
  line = '}'
else
  write (iu9, '(a)') '}'
endif

end subroutine write_this_cylindrical_map

!--------------------------------------------------------------------------------
! contains

subroutine write_this_grid_fieldmap (g_field, ele, iu9, line)

type (grid_field_struct) :: g_field
type (ele_struct) ele
integer iu9
character(*), optional :: line
character(200) string

!

write (iu9, '(a)') '{'
n = grid_field_dimension(g_field%geometry)
write (iu9, '(2x, 3a)')       'geometry            = ', trim(grid_field_geometry_name(g_field%geometry)), ','
if (g_field%master_parameter > 0) write (iu9, '(2x, 3a)') &
                              'master_parameter    = ', trim(attribute_name(ele, g_field%master_parameter)), ','
write (iu9, '(2x, 3a)')       'field_scale         = ', re_str(g_field%field_scale), ','
write (iu9, '(2x, 3a)')       'ele_anchor_pt       = ', trim(anchor_pt_name(g_field%ele_anchor_pt)), ','
write (iu9, '(2x, 3a)')       'field_type          = ', trim(em_field_type_name(g_field%field_type)), ','
write (iu9, '(2x, a, i0, a)') 'interpolation_order = ', g_field%interpolation_order, ','
write (iu9, '(2x, a, i0, a)') 'harmonic            = ', g_field%harmonic, ','
write (iu9, '(2x, 3a)')       'phi0_fieldmap       = ', re_str(g_field%phi0_fieldmap), ','
write (iu9, '(2x, 4a)')       'dr                  = ', trim(array_re_str(g_field%dr(1:n))), ','
write (iu9, '(2x, 4a)')       'r0                  = ', trim(array_re_str(g_field%r0)), ','
write (iu9, '(2x, a, l1, a)') 'curved_ref_frame    = ', g_field%curved_ref_frame, ','

end_str = '),'

do id1 = lbound(g_field%ptr%pt, 1), ubound(g_field%ptr%pt, 1)          
do id2 = lbound(g_field%ptr%pt, 2), ubound(g_field%ptr%pt, 2)
do id3 = lbound(g_field%ptr%pt, 3), ubound(g_field%ptr%pt, 3)

  if (all([id1, id2, id3] == ubound(g_field%ptr%pt))) end_str = ') &'

  select case (grid_field_dimension(g_field%geometry))
  case (1)
    write (string, '(2x, a, i0, 13a)') 'pt(', id1, ') = ('
  case (2)
    write (string, '(2x, 2(a, i0), 13a)') 'pt(', id1, ',', id2, ') = ('
  case (3)
    write (string, '(2x, 3(a, i0), 13a)') 'pt(', id1, ',', id2, ',', id3, ') = ('
  end select

  select case (g_field%field_type)
  case (mixed$)
    write (iu9, '(2x, a, 13a)') trim(string), &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(1))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(2))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(3))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(1))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(2))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(3))), end_str
  case (electric$)
    write (iu9, '(2x, a, 13a)') trim(string), &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(1))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(2))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(3))), end_str
  case (magnetic$)
    write (iu9, '(2x, a, 13a)') trim(string), &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(1))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(2))), ',', &
                                       trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(3))), end_str
  end select

enddo
enddo
enddo

! present(line) = T when single file is being constructed

if (present(line)) then
  line = '}'
else
  write (iu9, '(a)') '}'
endif

end subroutine write_this_grid_fieldmap

!--------------------------------------------------------------------------------
! contains

subroutine write_this_gen_grad_map_map (gg_map, ele, iu9, line)

type (gen_grad_map_struct), target :: gg_map
type (gen_grad1_struct), pointer :: gg
type (ele_struct) ele
integer iu9, ig, iz, n, id
character(*), optional :: line
character(40) fmt

!

write (iu9, '(a)') '{'
if (gg_map%master_parameter > 0) write (iu9, '(2x, 3a)') &
                              'master_parameter   = ', trim(attribute_name(ele, gg_map%master_parameter)), ','
write (iu9, '(2x, 3a)')       'field_scale        = ', re_str(gg_map%field_scale), ','
write (iu9, '(2x, 3a)')       'ele_anchor_pt      = ', trim(anchor_pt_name(gg_map%ele_anchor_pt)), ','
write (iu9, '(2x, 3a)')       'field_type         = ', trim(em_field_type_name(gg_map%field_type)), ','
write (iu9, '(2x, 3a)')       'dz                 = ', re_str(gg_map%dz), ','
write (iu9, '(2x, 4a)')       'r0                 = ', trim(array_re_str(gg_map%r0)), ','
write (iu9, '(2x, a, l1, a)') 'curved_ref_frame   = ', gg_map%curved_ref_frame, ','

do ig = 1, size(gg_map%gg)
  gg => gg_map%gg(ig)
  write (iu9, '(2x, a)') 'curve = {'
  write (iu9, '(4x, a, i0, a)') 'm = ', gg%m, ','
  write (iu9, '(4x, 3a)')       'kind = ', trim(expression_op_name(gg%sincos)), ','
  write (iu9, '(4x, a)')        'derivs = {'


  n = gg%n_deriv_max
  write (fmt, '(a, i0, a)') '(f13.5, a, ', n+1, 'es20.12, a)' 
  do iz = gg_map%iz0, gg_map%iz1-1
    write (iu9, fmt) iz*gg_map%dz, ':', gg%deriv(iz,0:n), ','
  enddo
  write (iu9, fmt) gg_map%iz1*gg_map%dz, ':', gg%deriv(gg_map%iz1,0:n)

  write (iu9, '(a)')   '    }'
  if (ig == size(gg_map%gg)) then
    write (iu9, '(a)') '  }'
  else
    write (iu9, '(a)') '  },'
  endif
enddo

! present(line) = T when single file is being constructed

if (present(line)) then
  line = '}'
else
  write (iu9, '(a)') '}'
endif

end subroutine write_this_gen_grad_map_map

end subroutine write_bmad_lattice_file
