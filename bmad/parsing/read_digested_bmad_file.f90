!+
! Subroutine read_digested_bmad_file (digested_file, lat, inc_version, err_flag, parser_calling, lat_files)
!
! Subroutine to read in a digested file. The subroutine will check that the version of the digested file 
! is up to date and that the digested file is current with respect to the original BMAD files that were used. 
! [See write_digested_bmad_file.]
!
! Note: This subroutine also reads in the common structures for bmad_parser2
!
! Input:
!   digested_file -- Character(*): Name of the digested file.
!
! Output:
!   lat             -- lat_struct: Output lattice structure
!   inc_version     -- integer: bmad version number stored in the lattice file.
!                       If the file is current this number should be the same 
!                       as the global parameter bmad_inc_version$. 
!                       Set to -1 if there is a read error.
!   err_flag        -- logical, optional: Set True if there is an error. False otherwise.
!   parser_calling  -- logical, optional: Is this routine being called from a parser routine (like bmad_parser)?
!                       Default is False. This argument determines what are considered errors. For example, a
!                       moved digested file is considered an error if this routine is called from a parser but
!                       not otherwise. The reason for this dichotomy is that a parser is able to reread the
!                       original lattice file.
!   lat_files(:)    -- character(400), optional, allocatable: List of Bmad lattice files that defined this lattice.
!-

subroutine read_digested_bmad_file (digested_file, lat, inc_version, err_flag, parser_calling, lat_files)

use ptc_interface_mod, dummy => read_digested_bmad_file
use bmad_parser_mod, dummy2 => read_digested_bmad_file
use wall3d_mod, dummy3 => read_digested_bmad_file

implicit none

type (lat_struct), target, intent(inout) :: lat
type (branch_struct), pointer :: branch
type (extra_parsing_info_struct) :: extra
type (bmad_common_struct) bmad_com_read
type (space_charge_common_struct) space_charge_com_read
type (ptc_common_struct) ptc_com_read
type (control_struct), pointer :: ctl

real(rp) value(num_ele_attrib$)

integer inc_version, d_unit, n_files, file_version, i, j, k, ix, ix_value(num_ele_attrib$)
integer stat_b(13), stat_b2, stat_b8, stat_b10, n_branch, n, nk, control_type, coupler_at
integer ierr, stat, ios, ios2, ios3, ios4, n_wall_section, garbage, j1, j2, io_err_level, n_custom, n_print
integer, allocatable :: index_list(:)

character(*) digested_file
character(*), optional, allocatable :: lat_files(:)
character(400) fname_read, fname_versionless, fname_full
character(400) input_file_name, full_digested_file, digested_prefix_in, digested_prefix_out
character(100), allocatable :: name_list(:)
character(*), parameter :: r_name = 'read_digested_bmad_file'

logical, optional :: err_flag, parser_calling
logical is_ok, parser_call
logical found_it, mode3, error, is_match, err, err_found

! Init all elements in lat

call init_bmad()

if (present(err_flag)) err_flag = .true.
err_found = .false.
parser_call = logic_option(.false., parser_calling)

if (parser_call) then
  io_err_level = s_info$
else
  io_err_level = s_error$
endif

if (digested_file == '') then ! For some reason the inquire statement will not catch this error.
  call out_io (io_err_level, r_name, 'DIGESTED FILE NAME IS BLANK!')
  return
endif

call init_lat (lat)

! Read the digested file.

d_unit = lunget()
inc_version = -1
lat%n_ele_track = 0

call fullfilename (digested_file, fname_full)
inquire (file = fname_full, name = full_digested_file)
call simplify_path (full_digested_file, full_digested_file)
open (unit = d_unit, file = full_digested_file, status = 'old',  form = 'unformatted', action = 'READ', err = 9000)

read (d_unit, err = 9010, end = 9010) n_files, file_version
if (present(lat_files)) call re_allocate (lat_files, n_files)

! Version is old

if (file_version < bmad_inc_version$) then
  call out_io (io_err_level, r_name, ['DIGESTED FILE VERSION OUT OF DATE \i0\ > \i0\ ' ],  &
                                i_array = [bmad_inc_version$, file_version ])
  close (d_unit)
  return
endif

if (file_version > bmad_inc_version$) then
  call out_io (io_err_level, r_name, &
       'DIGESTED FILE HAS VERSION: \i0\ ', &
       'GREATER THAN VERSION OF THIS PROGRAM: \i0\ ', &
       'WILL NOT USE THE DIGESTED FILE. YOU SHOULD RECOMPILE THIS PROGRAM.', &
       i_array = [file_version, bmad_inc_version$])
  close (d_unit)
  return
endif

! if the digested file is out of date then we still read in the file since
! we can possibly reuse the taylor series.

do i = 1, n_files
  stat_b = 0

  read (d_unit, err = 9020, end = 9020) fname_read, stat_b2, stat_b8, stat_b10
  if (present(lat_files)) lat_files(i) = fname_read
  if (.not. parser_call) cycle  ! Do not need to check if file moved in this case.

  ! Cannot use full file name to check if this is the original digested file since
  ! the path may change depending upon what system the program is running on and how
  ! things are mounted. So use stat() instead

  if (fname_read(1:10) == '!DIGESTED:') then
    fname_read = fname_read(11:)
    ierr = stat(full_digested_file, stat_b)
    ! Time stamp in file is created while file is being written to so is not accurate.
    is_match = (stat_b2 == stat_b(2))            !!!! .and. (stat_b10 == stat_b(10))
    j1 = len_trim(fname_read)
    j2 = len_trim(full_digested_file)
    do j = 0, min(j1, j2) - 1
      if (fname_read(j1-j:j1-j) /= full_digested_file(j2-j:j2-j)) exit
    enddo
    digested_prefix_in = fname_read(1:j1-j)
    digested_prefix_out = full_digested_file(1:j2-j)
    if (.not. is_match) then
      call out_io(s_info$, r_name, ' NOTE: MOVED DIGESTED FILE.')
      close (d_unit)
      return
    endif
    cycle
  endif

  call simplify_path (fname_read, fname_read)

  is_ok = .true.
  if (digested_prefix_in /= '') then
    if (index(fname_read, trim(digested_prefix_in)) == 1) then
      ix = len_trim(digested_prefix_in)
      fname_read = fname_read(ix+1:)
    else
      is_ok = .false.
    endif
  endif
  if (digested_prefix_out /= '') then
    fname_read = trim(digested_prefix_out) // trim(fname_read)
  endif
  ierr = stat(fname_read, stat_b)
  fname_versionless = fname_read
  is_match = (stat_b2 == stat_b(2)) .and. (stat_b10 == stat_b(10))

  inquire (file = fname_versionless, exist = found_it, name = fname_full)
  call simplify_path (fname_full, fname_full)
  if (.not. found_it .or. fname_read /= fname_full .or. .not. is_match) then
    call out_io(s_info$, r_name, 'NOTE: DIGESTED FILE OUT OF DATE.')
    close (d_unit)
    return
  endif
enddo

! we read (and write) the lat in pieces since it is
! too big to write in one piece

read (d_unit, err = 9030, end = 9030) lat%use_name, lat%machine, lat%lattice, lat%input_file_name, lat%title
read (d_unit, err = 9030, end = 9030) lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track
read (d_unit, err = 9030, end = 9030) lat%n_ele_track, lat%n_ele_max, lat%lord_state, lat%n_control_max, lat%n_ic_max
read (d_unit, err = 9030, end = 9030) lat%input_taylor_order, lat%photon_type, lat%ramper_slave_bookkeeping
read (d_unit, err = 9070, end = 9070) n_branch, lat%pre_tracker, n_custom, n_print

! Different compilers (EG ifort and gfortran) will produce different binary formats. 
! As a double check, check the version number again.

if (lat%version /= bmad_inc_version$) then
  call out_io (io_err_level, r_name, 'DIGESTED FILE BINARY FORMAT IS WRONG.', &
         '[CAN HAPPEN IF THE DIGESTED FILE IS CREATED WITH A PROGRAM COMPILED WITH A DIFFERENT COMPILER.]')
  close (d_unit)
  return
endif

! Global custom & print statements

if (n_custom > -1) then
  call re_allocate(lat%custom, n_custom)
  read (d_unit, err = 9070, end = 9070) lat%custom
endif

if (n_print > -1) then
  call re_allocate(lat%print_str, n_print)
  read (d_unit, err = 9070, end = 9070) lat%print_str
endif

! Defined constants and custom attributes

read (d_unit, err = 9035, end = 9035) n
allocate(index_list(n), name_list(n))
do i = 1, n
  read (d_unit, err = 9035, end = 9035) index_list(i), name_list(i)
enddo

read (d_unit, err = 9035, end = 9035) n
if (allocated(lat%constant)) deallocate(lat%constant)
allocate(lat%constant(n))
do i = 1, n
  read (d_unit, err = 9035, end = 9035) lat%constant(i)
enddo

! Allocate lat%ele, lat%control and lat%ic arrays

call allocate_lat_ele_array(lat, lat%n_ele_max+10)
call reallocate_control (lat, lat%n_control_max+10)

! Branches

do i = 0, lat%n_ele_max
  call read_this_ele(lat%ele(i), i, error)
  if (error) return
enddo

call allocate_branch_array (lat, n_branch)  ! Initial allocation

call read_this_wall3d (lat%branch(0)%wall3d, error)
if (error) return

read (d_unit, err = 9070, end = 9070) lat%branch(0)%name

do i = 1, n_branch
  branch => lat%branch(i)
  branch%ix_branch = i
  read (d_unit, err = 9070, end = 9070) branch%param
  read (d_unit, err = 9070, end = 9070) branch%name, branch%ix_from_branch, branch%ix_to_ele, &
                 branch%ix_from_ele, branch%n_ele_track, branch%n_ele_max

  call allocate_lat_ele_array (lat, branch%n_ele_max, i)
  do j = 0, branch%n_ele_max
    call read_this_ele (branch%ele(j), j, error)
    if (error) return
  enddo

  call read_this_wall3d (branch%wall3d, error)
  if (error) return

enddo

! read the control info, etc

do i = 1, lat%n_control_max
  ctl => lat%control(i)
  read (d_unit, err = 9040, end = 9040) n, nk, ctl%value, ctl%lord, ctl%slave, ctl%ix_attrib, ctl%attribute, ctl%slave_name
  if (n > 0) then
    allocate (ctl%stack(n))
    do j = 1, n
      read (d_unit, err = 9045, end = 9045) ctl%stack(j)
    enddo
  endif

  if (nk > 0) then
    allocate (ctl%y_knot(nk))
    read (d_unit, err = 9045, end = 9045) ctl%y_knot
  endif
enddo

do i = 1, lat%n_ic_max
  read (d_unit, err = 9050, end = 9050) lat%ic(i)
enddo

read (d_unit, err = 9060, end = 9060) lat%particle_start
read (d_unit, err = 9060, end = 9060) lat%beam_init

! Read extra state info.

read (d_unit, iostat = ios) found_it
if (found_it) then
  allocate (ptc_com_read%vertical_kick, ptc_com_read%old_integrator, ptc_com_read%exact_model, &
            ptc_com_read%exact_misalign, ptc_com_read%max_fringe_order)
  read (d_unit, iostat = ios) extra
  read (d_unit, iostat = ios2) bmad_com_read
  read (d_unit, iostat = ios3) space_charge_com_read
  read (d_unit, iostat = ios4) ptc_com_read%max_fringe_order, ptc_com_read%exact_model, ptc_com_read%exact_misalign, &
          ptc_com_read%vertical_kick, ptc_com_read%cut_factor, ptc_com_read%old_integrator, &
          ptc_com_read%use_orientation_patches, ptc_com_read%print_info_messages, ptc_com_read%translate_patch_drift_time


  if (ios /= 0 .or. ios2 /= 0 .or. ios3 /= 0 .or. ios4 /= 0) then
    call out_io (io_err_level, r_name, 'ERROR READING BMAD/SPACE_CHARGE/PTC COMMON PARAMETERS')
    close (d_unit)
    return
  endif
  if (extra%max_aperture_limit_set)               bmad_com%max_aperture_limit              = bmad_com_read%max_aperture_limit
  if (extra%default_ds_step_set)                  bmad_com%default_ds_step                 = bmad_com_read%default_ds_step
  if (extra%significant_length_set)               bmad_com%significant_length              = bmad_com_read%significant_length
  if (extra%rel_tol_tracking_set)                 bmad_com%rel_tol_tracking                = bmad_com_read%rel_tol_tracking
  if (extra%abs_tol_tracking_set)                 bmad_com%abs_tol_tracking                = bmad_com_read%abs_tol_tracking
  if (extra%rel_tol_adaptive_tracking_set)        bmad_com%rel_tol_adaptive_tracking       = bmad_com_read%rel_tol_adaptive_tracking
  if (extra%abs_tol_adaptive_tracking_set)        bmad_com%abs_tol_adaptive_tracking       = bmad_com_read%abs_tol_adaptive_tracking
  if (extra%init_ds_adaptive_tracking_set)        bmad_com%init_ds_adaptive_tracking       = bmad_com_read%init_ds_adaptive_tracking
  if (extra%min_ds_adaptive_tracking_set)         bmad_com%min_ds_adaptive_tracking        = bmad_com_read%min_ds_adaptive_tracking
  if (extra%fatal_ds_adaptive_tracking_set)       bmad_com%fatal_ds_adaptive_tracking      = bmad_com_read%fatal_ds_adaptive_tracking
  if (extra%autoscale_amp_abs_tol_set)            bmad_com%autoscale_amp_abs_tol           = bmad_com_read%autoscale_amp_abs_tol
  if (extra%autoscale_amp_rel_tol_set)            bmad_com%autoscale_amp_rel_tol           = bmad_com_read%autoscale_amp_rel_tol
  if (extra%autoscale_phase_tol_set)              bmad_com%autoscale_phase_tol             = bmad_com_read%autoscale_phase_tol
  if (extra%rf_phase_below_transition_ref_set)    bmad_com%rf_phase_below_transition_ref   = bmad_com_read%rf_phase_below_transition_ref
  if (extra%electric_dipole_moment_set)           bmad_com%electric_dipole_moment          = bmad_com_read%electric_dipole_moment
  if (extra%synch_rad_scale_set)                  bmad_com%synch_rad_scale                 = bmad_com_read%synch_rad_scale
  if (extra%taylor_order_set)                     bmad_com%taylor_order                    = bmad_com_read%taylor_order
  if (extra%d_orb_set)                            bmad_com%d_orb                           = bmad_com_read%d_orb
  if (extra%default_integ_order_set)              bmad_com%default_integ_order             = bmad_com_read%default_integ_order
  if (extra%runge_kutta_order_set)                bmad_com%runge_kutta_order               = bmad_com_read%runge_kutta_order
  if (extra%sr_wakes_on_set)                      bmad_com%sr_wakes_on                     = bmad_com_read%sr_wakes_on
  if (extra%lr_wakes_on_set)                      bmad_com%lr_wakes_on                     = bmad_com_read%lr_wakes_on
  if (extra%high_energy_space_charge_on_set)      bmad_com%high_energy_space_charge_on     = bmad_com_read%high_energy_space_charge_on
  if (extra%csr_and_space_charge_on_set)          bmad_com%csr_and_space_charge_on         = bmad_com_read%csr_and_space_charge_on
  if (extra%spin_tracking_on_set)                 bmad_com%spin_tracking_on                = bmad_com_read%spin_tracking_on
  if (extra%spin_sokolov_ternov_flipping_on_set)  bmad_com%spin_sokolov_ternov_flipping_on = bmad_com_read%spin_sokolov_ternov_flipping_on
  if (extra%radiation_damping_on_set)             bmad_com%radiation_damping_on            = bmad_com_read%radiation_damping_on
  if (extra%radiation_zero_average_set)           bmad_com%radiation_zero_average          = bmad_com_read%radiation_zero_average
  if (extra%radiation_fluctuations_on_set)        bmad_com%radiation_fluctuations_on       = bmad_com_read%radiation_fluctuations_on
  if (extra%conserve_taylor_maps_set)             bmad_com%conserve_taylor_maps            = bmad_com_read%conserve_taylor_maps
  if (extra%absolute_time_tracking_set)           bmad_com%absolute_time_tracking          = bmad_com_read%absolute_time_tracking
  if (extra%absolute_time_ref_shift_set)          bmad_com%absolute_time_ref_shift         = bmad_com_read%absolute_time_ref_shift
  if (extra%convert_to_kinetic_momentum_set)      bmad_com%convert_to_kinetic_momentum     = bmad_com_read%convert_to_kinetic_momentum
  if (extra%aperture_limit_on_set)                bmad_com%aperture_limit_on               = bmad_com_read%aperture_limit_on
  if (extra%sad_eps_scale_set)                    bmad_com%sad_eps_scale                   = bmad_com_read%sad_eps_scale
  if (extra%sad_amp_max_set)                      bmad_com%sad_amp_max                     = bmad_com_read%sad_amp_max
  if (extra%sad_n_div_max_set)                    bmad_com%sad_n_div_max                   = bmad_com_read%sad_n_div_max
  if (extra%max_num_runge_kutta_step_set)         bmad_com%max_num_runge_kutta_step        = bmad_com_read%max_num_runge_kutta_step
  if (extra%debug_set)                            bmad_com%debug                           = bmad_com_read%debug

  if (extra%ds_track_step_set)                    space_charge_com%ds_track_step                    = space_charge_com_read%ds_track_step
  if (extra%dt_track_step_set)                    space_charge_com%dt_track_step                    = space_charge_com_read%dt_track_step
  if (extra%cathode_strength_cutoff_set)          space_charge_com%cathode_strength_cutoff          = space_charge_com_read%cathode_strength_cutoff
  if (extra%sc_rel_tol_tracking_set)              space_charge_com%rel_tol_tracking                 = space_charge_com_read%rel_tol_tracking
  if (extra%sc_abs_tol_tracking_set)              space_charge_com%abs_tol_tracking                 = space_charge_com_read%abs_tol_tracking
  if (extra%beam_chamber_height_set)              space_charge_com%beam_chamber_height              = space_charge_com_read%beam_chamber_height
  if (extra%lsc_sigma_cutoff_set)                 space_charge_com%lsc_sigma_cutoff                 = space_charge_com_read%lsc_sigma_cutoff
  if (extra%particle_sigma_cutoff_set)            space_charge_com%particle_sigma_cutoff            = space_charge_com_read%particle_sigma_cutoff
  if (extra%space_charge_mesh_size_set)           space_charge_com%space_charge_mesh_size           = space_charge_com_read%space_charge_mesh_size
  if (extra%csr3d_mesh_size_set)                  space_charge_com%csr3d_mesh_size                  = space_charge_com_read%csr3d_mesh_size
  if (extra%n_bin_set)                            space_charge_com%n_bin                            = space_charge_com_read%n_bin
  if (extra%particle_bin_span_set)                space_charge_com%particle_bin_span                = space_charge_com_read%particle_bin_span
  if (extra%n_shield_images_set)                  space_charge_com%n_shield_images                  = space_charge_com_read%n_shield_images
  if (extra%sc_min_in_bin_set)                    space_charge_com%sc_min_in_bin                    = space_charge_com_read%sc_min_in_bin
  if (extra%lsc_kick_transverse_dependence_set)   space_charge_com%lsc_kick_transverse_dependence   = space_charge_com_read%lsc_kick_transverse_dependence
  if (extra%sc_debug_set)                         space_charge_com%debug                            = space_charge_com_read%debug
  if (extra%diagnostic_output_file_set)           space_charge_com%diagnostic_output_file           = space_charge_com_read%diagnostic_output_file

  if (extra%use_orientation_patches_set)          ptc_com%use_orientation_patches                   = ptc_com_read%use_orientation_patches
  if (extra%print_info_messages_set)              ptc_com%print_info_messages                       = ptc_com_read%print_info_messages
  if (extra%cut_factor_set)                       ptc_com%cut_factor                                = ptc_com_read%cut_factor
  if (extra%max_fringe_order_set)                 ptc_com%max_fringe_order                          = ptc_com_read%max_fringe_order
  if (extra%vertical_kick_set)                    ptc_com%vertical_kick                             = ptc_com_read%vertical_kick
  if (extra%old_integrator_set)                   ptc_com%old_integrator                            = ptc_com_read%old_integrator
  if (extra%exact_model_set)                      ptc_com%exact_model                               = ptc_com_read%exact_model
  if (extra%exact_misalign_set)                   ptc_com%exact_misalign                            = ptc_com_read%exact_misalign
  if (extra%translate_patch_drift_time_set)       ptc_com%translate_patch_drift_time                = ptc_com_read%translate_patch_drift_time

  if (extra%undeterministic_ran_function_called) err_found = .true.  ! So lattice will be reparsed

  if (extra%ran_seed /= 0) then
    call ran_default_state (set_state = extra%ran_state) ! Get random state.
  endif
endif

! Setup any attribute aliases in the global attribute name table.
! This is done last in read_digested bmad_file so as to not to pollute the name table if 
! there is an error.

if (.not. err_found) then
  do i = 1, size(index_list)
    call set_custom_attribute_name(name_list(i), err, index_list(i))
    if (err) err_found = .true.
  enddo
endif

! And finish

close (d_unit)

lat%param%stable = .true.  ! Assume this 
inc_version = file_version

if (present(err_flag)) err_flag = err_found

if (.not. err_found .and. allocated(lat%print_str)) then
  do i = 1, size(lat%print_str)
    call out_io (s_important$, r_name, 'Message in Lattice File: ' // lat%print_str(i))
  enddo
endif

if (.not. err_found) then
  if (lat%input_taylor_order /= 0) ptc_private%taylor_order_saved = lat%input_taylor_order
  call set_ptc (1.0e12_rp, lat%param%particle)  ! Energy value used does not matter here
  call parser_init_custom_elements (lat)
endif

call create_lat_ele_nametable (lat, lat%nametable)
call ramper_slave_setup(lat, .true.)

return

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

9000  continue
if (.not. parser_call) then
  call out_io (io_err_level, r_name, 'DIGESTED FILE DOES NOT EXIST: ' // trim(full_digested_file))
endif
close (d_unit)
return

9010  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE VERSION.')
close (d_unit)
return

9020  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE FILE AND DATE.')
close (d_unit)
return

9025  continue
call out_io(io_err_level, r_name, 'ERROR READING BMAD_COM COMMON BLOCK.')
close (d_unit)
return

9030  continue
 call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE LATTICE GLOBALS.')
close (d_unit)
return

9035  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE GENERAL PARAMETER NAME LIST.')
close (d_unit)
return

9040  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE CONTROL.')
error = .true.
close (d_unit)
return

9045  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE CONTROL STACK.')
error = .true.
close (d_unit)
return

9050  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE IC.')
close (d_unit)
return

9060  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED PARTICLE_START/BEAM_INIT.')
close (d_unit)
return

9070  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE BRANCH DATA.')
close (d_unit)
return

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
contains

subroutine read_this_ele (ele, ix_ele_in, error)

type (ele_struct), target :: ele
type (photon_element_struct), pointer :: ph
type (photon_reflect_table_struct), pointer :: prt
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_struct), pointer :: ct_map
type (gen_grad_map_struct), pointer :: gg_map
type (gen_grad1_struct), pointer :: ggcoef
type (grid_field_struct), pointer :: g_field
type (ac_kicker_struct), pointer :: ac
type (wake_struct), pointer :: wake
type (converter_distribution_struct), pointer :: c_dist
type (converter_prob_pc_r_struct), pointer :: ppcr
type (converter_direction_out_struct), pointer :: c_dir
type (control_ramp1_struct), pointer ::rmp
type (wake_sr_z_long_struct), pointer :: srz

integer i, j, lb1, lb2, lb3, ub1, ub2, ub3, n_cyl, n_cart, n_gen, n_grid, ix_ele, ix_branch, ix_wall3d
integer i_min(3), i_max(3), ix_ele_in, ix_t(6), ios, k_max, ix_e, n_angle, n_energy
integer ix_r, ix_s, n_var, ix_d, ix_m, idum, n_cus, ix_convert, ix_c, nix
integer ix_sr_long, ix_sr_trans, ix_sr_z, ix_lr_mode, ix_wall3d_branch, ix_st(0:3)
integer i0, i1, j0, j1, j2, ix_ptr, lb(3), ub(3), nt, n0, n1, n2, nn(7), ne, nr, ns, nc, n_foil

logical error, is_alloc_disp, is_alloc_seg, is_alloc_h_mis, is_alloc_pix, is_alloc_ref_sigma, is_alloc_ref_pi, is_alloc_eprob
logical ac_kicker_alloc, rad_map_alloc

!

error = .true.

read (d_unit, err = 9100, end = 9100) &
        mode3, ix_r, ix_s, ix_wall3d_branch, ac_kicker_alloc, rad_map_alloc, &
        ix_convert, ix_d, ix_m, ix_t, ix_st, ix_e, ix_sr_long, ix_sr_trans, ix_sr_z, &
        ix_lr_mode, ix_wall3d, ix_c, n_cart, n_cyl, n_gen, n_grid, n_foil, n_cus, ix_convert

read (d_unit, err = 9100, end = 9100) &
        ele%name, ele%type, ele%alias, ele%component_name, ele%x, ele%y, &
        ele%a, ele%b, ele%z, ele%vec0, ele%mat6, ele%spin_q, &
        ele%c_mat, ele%gamma_c, ele%s_start, ele%s, ele%key, ele%floor, &
        ele%is_on, ele%sub_key, ele%lord_status, ele%slave_status, &
        ele%n_slave, ele%n_slave_field, ele%ix1_slave, ele%n_lord, ele%n_lord_field, ele%n_lord_ramper, &
        ele%ic1_lord, ele%ix_pointer, ele%ixx, &
        ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
        ele%spin_tracking_method, ele%symplectify, ele%mode_flip, &
        ele%multipoles_on, ele%taylor_map_includes_offsets, ele%Field_master, &
        ele%logic, ele%field_calc, ele%aperture_at, ele%spin_taylor_ref_orb_in, &
        ele%aperture_type, ele%csr_method, ele%space_charge_method, ele%orientation, &
        ele%map_ref_orb_in, ele%map_ref_orb_out, ele%time_ref_orb_in, ele%time_ref_orb_out, &
        ele%offset_moves_aperture, ele%ix_branch, ele%ref_time, ele%scale_multipoles, &
        ele%bookkeeping_state, ele%ptc_integration_type, ele%ref_species

call init_multipole_cache(ele)

! Decompress value array

read (d_unit, err = 9110, end = 9110) k_max
read (d_unit, err = 9120, end = 9120) ix_value(1:k_max), value(1:k_max)
do k = 1, k_max
  ele%value(ix_value(k)) = value(k)
enddo

! Control vars

if (ix_c /= 0) then
  allocate (ele%control)
  read (d_unit, err = 9120, end = 9120) n_var, nk, nr, nix

  if (nk > -1) then
    allocate(ele%control%x_knot(nk))
    read (d_unit, err = 9120, end = 9120) ele%control%x_knot
  endif

  if (n_var > -1) then
    allocate (ele%control%var(n_var))
    do i = 1, n_var
      read (d_unit, err = 9120, end = 9120) ele%control%var(i)
    enddo
  endif

  if (nr > -1) then
    allocate(ele%control%ramp(nr))
    do i = 1, nr
      rmp => ele%control%ramp(i)
      read (d_unit, err = 9040, end = 9040) rmp%slave_name, n, nk, rmp%attribute, rmp%is_controller
      if (n > 0) then
        allocate (rmp%stack(n))
        do j = 1, n
          read (d_unit, err = 9045, end = 9045) rmp%stack(j)
        enddo
      endif

      if (nk > 0) then
        allocate (rmp%y_knot(nk))
        read (d_unit, err = 9045, end = 9045) rmp%y_knot
      endif
    enddo
  endif

  if (nix > -1) then
    allocate(ele%control%ramper_lord(nix))
    read (d_unit, err = 9120, end = 9120) ele%control%ramper_lord%ix_ele
    read (d_unit, err = 9120, end = 9120) ele%control%ramper_lord%ix_con
    do i = 1, nix
      ele%control%ramper_lord(nix)%attrib_ptr => null()
    enddo
  endif
endif

! AC_kicker

if (ac_kicker_alloc) then
  allocate (ele%ac_kick)
  ac => ele%ac_kick
  read (d_unit, err = 9130, end = 9130) n1, n2
  if (n1 > -1) then
    allocate (ac%amp_vs_time(n1))
    do n = lbound(ac%amp_vs_time, 1), ubound(ac%amp_vs_time, 1)
      read (d_unit, err = 9130, end = 9130) ac%amp_vs_time(n)
    enddo
  endif

  if (n2 > -1) then
    allocate(ac%frequency(n2))
    do n = lbound(ac%frequency, 1), ubound(ac%frequency, 1)
      read (d_unit, err = 9130, end = 9130) ac%frequency(n)
    enddo
  endif
endif

! Foil

if (n_foil > 0) then
  allocate(ele%foil)
  allocate(ele%foil%material(n_foil))
  do n = 1, n_foil
    read (d_unit, err = 9120, end = 9120) ele%foil%material(n)
  enddo
endif

! Converter

if (ix_convert == 1) then
  allocate (ele%converter)
  read (d_unit, err = 9120, end = 9120) ele%converter%species_out, ele%converter%material_type, ns
  allocate (ele%converter%dist(ns))
  do n = 1, size(ele%converter%dist)
    c_dist => ele%converter%dist(n)
    read (d_unit, err = 9120, end = 9120) c_dist%thickness, ns
    allocate (c_dist%sub_dist(ns))
    do j = 1, size(c_dist%sub_dist)
      read (d_unit, err = 9120, end = 9120) c_dist%sub_dist(j)%pc_in
      ppcr => c_dist%sub_dist(j)%prob_pc_r
      read (d_unit, err = 9120, end = 9120) ppcr%integrated_prob, ne, nr
      allocate (ppcr%pc_out(ne), ppcr%r(nr), ppcr%prob(ne,nr), ppcr%spin_z(ne,nr))
      read (d_unit, err = 9120, end = 9120) ppcr%pc_out
      read (d_unit, err = 9120, end = 9120) ppcr%r
      read (d_unit, err = 9120, end = 9120) ppcr%prob
      read (d_unit, err = 9120, end = 9120) ppcr%spin_z
      c_dir => c_dist%sub_dist(j)%dir_out
      read (d_unit, err = 9120, end = 9120) nn
      allocate (c_dir%beta%fit_1D_r(nn(1)), c_dir%alpha_x%fit_1D_r(nn(2)), c_dir%alpha_y%fit_1D_r(nn(3)), &
                c_dir%c_x%fit_1D_r(nn(4)), c_dir%dxds_min%fit_1D_r(nn(5)), c_dir%dxds_max%fit_1D_r(nn(6)), &
                c_dir%dyds_max%fit_1D_r(nn(7)))
      read (d_unit, err = 9120, end = 9120) c_dir%beta%fit_1d_r, c_dir%beta%fit_2d_pc, c_dir%beta%fit_2d_r, c_dir%beta%c0
      read (d_unit, err = 9120, end = 9120) c_dir%alpha_x%fit_1d_r, c_dir%alpha_x%fit_2d_pc, c_dir%alpha_x%fit_2d_r, c_dir%alpha_x%c0
      read (d_unit, err = 9120, end = 9120) c_dir%alpha_y%fit_1d_r, c_dir%alpha_y%fit_2d_pc, c_dir%alpha_y%fit_2d_r, c_dir%alpha_y%c0
      read (d_unit, err = 9120, end = 9120) c_dir%c_x%fit_1d_r, c_dir%c_x%fit_2d_pc, c_dir%c_x%fit_2d_r, c_dir%c_x%c0
      read (d_unit, err = 9120, end = 9120) c_dir%dxds_min%fit_1d_r, c_dir%dxds_min%fit_2d_pc, c_dir%dxds_min%fit_2d_r, c_dir%dxds_min%c0
      read (d_unit, err = 9120, end = 9120) c_dir%dxds_max%fit_1d_r, c_dir%dxds_max%fit_2d_pc, c_dir%dxds_max%fit_2d_r, c_dir%dxds_max%c0
      read (d_unit, err = 9120, end = 9120) c_dir%dyds_max%fit_1d_r, c_dir%dyds_max%fit_2d_pc, c_dir%dyds_max%fit_2d_r, c_dir%dyds_max%c0
    enddo
  enddo
endif

! Cartesian_map

if (n_cart > 0) then
  allocate (ele%cartesian_map(n_cart))

  do i = 1, n_cart
    ct_map => ele%cartesian_map(i)
    read (d_unit, err = 9120, end = 9120) ct_map%field_scale, ct_map%master_parameter, ct_map%ele_anchor_pt, ct_map%field_type, ct_map%r0
    read (d_unit, err = 9120, end = 9120) ix_ele, ix_branch, ix_ptr, n

    if (ix_ele > 0) then
      ele%cartesian_map(i)%ptr => lat%branch(ix_branch)%ele(ix_ele)%cartesian_map(ix_ptr)%ptr
      ele%cartesian_map(i)%ptr%n_link = ele%cartesian_map(i)%ptr%n_link + 1
    else
      allocate (ele%cartesian_map(i)%ptr)
      read (d_unit, err = 9120, end = 9120) ct_map%ptr%file
      allocate (ct_map%ptr%term(n))
      do j = 1, n
        read (d_unit, err = 9120, end = 9120) ct_map%ptr%term(j)
      enddo
    endif
  enddo
endif

! Cylindrical map

if (n_cyl > 0) then
  allocate (ele%cylindrical_map(n_cyl))

  do i = 1, n_cyl
    cl_map => ele%cylindrical_map(i)
    read (d_unit, err = 9120, end = 9120) cl_map%field_scale, cl_map%master_parameter, cl_map%harmonic, &
                cl_map%phi0_fieldmap, cl_map%theta0_azimuth, cl_map%ele_anchor_pt, cl_map%m, cl_map%dz, cl_map%r0
    read (d_unit, err = 9120, end = 9120) ix_ele, ix_branch, ix_ptr, n

    if (ix_ele > 0) then
      ele%cylindrical_map(i)%ptr => lat%branch(ix_branch)%ele(ix_ele)%cylindrical_map(ix_ptr)%ptr
      ele%cylindrical_map(i)%ptr%n_link = ele%cylindrical_map(i)%ptr%n_link + 1
    else
      allocate (ele%cylindrical_map(i)%ptr)
      read (d_unit, err = 9120, end = 9120) cl_map%ptr%file
      allocate (cl_map%ptr%term(n))
      do j = 1, n
        read (d_unit, err = 9120, end = 9120) cl_map%ptr%term(j)
      enddo
    endif
  enddo
endif

! Gen_grad_field

if (n_gen > 0) then
  allocate (ele%gen_grad_map(n_gen))

  do i = 1, n_gen
    gg_map => ele%gen_grad_map(i)

    read (d_unit, err = 9120, end = 9120) gg_map%field_scale, gg_map%master_parameter, gg_map%curved_ref_frame, &
                            gg_map%ele_anchor_pt, gg_map%field_type, gg_map%dz, gg_map%r0, ns, gg_map%iz0, gg_map%iz1
    allocate (gg_map%gg(ns))
    n0 = gg_map%iz0;  n1 = gg_map%iz1

    do j = 1, size(gg_map%gg)
      ggcoef => gg_map%gg(j)
      read (d_unit, err = 9120, end = 9120) ggcoef%m, ggcoef%sincos, ggcoef%n_deriv_max, lb2
      allocate (ggcoef%deriv(n0:n1, 0:lb2))
      do k = n0, n1
        read (d_unit, err = 9120, end = 9120) ggcoef%deriv(k,:)
      enddo
    enddo
  enddo
endif

! Grid_field

if (n_grid > 0) then
  allocate (ele%grid_field(n_grid))

  do i = 1, n_grid
    g_field => ele%grid_field(i)
    read (d_unit, err = 9120, end = 9120) g_field%field_scale, g_field%master_parameter, &
                g_field%ele_anchor_pt, g_field%phi0_fieldmap, g_field%dr, &
                g_field%r0, g_field%harmonic, g_field%geometry, &
                g_field%curved_ref_frame, g_field%field_type, g_field%interpolation_order
    read (d_unit, err = 9120, end = 9120) ix_ele, ix_branch, ix_ptr, lb, ub

    if (ix_ele > 0) then
      ele%grid_field(i)%ptr => lat%branch(ix_branch)%ele(ix_ele)%grid_field(ix_ptr)%ptr
      ele%grid_field(i)%ptr%n_link = ele%grid_field(i)%ptr%n_link + 1
    else
      allocate (ele%grid_field(i)%ptr)
      read (d_unit, err = 9120, end = 9120) g_field%ptr%file
      allocate (g_field%ptr%pt(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
      do j = lb(3), ub(3)
        read (d_unit, err = 9120, end = 9120) g_field%ptr%pt(:, :, j)
      enddo
    endif
  enddo
endif

! Mode3

if (mode3) then
  allocate(ele%mode3)
  read (d_unit, err = 9150, end = 9150) ele%mode3
endif

if (ix_r /= 0) then
  read (d_unit, err = 9350, end = 9350) i_min, i_max
  allocate (ele%r(i_min(1):i_max(1), i_min(2):i_max(2), i_min(3):i_max(3)))
  do i = i_min(3), i_max(3)
    read (d_unit, err = 9400, end = 9400) ele%r(:,:,i)
  enddo
endif

if (n_cus /= 0) then
  allocate (ele%custom(n_cus))
  read (d_unit, err = 9410, end = 9410) ele%custom(:)
endif

if (ix_s /= 0) then
  allocate (ele%photon)
  ph => ele%photon
  read (d_unit, err = 9360, end = 9360) ph%target, ph%material, ph%curvature, &
    ph%displacement%active, ph%displacement%dr, ph%displacement%r0, is_alloc_disp, &
    ph%h_misalign%active, ph%h_misalign%dr, ph%h_misalign%r0, is_alloc_h_mis, &
    ph%segmented%active, ph%segmented%dr, ph%segmented%r0, is_alloc_seg, &
    ph%pixel%dr, ph%pixel%r0, is_alloc_pix, &
    is_alloc_eprob, is_alloc_ref_sigma, is_alloc_ref_pi
         

  if (is_alloc_disp) then
    read (d_unit, err = 9361, end = 9361) i0, j0, i1, j1
    allocate(ph%displacement%pt(i0:i1, j0:j1))
    do i = lbound(ph%displacement%pt, 1), ubound(ph%displacement%pt, 1)
    do j = lbound(ph%displacement%pt, 2), ubound(ph%displacement%pt, 2)
      read (d_unit, err = 9362, end = 9362) ph%displacement%pt(i,j)
    enddo
    enddo
  endif

  if (is_alloc_seg) then
    read (d_unit, err = 9361, end = 9361) i0, j0, i1, j1
    allocate(ph%segmented%pt(i0:i1, j0:j1))
    do i = lbound(ph%segmented%pt, 1), ubound(ph%segmented%pt, 1)
    do j = lbound(ph%segmented%pt, 2), ubound(ph%segmented%pt, 2)
      read (d_unit, err = 9362, end = 9362) ph%segmented%pt(i,j)
    enddo
    enddo
  endif

  if (is_alloc_h_mis) then
    read (d_unit, err = 9361, end = 9361) i0, j0, i1, j1
    allocate(ph%h_misalign%pt(i0:i1, j0:j1))
    do i = lbound(ph%h_misalign%pt, 1), ubound(ph%h_misalign%pt, 1)
    do j = lbound(ph%h_misalign%pt, 2), ubound(ph%h_misalign%pt, 2)
      read (d_unit, err = 9362, end = 9362) ph%h_misalign%pt(i,j)
    enddo
    enddo
  endif

  if (is_alloc_pix) then
    read (d_unit, err = 9361, end = 9361) i0, j0, i1, j1
    allocate(ph%pixel%pt(i0:i1, j0:j1))
    ! Note: At startup detectors do not have any grid data that needs saving
  endif

  if (is_alloc_ref_sigma) then
    prt => ph%reflectivity_table_sigma
    read (d_unit, err = 9370, end = 9370) n_energy, n_angle
    allocate(prt%angle(n_angle), prt%energy(n_energy))
    allocate(prt%p_reflect(n_angle,n_energy), prt%bragg_angle(n_energy))
    read (d_unit, err = 9370, end = 9370) prt%angle
    read (d_unit, err = 9370, end = 9370) prt%energy
    read (d_unit, err = 9370, end = 9370) prt%bragg_angle
    do j = 1, n_energy
      read (d_unit, err = 9370, end = 9370) prt%p_reflect(:,j)
    enddo
  endif

  if (is_alloc_ref_pi) then
    prt => ph%reflectivity_table_pi
    read (d_unit, err = 9380, end = 9380) n_energy, n_angle
    allocate(prt%angle(n_angle), prt%energy(n_energy), prt%p_reflect_scratch(n_angle))
    allocate(prt%p_reflect(n_angle,n_energy), prt%int1(n_energy))
    read (d_unit, err = 9380, end = 9380) prt%angle
    read (d_unit, err = 9380, end = 9380) prt%energy
    read (d_unit, err = 9380, end = 9380) prt%bragg_angle
    do j = 1, n_energy
      read (d_unit, err = 9380, end = 9380) prt%p_reflect(:,j)
    enddo
  endif

  if (is_alloc_ref_sigma .and. is_alloc_ref_pi) then
    ph%reflectivity_table_type = polarized$
  elseif (is_alloc_ref_sigma .or. is_alloc_ref_pi) then
    ph%reflectivity_table_type = unpolarized$
  endif

  if (is_alloc_eprob) then
    read (d_unit, err = 9390, end = 9390) n
    allocate (ph%init_energy_prob(n), ph%integrated_init_energy_prob(n))
    read (d_unit, err = 9390, end = 9390) ph%init_energy_prob
    read (d_unit, err = 9390, end = 9390) ph%integrated_init_energy_prob
  endif
endif

if (ix_d /= 0) then
  allocate (ele%descrip)
  read (d_unit, err = 9500, end = 9500) ele%descrip
endif

if (ix_m /= 0) then
  call multipole_init (ele, magnetic$)
  read (d_unit, err = 9600, end = 9600) ele%a_pole, ele%b_pole
endif
  
if (ix_e /= 0) then
  call multipole_init (ele, electric$)
  read (d_unit, err = 9600, end = 9600) ele%a_pole_elec, ele%b_pole_elec
endif

do j = 1, size(ele%taylor)
  if (ix_t(j) == -1) cycle
  read (d_unit, err = 9650, end = 9650) ele%taylor(j)%ref
  allocate (ele%taylor(j)%term(ix_t(j)))
  do k = 1, ix_t(j)
    read (d_unit, err = 9700, end = 9700) ele%taylor(j)%term(k)
  enddo
enddo

do i = 0, 3
  if (ix_st(i) == -1) cycle
  read (d_unit, err = 9650, end = 9650) ele%spin_taylor(i)%ref
  allocate (ele%spin_taylor(i)%term(ix_st(i)))
  do k = 1, ix_st(i)
    read (d_unit, err = 9700, end = 9700) ele%spin_taylor(i)%term(k)
  enddo
enddo

! If ix_lr_mode is negative then it is a pointer to a previously read wake. 
! See write_digested_bmad_file.

if (ix_sr_long /= 0 .or. ix_sr_trans /= 0 .or. ix_sr_z /= 0 .or. ix_lr_mode /= 0) then
  if (ix_lr_mode < 0) then
    call transfer_wake (ele%branch%ele(abs(ix_lr_mode))%wake, ele%wake)

  else
    call init_wake (ele%wake, ix_sr_long, ix_sr_trans, ix_sr_z, ix_lr_mode)
    wake => ele%wake
    read (d_unit, err = 9800, end = 9800) wake%sr%z_ref_long, wake%sr%z_ref_trans, wake%sr%z_max, wake%sr%scale_with_length, wake%sr%amp_scale, wake%sr%z_scale

    do i = 1, size(wake%sr%long)
      read (d_unit, err = 9800, end = 9800) wake%sr%long(i)
    enddo

    do i = 1, size(wake%sr%trans)
      read (d_unit, err = 9800, end = 9800) wake%sr%trans(i)
    enddo

    srz => wake%sr%z_long
    read (d_unit, err = 9800, end = 9800) srz%smoothing_sigma, srz%position_dependence, srz%dz, srz%z0, srz%time_based
    do i = 1, size(srz%w)
      read (d_unit, err = 9800, end = 9800) srz%w(i), srz%fw(i)
    enddo

    read (d_unit, err = 9800, end = 9800) wake%lr%t_ref, wake%lr%freq_spread, wake%lr%self_wake_on, wake%lr%amp_scale, wake%lr%time_scale

    do i = 1, size(wake%lr%mode)
      read (d_unit, err = 9800, end = 9800) wake%lr%mode(i)
    enddo
  endif
endif

if (ix_wall3d > 0) then
  call read_this_wall3d (ele%wall3d, error)
  if (error) return
elseif (ix_wall3d < 0) then
  read (d_unit, err = 9900, end = 9900) idum
  ele%wall3d => lat%branch(ix_wall3d_branch)%ele(abs(ix_wall3d))%wall3d
  if (.not. associated(ele%wall3d)) then
    call out_io(io_err_level, r_name, 'ERROR IN WALL3D INIT.')
    close (d_unit)
    return
  endif
  ele%wall3d%n_link = ele%wall3d%n_link + 1
else
  read (d_unit, err = 9900, end = 9900) idum
endif

!

if (rad_map_alloc) then
  allocate (ele%rad_map)
  read (d_unit, err = 9900, end = 9900) ele%rad_map%rm0
  read (d_unit, err = 9900, end = 9900) ele%rad_map%rm1, ele%rad_map%stale
endif

!

ele%old_value = ele%value

error = .false.

return

!--------------------------------------------------------------

9040  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE CONTROL.')
error = .true.
close (d_unit)
return

9045  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE CONTROL STACK.')
error = .true.
close (d_unit)
return

9100  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9110  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING K_MAX OF ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9120  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING VALUES OF ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9130  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING AC_KICKER VALUES OF ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9140  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING EM_FIELD COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9150  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING MODE3 COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9350  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %R ARRAY SIZE: ' // ele%name)
close (d_unit)
return

9360  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9361  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%GRID BOUNDS FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9362  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%SURFACE%GRID FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9370  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%REFLECTIVITY_TABLE_SIGMA FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9380  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%REFLECTIVITY_TABLE_PI FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9390  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%INTEGRATED_INIT_ENERGY_PROB FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9400  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING R COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9410  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING CUSTOM COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9500  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING DESCRIP COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9600  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING AN,BN COMPONENTS FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9650  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %TAYLOR(:)%REF FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9700  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING TAYLOR COMPONENTS FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9800  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9900  continue
call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING IDUM1 FOR ELEMENT: ' // ele%name)
close (d_unit)
return

end subroutine read_this_ele 

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! contains

subroutine read_this_wall3d (wall3d, error)

type (wall3d_struct), pointer :: wall3d(:)

integer i, j, k, n_wall, n_wall_section, ios
logical error

!

error = .true.

read (d_unit, iostat = ios) n_wall
if (n_wall > 0) allocate(wall3d(n_wall))

do i = 1, n_wall
  read (d_unit, iostat = ios) n_wall_section, wall3d(i)%type, &
          wall3d(i)%ele_anchor_pt, wall3d(i)%superimpose, &
          wall3d(i)%thickness, wall3d(i)%clear_material, wall3d(i)%opaque_material

  if (n_wall_section == 0) then
    error = .false.
    return
  endif

  if (ios /= 0) then
     call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', &
                                     'ERROR READING WALL3D N_WALL_SECTION NUMBER')
    close (d_unit)
    return
  endif

  call re_allocate (wall3d(i)%section, n_wall_section)

  do j = 1, n_wall_section
    call read_this_wall3d_section (wall3d(i)%section(j))
  enddo
enddo

error = .false.

end subroutine read_this_wall3d

!-----------------------------------------------
! contains

subroutine read_this_wall3d_section (sec)

type (wall3d_section_struct), target :: sec
type (wall3d_vertex_struct), pointer :: v
integer nv

!

read (d_unit, iostat = ios) sec%name, sec%material, sec%type, sec%n_vertex_input, &
                   sec%ix_ele, sec%ix_branch, sec%patch_in_region, &
                   sec%thickness, sec%s, sec%r0, sec%dx0_ds, sec%dy0_ds, &
                   sec%x0_coef, sec%y0_coef, sec%dr_ds, sec%p1_coef, sec%p2_coef, nv
if (ios /= 0) then
  call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', 'ERROR READING WALL3D SECTION')
  close (d_unit)
  return
endif

allocate(sec%v(nv))
do k = 1, nv
  read (d_unit, iostat = ios) sec%v(k)
  if (ios /= 0) then
    call out_io(io_err_level, r_name, 'ERROR READING DIGESTED FILE.', 'ERROR READING WALL3D VERTEX')
    close (d_unit)
    return
  endif
enddo

end subroutine read_this_wall3d_section

end subroutine read_digested_bmad_file
