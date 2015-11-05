!+
! Subroutine read_digested_bmad_file (digested_file, lat, inc_version, err_flag)
!
! Subroutine to read in a digested file. The subroutine will check that
! the version of the digested file is up to date and that the digested file
! is current with respect to the original BMAD files that were used. [See
! write_digested_bmad_file.]
!
! Note: This subroutine also reads in the common structures for bmad_parser2
!
! Modules Needed:
!   use bmad
!
! Input:
!   digested_file -- Character(*): Name of the digested file.
!
! Output:
!   lat         -- lat_struct: Output lattice structure
!   inc_version -- Integer: bmad_inc_version number stored in the lattice file.
!                   If the file is current this number should be the same 
!                   as the global parameter bmad_inc_version$. 
!                   Set to -1 if there is a read error.
!   err_flag    -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine read_digested_bmad_file (digested_file, lat, inc_version, err_flag)

use ptc_interface_mod, dummy => read_digested_bmad_file

implicit none

type (lat_struct), target, intent(inout) :: lat
type (branch_struct), pointer :: branch
type (extra_parsing_info_struct) :: extra
type (ptc_parameter_struct) ptc_param
type (control_struct), pointer :: c
type (bmad_common_struct) bmad_com_read

real(rp) value(num_ele_attrib$)

integer inc_version, d_unit, n_files, file_version, i, j, k, ix, ix_value(num_ele_attrib$)
integer stat_b(13), stat_b2, stat_b8, stat_b10, n_branch, n, control_type, coupler_at, idum1
integer ierr, stat, ios, ios2, n_wall_section, garbage, idum2, j1, j2

character(*) digested_file
character(200) fname_read, fname_versionless, fname_full
character(200) input_file_name, full_digested_file, digested_prefix_in, digested_prefix_out
character(200), allocatable :: file_names(:)
character(60) p_name
character(25) :: r_name = 'read_digested_bmad_file'

logical, optional :: err_flag
logical is_ok, l_dummy
logical found_it, can_read_this_old_version, mode3, error, is_match, err, err_found

! init all elements in lat

if (present(err_flag)) err_flag = .true.
err_found = .false.

call init_lat (lat)

! Read the digested file.
! Some old versions can be read even though they are not the current version.

d_unit = lunget()
inc_version = -1
lat%n_ele_track = 0

call fullfilename (digested_file, full_digested_file)
inquire (file = full_digested_file, name = full_digested_file)
call simplify_path (full_digested_file, full_digested_file)
open (unit = d_unit, file = full_digested_file, status = 'old',  &
                     form = 'unformatted', action = 'READ', err = 9000)

read (d_unit, err = 9010) n_files, file_version
allocate (file_names(n_files))

! Version is old but recent enough to read.

can_read_this_old_version = .false.

if (file_version < bmad_inc_version$) then
  call out_io (s_info$, r_name, ['DIGESTED FILE VERSION OUT OF DATE \i0\ > \i0\ ' ],  &
                                i_array = [bmad_inc_version$, file_version ])
  if (can_read_this_old_version) then 
    err_found = .true.
  else
    close (d_unit)
    return
  endif
endif

if (file_version > bmad_inc_version$) then
  call out_io (s_info$, r_name, &
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

  read (d_unit, err = 9020, end = 9020) fname_read, stat_b2, stat_b8, stat_b10, idum1, idum2
  file_names(i) = fname_read

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
      if (.not. err_found) call out_io(s_info$, r_name, ' NOTE: MOVED DIGESTED FILE.')
      err_found = .true.
    endif
    cycle
  endif

  if (fname_read(1:7) == '!PRINT:') cycle  ! Only print at end if no errors.

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
    if (.not. err_found) call out_io(s_info$, r_name, 'NOTE: DIGESTED FILE OUT OF DATE.')
    err_found = .true.
  endif

enddo

! we read (and write) the lat in pieces since it is
! too big to write in one piece

read (d_unit, err = 9030)  &   
        lat%use_name, lat%lattice, lat%input_file_name, lat%title, &
        lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track, &
        lat%n_ele_track, lat%n_ele_max, lat%lord_state, lat%n_control_max, lat%n_ic_max, &
        lat%input_taylor_order, lat%absolute_time_tracking, l_dummy, lat%pre_tracker, lat%photon_type
read (d_unit, err = 9070) n_branch

! custom attribute names

read (d_unit, err = 9035) n
if (n > 0) then
  allocate (lat%attribute_alias(n))
  do i = 1, n
    read (d_unit, err = 9035) lat%attribute_alias(i)
  enddo
endif

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

do i = 1, n_branch
  branch => lat%branch(i)
  branch%ix_branch = i
  read (d_unit, err = 9070) branch%param
  read (d_unit, err = 9070) branch%name, branch%ix_from_branch, &
                 branch%ix_from_ele, branch%n_ele_track, branch%n_ele_max, idum1

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
  c => lat%control(i)
  read (d_unit, err = 9040) n, c%ix_lord, c%slave, c%ix_attrib
  if (n > 0) then
    allocate (c%stack(n))
    do j = 1, n
      read (d_unit, err = 9045) c%stack(j)
    enddo
  endif
enddo

do i = 1, lat%n_ic_max
  read (d_unit, err = 9050) lat%ic(i)
enddo

read (d_unit, err = 9060) lat%beam_start

! Read PTC info

read (d_unit, iostat = ios) ptc_param
if (ios /= 0) then
  call out_io(s_error$, r_name, 'ERROR READING PTC PARAMETERS.')
  close (d_unit)
  return
endif

if (ios == 0) then
  call set_ptc (exact_modeling = ptc_param%exact_model, exact_misalign = ptc_param%exact_misalign)
endif

! Read extra state info.

read (d_unit, iostat = ios) found_it
if (found_it) then
  read (d_unit, iostat = ios) extra
  read (d_unit, iostat = ios2) bmad_com_read
  if (ios /= 0 .or. ios2 /= 0) then
    call out_io (s_error$, r_name, 'ERROR READING BMAD COMMON PARAMETERS')
    close (d_unit)
    return
  endif
  if (extra%max_aperture_limit_set)             bmad_com%max_aperture_limit              = bmad_com_read%max_aperture_limit
  if (extra%default_ds_step_set)                bmad_com%default_ds_step                 = bmad_com_read%default_ds_step
  if (extra%significant_length_set)             bmad_com%significant_length              = bmad_com_read%significant_length
  if (extra%rel_tol_tracking_set)               bmad_com%rel_tol_tracking                = bmad_com_read%rel_tol_tracking
  if (extra%abs_tol_tracking_set)               bmad_com%abs_tol_tracking                = bmad_com_read%abs_tol_tracking
  if (extra%rel_tol_adaptive_tracking_set)      bmad_com%rel_tol_adaptive_tracking       = bmad_com_read%rel_tol_adaptive_tracking
  if (extra%abs_tol_adaptive_tracking_set)      bmad_com%abs_tol_adaptive_tracking       = bmad_com_read%abs_tol_adaptive_tracking
  if (extra%init_ds_adaptive_tracking_set)      bmad_com%init_ds_adaptive_tracking       = bmad_com_read%init_ds_adaptive_tracking
  if (extra%min_ds_adaptive_tracking_set)       bmad_com%min_ds_adaptive_tracking        = bmad_com_read%min_ds_adaptive_tracking
  if (extra%fatal_ds_adaptive_tracking_set)     bmad_com%fatal_ds_adaptive_tracking      = bmad_com_read%fatal_ds_adaptive_tracking
  if (extra%electric_dipole_moment_set)         bmad_com%electric_dipole_moment          = bmad_com_read%electric_dipole_moment
  if (extra%taylor_order_set)                   bmad_com%taylor_order                    = bmad_com_read%taylor_order
  if (extra%default_integ_order_set)            bmad_com%default_integ_order             = bmad_com_read%default_integ_order
  if (extra%ptc_max_fringe_order_set)           bmad_com%ptc_max_fringe_order            = bmad_com_read%ptc_max_fringe_order
  if (extra%use_hard_edge_drifts_set)           bmad_com%use_hard_edge_drifts            = bmad_com_read%use_hard_edge_drifts
  if (extra%sr_wakes_on_set)                    bmad_com%sr_wakes_on                     = bmad_com_read%sr_wakes_on
  if (extra%lr_wakes_on_set)                    bmad_com%lr_wakes_on                     = bmad_com_read%lr_wakes_on
  if (extra%mat6_track_symmetric_set)           bmad_com%mat6_track_symmetric            = bmad_com_read%mat6_track_symmetric
  if (extra%auto_bookkeeper_set)                bmad_com%auto_bookkeeper                 = bmad_com_read%auto_bookkeeper
  if (extra%space_charge_on_set)                bmad_com%space_charge_on                 = bmad_com_read%space_charge_on
  if (extra%coherent_synch_rad_on_set)          bmad_com%coherent_synch_rad_on           = bmad_com_read%coherent_synch_rad_on
  if (extra%spin_tracking_on_set)               bmad_com%spin_tracking_on                = bmad_com_read%spin_tracking_on
  if (extra%radiation_damping_on_set)           bmad_com%radiation_damping_on            = bmad_com_read%radiation_damping_on
  if (extra%radiation_fluctuations_on_set)      bmad_com%radiation_fluctuations_on       = bmad_com_read%radiation_fluctuations_on
  if (extra%conserve_taylor_maps_set)           bmad_com%conserve_taylor_maps            = bmad_com_read%conserve_taylor_maps
  if (extra%absolute_time_tracking_default_set) bmad_com%absolute_time_tracking_default  = bmad_com_read%absolute_time_tracking_default
  if (extra%convert_to_kinetic_momentum_set)    bmad_com%convert_to_kinetic_momentum     = bmad_com_read%convert_to_kinetic_momentum
  if (extra%aperture_limit_on_set)              bmad_com%aperture_limit_on               = bmad_com_read%aperture_limit_on

endif

! Setup any attribute aliases in the global attribute name table.
! This is done last in read_digested bmad_file so as to not to pollute the name table if 
! there is an error.

if (.not. err_found .and. allocated(lat%attribute_alias)) then
  do i = 1, size(lat%attribute_alias)
    p_name = lat%attribute_alias(i)
    ix = index(lat%attribute_alias(i), '=')
    call set_attribute_alias(p_name(1:ix-1), p_name(ix+1:), err)
    if (err) err_found = .true.
  enddo
endif

! And finish

close (d_unit)

lat%param%stable = .true.  ! Assume this 
inc_version = file_version

if (present(err_flag)) err_flag = err_found

if (.not. err_found) then
  do i = 1, size(file_names)
    fname_read = file_names(i)
    if (fname_read(1:7) /= '!PRINT:') cycle
    call out_io (s_dwarn$, r_name, 'Print Message in Lattice File: ' // fname_read(8:))
  enddo
endif

return

!------------------

9000  continue
!! call out_io (s_warn$, r_name, 'Digested file does not exist: ' // trim(full_digested_file))
close (d_unit)
return

!--------------------------------------------------------------

9010  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE VERSION.')
close (d_unit)
return

!--------------------------------------------------------------

9020  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE FILE AND DATE.')
close (d_unit)
return

!--------------------------------------------------------------

9025  continue
call out_io(s_error$, r_name, 'ERROR READING BMAD_COM COMMON BLOCK.')
close (d_unit)
return

!--------------------------------------------------------------

9030  continue
 call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE LATTICE GLOBALS.')
close (d_unit)
return

!--------------------------------------------------------------

9035  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE GENERAL PARAMETER NAME LIST.')
close (d_unit)
return

!--------------------------------------------------------------

9040  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE CONTROL.')
close (d_unit)
return

!--------------------------------------------------------------

9045  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE CONTROL STACK.')
close (d_unit)
return

!--------------------------------------------------------------

9050  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE IC.')
close (d_unit)
return

!--------------------------------------------------------------

9060  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED BEAM_INIT.')
close (d_unit)
return

!--------------------------------------------------------------

9070  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE BRANCH DATA.')
close (d_unit)
return

!--------------------------------------------------------------

9080  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED RANDOM NUMBER STATE.')
close (d_unit)
return

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
contains

subroutine read_this_ele (ele, ix_ele_in, error)

type (ele_struct), target :: ele
type (em_field_mode_struct), pointer :: mode
type (photon_surface_struct), pointer :: surf
type (surface_grid_pt_struct), pointer :: s_pt

real(rp) rdum

integer i, j, lb1, lb2, lb3, ub1, ub2, ub3, nf, ng, ix_ele, ix_branch, ix_wall3d
integer n_em_field_mode, i_min(3), i_max(3), ix_ele_in, ix_t(6), ios, k_max, ix_e
integer ix_wig, ix_r, ix_s, ix_wig_branch, idum1, idum2, idum3, n_var, ix_d, ix_m
integer ix_sr_long, ix_sr_trans, ix_lr, ix_wall3d_branch, ix_st(3,3)
integer i0, i1, j0, j1, j2

logical error, is_alloc_pt

!

error = .true.

read (d_unit, err = 9100, end = 9100) &
        mode3, ix_wig, ix_wig_branch, ix_r, ix_s, ix_wall3d_branch, &
        idum2, idum3, ix_d, ix_m, ix_t, ix_st, ix_e, ix_sr_long, ix_sr_trans, &
        ix_lr, ix_wall3d, n_em_field_mode, n_var
read (d_unit, err = 9100, end = 9100) &
        ele%name, ele%type, ele%alias, ele%component_name, ele%x, ele%y, &
        ele%a, ele%b, ele%z, rdum, ele%vec0, ele%mat6, &
        ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
        ele%is_on, ele%sub_key, ele%lord_status, ele%slave_status, &
        ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
        ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
        ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
        ele%spin_tracking_method, ele%symplectify, ele%mode_flip, &
        ele%multipoles_on, ele%taylor_map_includes_offsets, ele%Field_master, &
        ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
        ele%aperture_type, ele%csr_calc_on, ele%orientation, &
        ele%map_ref_orb_in, ele%map_ref_orb_out, ele%offset_moves_aperture, &
        ele%ix_branch, ele%ref_time, ele%scale_multipoles, idum1, &
        idum2, ele%bookkeeping_state, ele%ptc_integration_type

! Decompress value array

read (d_unit, err = 9110) k_max
read (d_unit, err = 9120) ix_value(1:k_max), value(1:k_max)
do k = 1, k_max
  ele%value(ix_value(k)) = value(k)
enddo

! Control vars

if (n_var /= 0) then
  allocate (ele%control_var(n_var))
  do i = 1, n_var
    read (d_unit) ele%control_var(i)
  enddo
endif

! EM field def

call init_em_field (ele%em_field, n_em_field_mode)

if (n_em_field_mode > 0) then

  read (d_unit) ele%em_field%mode_to_autoscale

  do i = 1, n_em_field_mode
    mode => ele%em_field%mode(i)
    read (d_unit, err = 9140) nf, ng, ix_ele, ix_branch, mode%harmonic, mode%f_damp, mode%phi0_ref, &
                         mode%stored_energy, mode%m, mode%phi0_azimuth, mode%field_scale, mode%master_scale

    if (ix_ele /= 0) then
      mode%map  => lat%branch(ix_branch)%ele(ix_ele)%em_field%mode(i)%map
      if (associated(mode%map)) mode%map%n_link = mode%map%n_link + 1
      mode%grid => lat%branch(ix_branch)%ele(ix_ele)%em_field%mode(i)%grid
      if (associated(mode%grid)) mode%grid%n_link = mode%grid%n_link + 1
      cycle
    endif

    if (nf > 0) then
      allocate (mode%map)
      allocate (mode%map%term(nf))
      read (d_unit, err = 9140) mode%map%file, mode%map%dz, mode%map%ele_anchor_pt
      read (d_unit, err = 9140) mode%map%term
    endif

    if (ng > 0) then
      allocate (mode%grid)
      read (d_unit, err = 9140) lb1, ub1, lb2, ub2, lb3, ub3, mode%grid%type, mode%grid%file, &
                                mode%grid%dr, mode%grid%r0, mode%grid%ele_anchor_pt
      allocate (mode%grid%pt(lb1:ub1, lb2:ub2, lb3:ub3))
      do j = lb3, ub3
        read (d_unit, err = 9140) mode%grid%pt(:,:,j)
      enddo
    endif

  enddo
endif

! Mode3

if (mode3) then
  allocate(ele%mode3)
  read (d_unit, err = 9150) ele%mode3
endif

if (ix_wig > 0) then
  allocate (ele%wig)
  ele%wig%n_link = 1
  allocate (ele%wig%term(ix_wig))
  do j = 1, ix_wig
    read (d_unit, err = 9200) ele%wig%term(j)
  enddo
elseif (ix_wig < 0) then  ! Points to another wiggler with same field
  ele%wig => lat%branch(ix_wig_branch)%ele(abs(ix_wig))%wig
  ele%wig%n_link = ele%wig%n_link + 1
endif

if (ix_r /= 0) then
  read (d_unit, err = 9350) i_min, i_max
  allocate (ele%r(i_min(1):i_max(1), i_min(2):i_max(2), i_min(3):i_max(3)))
  do i = i_min(3), i_max(3)
    read (d_unit, err = 9400) ele%r(:,:,i)
  enddo
endif

if (ix_s /= 0) then
  allocate (ele%photon)
  surf => ele%photon%surface
  read (d_unit, err = 9360) ele%photon%target, ele%photon%material, &
         surf%curvature_xy, surf%has_curvature, &
         surf%grid%type, surf%grid%dr, surf%grid%r0, surf%segment, is_alloc_pt
  if (is_alloc_pt) then
    read (d_unit, err = 9361) i0, j0, i1, j1
    allocate(surf%grid%pt(i0:i1, j0:j1))
    ! Detectors do not have any grid data that needs saving
    if (ele%key /= detector$) then
      do i = lbound(surf%grid%pt, 1), ubound(surf%grid%pt, 1)
      do j = lbound(surf%grid%pt, 2), ubound(surf%grid%pt, 2)
        read (d_unit, err = 9362) surf%grid%pt(i,j)
      enddo
      enddo
    endif
  endif

endif

if (ix_d /= 0) then
  allocate (ele%descrip)
  read (d_unit, err = 9500) ele%descrip
endif

if (ix_m /= 0) then
  call multipole_init (ele)
  read (d_unit, err = 9600) ele%a_pole, ele%b_pole
endif
  
if (ix_e /= 0) then
  call elec_multipole_init (ele)
  read (d_unit, err = 9600) ele%a_pole_elec, ele%b_pole_elec
endif

do j = 1, size(ele%taylor)
  if (ix_t(j) == -1) cycle
  read (d_unit, err = 9650) ele%taylor(j)%ref
  allocate (ele%taylor(j)%term(ix_t(j)))
  do k = 1, ix_t(j)
    read (d_unit, err = 9700) ele%taylor(j)%term(k)
  enddo
enddo

do i = 1, 3; do j = 1, 3
  if (ix_st(i,j) == -1) cycle
  read (d_unit, err = 9650) ele%spin_taylor(i,j)%ref
  allocate (ele%spin_taylor(i,j)%term(ix_t(j)))
  do k = 1, ix_st(i,j)
    read (d_unit, err = 9700) ele%spin_taylor(i,j)%term(k)
  enddo
enddo; enddo

! If ix_lr is negative then it is a pointer to a previously read wake. 
! See write_digested_bmad_file.

if (ix_sr_long /= 0 .or. ix_sr_trans /= 0 .or. ix_lr /= 0) then
  if (ix_lr < 0) then
    call transfer_wake (ele%branch%ele(abs(ix_lr))%wake, ele%wake)

  else
    call init_wake (ele%wake, ix_sr_long, ix_sr_trans, ix_lr)
    read (d_unit, err = 9800) ele%wake%sr_file
    read (d_unit, err = 9840) ele%wake%sr_long%mode
    read (d_unit, err = 9850) ele%wake%sr_trans%mode
    read (d_unit, err = 9820) ele%wake%lr_file
    read (d_unit, err = 9830) ele%wake%lr
    read (d_unit, err = 9860) ele%wake%z_sr_max
  endif
endif

if (ix_wall3d > 0) then
  call read_this_wall3d (ele%wall3d, error)
  if (error) return
elseif (ix_wall3d < 0) then
  read (d_unit, err = 9900) idum1
  ele%wall3d => lat%branch(ix_wall3d_branch)%ele(abs(ix_wall3d))%wall3d
  if (.not. associated(ele%wall3d)) then
    call out_io(s_error$, r_name, 'ERROR IN WALL3D INIT.')
    close (d_unit)
    return
  endif
  ele%wall3d%n_link = ele%wall3d%n_link + 1
else
  read (d_unit, err = 9900) idum1
endif

!

ele%old_value = ele%value

error = .false.

return

!--------------------------------------------------------------

9100  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9110  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING K_MAX OF ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9120  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING VALUES OF ELEMENT # \i0\ ', &
                                  i_array = [ix_ele_in])
close (d_unit)
return

9140  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING EM_FIELD COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9150  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING MODE3 COMPONENT FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9200  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING WIGGLER TERM FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9350  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %R ARRAY SIZE: ' // ele%name)
close (d_unit)
return

9360  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9361  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%GRID BOUNDS FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9362  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %PHOTON%SURFACE%GRID FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9400  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING R TERM FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9500  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING DESCRIP TERM FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9600  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING AN,BN TERMS FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9650  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING %TAYLOR(:)%REF FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9700  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING TAYLOR TERMS FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9800  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%SR_FILE FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9820  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%LR_FILE FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9830  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%LR FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9840  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_long%mode FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9850  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_trans%mode FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9860  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%Z_SR_MAX FOR ELEMENT: ' // ele%name)
close (d_unit)
return

9900  continue
call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
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
  read (d_unit, iostat = ios) n_wall_section
  if (n_wall_section == 0) then
    error = .false.
    return
  endif

  if (ios /= 0) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                     'ERROR READING WALL3D N_WALL_SECTION NUMBER')
    close (d_unit)
    return
  endif

  read (d_unit, iostat = ios) wall3d(i)%ele_anchor_pt, wall3d(i)%superimpose, &
           wall3d(i)%thickness, wall3d(i)%clear_material, wall3d(i)%opaque_material

  if (ios /= 0) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                     'ERROR READING WALL PRIORITY')
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
                   sec%thickness, sec%s, sec%r0, sec%x_safe, sec%y_safe, &
                   sec%dx0_ds, sec%dy0_ds, sec%x0_coef, sec%y0_coef, sec%dr_ds, sec%p1_coef, &
                   sec%p2_coef, nv
if (ios /= 0) then
  call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', 'ERROR READING WALL3D SECTION')
  close (d_unit)
  return
endif

allocate(sec%v(nv))
do k = 1, nv
  read (d_unit, iostat = ios) sec%v(k)
  if (ios /= 0) then
    call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', 'ERROR READING WALL3D VERTEX')
    close (d_unit)
    return
  endif
enddo

end subroutine read_this_wall3d_section

end subroutine read_digested_bmad_file
