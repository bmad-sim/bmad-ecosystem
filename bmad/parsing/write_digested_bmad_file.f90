!+
! Subroutine write_digested_bmad_file (digested_name, lat, n_files, file_names, extra, err_flag)
!
! Subroutine to write a digested file. The names of the original files used
! to create the LAT structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! Input:
!   digested_name -- Character(*): Name for the digested file.
!   lat           -- lat_struct: Input lat structure.
!   n_files       -- Integer, optional: Number of original files
!   file_names(:) -- Character(*), optional: Names of the original 
!                     files used to create the lat structure.
!   extra         -- extra_parsing_info_struct, optional: Extra info that can
!                     be stored in the digested file.
!
! Output:
!   err_flag -- Logical, optional: Set True if there is a problem. EG: No write permission.
!                 Set False if everything is OK.
!-

subroutine write_digested_bmad_file (digested_name, lat,  n_files, file_names, extra, err_flag)

use equality_mod, only: operator(==)
use ptc_interface_mod, dummy => write_digested_bmad_file

implicit none

type (lat_struct), target, intent(in) :: lat
type (branch_struct), pointer :: branch
type (extra_parsing_info_struct), optional :: extra
type (wake_struct), pointer :: wake
type (control_struct), pointer :: ctl

real(rp) value(num_ele_attrib$)

integer, intent(in), optional :: n_files
integer d_unit, i, j, k, n, nk, n_file, ix_value(num_ele_attrib$), ierr
integer stat_b(24), stat, n_wake, n_wall_section, n_custom, n_print
integer, allocatable :: ix_ele_wake(:)
integer, allocatable :: index_list(:)

character(*) digested_name
character(*), optional :: file_names(:)
character(200) fname, full_digested_name
character(100), allocatable :: name_list(:)
character(*), parameter :: r_name = 'write_digested_bmad_file'
character(30) time_stamp

logical, optional :: err_flag
logical is_open

!! external stat ! Removed because it caused a mac link problem. DCS.

! Write input file names to the digested file
! The idea is that read_digested_bmad_file can look at these files and see
! if a file has been modified since the digested file has been created.
! Additionally record if one of the random number functions was called.

if (present(err_flag)) err_flag = .true.
n_file = 0
if (present(n_files)) n_file = n_files

d_unit = lunget()

call fullfilename (digested_name, fname)
inquire (file = fname, name = full_digested_name)
call simplify_path (full_digested_name, full_digested_name)
open (unit = d_unit, file = full_digested_name, form = 'unformatted', err = 9000)

write (d_unit, err = 9000) n_file+1, bmad_inc_version$

! Write digested file name

stat_b = 0
ierr = stat(full_digested_name, stat_b)

fname = '!DIGESTED:' // full_digested_name
write (d_unit) fname, stat_b(2), stat_b(8), stat_b(10) ! Inode, Size, Modification date
 
! write other file names.
! file names starting with '!' are not true file names but information to be stored in file.

do j = 1, n_file
  stat_b = 0
  if (file_names(j)(1:1) /= '!') then  
    call simplify_path (file_names(j), fname)
    ierr = stat(fname, stat_b)
  endif
  write (d_unit) fname, stat_b(2), stat_b(8), stat_b(10) ! Inode, Size, Modification date
enddo

! Write the lat structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.
! Note: Set lat%ramper_slave_bookkeeping_done = False since ramper pointer not in digested file.
n_custom = -1
if (allocated(lat%custom)) n_custom = size(lat%custom)
n_print = -1
if (allocated(lat%print_str)) n_print = size(lat%print_str)
write (d_unit) lat%use_name, lat%machine, lat%lattice, lat%input_file_name, lat%title
write (d_unit) lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track
write (d_unit) lat%n_ele_track, lat%n_ele_max, lat%lord_state, lat%n_control_max, lat%n_ic_max
write (d_unit) lat%input_taylor_order, lat%photon_type, .false.
write (d_unit) ubound(lat%branch, 1), lat%pre_tracker, n_custom, n_print

! Global custom

if (n_custom /= -1) write(d_unit) lat%custom
if (n_print /= -1) write (d_unit) lat%print_str

! Defined constants and custom attributes

call custom_ele_attrib_name_list(index_list, name_list)
write (d_unit) size(index_list)
do i = 1, size(index_list)
  write (d_unit) index_list(i), name_list(i)
enddo

if (allocated(lat%constant)) then
  write (d_unit) size(lat%constant)
  do i = 1, size(lat%constant)
    write (d_unit) lat%constant(i)
  enddo
else
  write (d_unit) 0
endif

! Branches

allocate (ix_ele_wake(100))

n_wake = 0  ! number of wakes written to the digested file for this branch.
do i = 0, lat%n_ele_max
  call write_this_ele (lat%ele(i))
enddo

call write_this_wall3d (lat%branch(0)%wall3d, associated(lat%branch(0)%wall3d))

write (d_unit) lat%branch(0)%name

do i = 1, ubound(lat%branch, 1)
  n_wake = 0  ! number of wakes written to the digested file for this branch.
  branch => lat%branch(i)
  write (d_unit) branch%param
  write (d_unit) branch%name, branch%ix_from_branch, branch%ix_to_ele, &
                 branch%ix_from_ele, branch%n_ele_track, branch%n_ele_max
  do j = 0, branch%n_ele_max
    call write_this_ele(branch%ele(j))
  enddo
  call write_this_wall3d (branch%wall3d, associated(branch%wall3d))
enddo

! write the control info, etc

do i = 1, lat%n_control_max
  ctl => lat%control(i)
  n = 0
  if (allocated(ctl%stack)) n = size(ctl%stack)
  nk = 0
  if (allocated(ctl%y_knot)) nk = size(ctl%y_knot)
  write (d_unit) n, nk, ctl%value, ctl%lord, ctl%slave, ctl%ix_attrib, ctl%attribute, ctl%slave_name
  do j = 1, n
    write (d_unit) ctl%stack(j)
  enddo
  if (nk /= 0) write (d_unit) ctl%y_knot
enddo

do i = 1, lat%n_ic_max
  write (d_unit) lat%ic(i)
enddo

write (d_unit) lat%particle_start
write (d_unit) lat%beam_init

! Write common stuff


if (present(extra)) then
  call ran_default_state (get_state = extra%ran_state) ! Get random state.
  write (d_unit) .true.
  write (d_unit) extra
  write (d_unit) bmad_com
  write (d_unit) space_charge_com
  write (d_unit) ptc_com%max_fringe_order, ptc_com%exact_model, ptc_com%exact_misalign, &
          ptc_com%vertical_kick, ptc_com%cut_factor, ptc_com%old_integrator, &
          ptc_com%use_orientation_patches, ptc_com%print_info_messages, ptc_com%exact_misalign
else
  write (d_unit) .false.
endif

! End stuff

close (d_unit)

if (present(err_flag)) err_flag = .false.

return

!------------------------------------------------------
! Errors

9000  continue
call out_io (s_warn$, r_name, &
               'NOTE: CANNOT OPEN FILE FOR OUTPUT:', &
               '    ' // trim(digested_name), &
               '     [This does not affect program operation]')
inquire (d_unit, opened = is_open)
if (is_open) close (d_unit)
return

!-------------------------------------------------------------------------------------
contains

subroutine write_this_ele (ele)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele2
type (wake_struct), pointer :: wake
type (photon_element_struct), pointer :: ph
type (photon_reflect_table_struct), pointer :: prt
type (surface_grid_pt_struct), pointer :: s_pt
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_struct), pointer :: ct_map
type (gen_grad_map_struct), pointer :: gg_map
type (gen_grad1_struct), pointer :: ggcoef
type (grid_field_struct), pointer :: g_field
type (ac_kicker_struct), pointer :: ac
type (converter_distribution_struct), pointer :: c_dist
type (converter_prob_pc_r_struct), pointer :: ppcr
type (converter_direction_out_struct), pointer :: c_dir
type (control_ramp1_struct), pointer ::rmp

integer ix_wall3d, ix_r, ix_d, ix_m, ix_e, ix_t(6), ix_st(0:3), ie, ib, ix_wall3d_branch
integer ix_sr_long, ix_sr_trans, ix_sr_time, ix_lr_mode, ie_max, ix_s, n_var, ix_ptr, im, n1, n2
integer i, j, k, n, nr, n_gen, n_grid, n_cart, n_cyl, ix_ele, ix_c, ix_branch
integer n_cus, ix_convert, n_energy, n_angle, n_foil

logical write_wake, mode3

!

ix_d = 0; ix_m = 0; ix_e = 0; ix_t = -1; ix_r = 0; ix_s = 0
ix_sr_long = 0; ix_sr_trans = 0; ix_sr_time, ix_lr_mode = 0; ix_st = -1
mode3 = .false.; ix_wall3d = 0; ix_convert = 0; ix_c = 0
n_cart = 0; n_gen = 0; n_grid = 0; n_cyl = 0; n_cus = 0; n_foil = 0

if (associated(ele%mode3))             mode3 = .true.
if (associated(ele%cartesian_map))     n_cart = size(ele%cartesian_map)
if (associated(ele%cylindrical_map))   n_cyl = size(ele%cylindrical_map)
if (associated(ele%gen_grad_map))      n_gen = size(ele%gen_grad_map)
if (associated(ele%grid_field))        n_grid = size(ele%grid_field)
if (associated(ele%custom))            n_cus = size(ele%custom)
if (associated(ele%foil))              n_foil = size(ele%foil%material)
if (associated(ele%converter))         ix_convert = 1
if (associated(ele%r))                 ix_r = 1
if (associated(ele%photon))            ix_s = 1
if (associated(ele%descrip))           ix_d = 1
if (associated(ele%a_pole))            ix_m = 1
if (associated(ele%a_pole_elec))       ix_e = 1
if (associated(ele%control))           ix_c = 1
do n = 1, size(ele%taylor)
  if (associated(ele%taylor(n)%term)) ix_t(n) = size(ele%taylor(n)%term)
enddo
do i = 0, 3
  if (associated(ele%spin_taylor(i)%term)) ix_st(i) = size(ele%spin_taylor(i)%term)
enddo
if (associated(ele%wall3d))     ix_wall3d = size(ele%wall3d)

! Since some large lattices with a large number of wakes can take a lot of time writing the wake info,
! we only write a wake when needed and ix_lr_mode serves as a pointer to a previously written wake.

write_wake = .true.
wake => ele%wake
if (associated(wake)) then
  do j = 1, n_wake
    if (.not. ele%branch%ele(ix_ele_wake(j))%wake == wake) cycle
    write_wake = .false.
    ix_lr_mode = -ix_ele_wake(j)        
  enddo

  if (write_wake) then
    if (allocated(wake%sr%long))      ix_sr_long    = size(wake%sr%long)
    if (allocated(wake%sr%trans))     ix_sr_trans   = size(wake%sr%trans)
    if (allocated(wake%sr%time))      ix_sr_time    = size(wake%sr%time)
    if (allocated(wake%lr%mode))      ix_lr_mode    = size(wake%lr%mode)
    n_wake = n_wake + 1
    if (n_wake > size(ix_ele_wake)) call re_allocate(ix_ele_wake, 2*size(ix_ele_wake))
    ix_ele_wake(n_wake) = ele%ix_ele
  endif
endif

! Wall3d
! Go through all the elements before the current element and if some element has the same
! wall3d model then just set ix_wall3d, ix_wall3d_branch to point to this element.
! Negative sign for ix_wall3d will signal to read_digested_bmad_file that this is the situation.

ix_wall3d_branch = 0

if (associated(ele%wall3d)) then
  wall3d_branch_loop: do ib = 0, ele%ix_branch
    ie_max = lat%branch(ib)%n_ele_max
    if (ib == ele%ix_branch) ie_max = ele%ix_ele - 1
    wall_ele_loop: do ie = 1, ie_max
      ele2 => lat%branch(ib)%ele(ie)
      if (.not. associated(ele2%wall3d)) cycle
      if (size(ele2%wall3d) /= size(ele%wall3d)) cycle
      do j = 1, size(ele%wall3d)
        if (size(ele2%wall3d(j)%section) /= size(ele%wall3d(j)%section)) cycle wall_ele_loop
        if (.not. all(ele2%wall3d(j)%section == ele%wall3d(j)%section)) cycle wall_ele_loop
      enddo
      ix_wall3d = -ie
      ix_wall3d_branch = ib
      exit wall3d_branch_loop
    enddo wall_ele_loop
  enddo wall3d_branch_loop
endif

! Now write the element info. 
! The last zero is for future use.

write (d_unit) mode3, ix_r, ix_s, ix_wall3d_branch, associated(ele%ac_kick), associated(ele%rad_map), &
          ix_convert, ix_d, ix_m, ix_t, ix_st, ix_e, ix_sr_long, ix_sr_trans, &
          ix_lr_mode, ix_wall3d, ix_c, n_cart, n_cyl, n_gen, n_grid, n_foil, n_cus, ix_convert

write (d_unit) &
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

! This compresses the ele%value array

k = 0
do j = 1, size(ele%value)
  if (ele%value(j) == 0) cycle
  k = k + 1
  value(k) = ele%value(j)
  ix_value(k) = j
  enddo
write (d_unit) k
write (d_unit) ix_value(1:k), value(1:k)

! Control vars

if (ix_c == 1) then
  n_var = -1; nk = -1; nr = -1; ix_r = -1
  if (allocated(ele%control%var)) n_var = size(ele%control%var)
  if (allocated(ele%control%x_knot)) nk = size(ele%control%x_knot)
  if (allocated(ele%control%ramp)) nr = size(ele%control%ramp)
  if (allocated(ele%control%ramper_lord)) ix_r = size(ele%control%ramper_lord)
  write (d_unit) n_var, nk, nr, ix_r

  if (nk > -1) write (d_unit) ele%control%x_knot

  do i = 1, n_var
    write (d_unit) ele%control%var(i)
  enddo

  do i = 1, nr
    rmp => ele%control%ramp(i)
    n = 0
    if (allocated(rmp%stack)) n = size(rmp%stack)
    nk = 0
    if (allocated(rmp%y_knot)) nk = size(rmp%y_knot)
    write (d_unit) rmp%slave_name, n, nk, rmp%attribute, rmp%is_controller
    do j = 1, n
      write (d_unit) rmp%stack(j)
    enddo
    if (nk /= 0) write (d_unit) rmp%y_knot
  enddo

  if (ix_r > -1) then
    write (d_unit) ele%control%ramper_lord%ix_ele
    write (d_unit) ele%control%ramper_lord%ix_con
  endif
endif

! AC_kicker

if (associated(ele%ac_kick)) then
  ac => ele%ac_kick
  n1 = -1; n2 = -1
  if (allocated(ac%amp_vs_time)) n1 = size(ac%amp_vs_time)
  if (allocated(ac%frequency)) n2 = size(ac%frequency)
  write (d_unit) n1, n2

  if (allocated(ac%amp_vs_time)) then
    do n = lbound(ac%amp_vs_time, 1), ubound(ac%amp_vs_time, 1)
      write (d_unit) ac%amp_vs_time(n)
    enddo
  endif

  if (allocated(ac%frequency)) then
    do n = lbound(ac%frequency, 1), ubound(ac%frequency, 1)
      write (d_unit) ac%frequency(n)
    enddo
  endif
endif

! Foil

if (associated(ele%foil)) then
  do n = 1, size(ele%foil%material)
    write (d_unit) ele%foil%material(n)
  enddo
endif

! Converter

if (associated(ele%converter)) then
  write (d_unit) ele%converter%species_out, ele%converter%material_type, size(ele%converter%dist)
  do n = 1, size(ele%converter%dist)
    c_dist => ele%converter%dist(n)
    write (d_unit) c_dist%thickness, size(c_dist%sub_dist)
    do j = 1, size(c_dist%sub_dist)
      write (d_unit) c_dist%sub_dist(j)%pc_in
      ppcr => c_dist%sub_dist(j)%prob_pc_r
      write (d_unit) ppcr%integrated_prob, size(ppcr%pc_out), size(ppcr%r)
      write (d_unit) ppcr%pc_out
      write (d_unit) ppcr%r
      write (d_unit) ppcr%prob
      write (d_unit) ppcr%spin_z
      c_dir => c_dist%sub_dist(j)%dir_out
      write (d_unit) [size(c_dir%beta%fit_1D_r), size(c_dir%alpha_x%fit_1D_r), size(c_dir%alpha_y%fit_1D_r), &
                  size(c_dir%c_x%fit_1D_r), size(c_dir%dxds_min%fit_1D_r), size(c_dir%dxds_max%fit_1D_r), &
                  size(c_dir%dyds_max%fit_1D_r)]
      write (d_unit) c_dir%beta%fit_1d_r, c_dir%beta%fit_2d_pc, c_dir%beta%fit_2d_r, c_dir%beta%c0
      write (d_unit) c_dir%alpha_x%fit_1d_r, c_dir%alpha_x%fit_2d_pc, c_dir%alpha_x%fit_2d_r, c_dir%alpha_x%c0
      write (d_unit) c_dir%alpha_y%fit_1d_r, c_dir%alpha_y%fit_2d_pc, c_dir%alpha_y%fit_2d_r, c_dir%alpha_y%c0
      write (d_unit) c_dir%c_x%fit_1d_r, c_dir%c_x%fit_2d_pc, c_dir%c_x%fit_2d_r, c_dir%c_x%c0
      write (d_unit) c_dir%dxds_min%fit_1d_r, c_dir%dxds_min%fit_2d_pc, c_dir%dxds_min%fit_2d_r, c_dir%dxds_min%c0
      write (d_unit) c_dir%dxds_max%fit_1d_r, c_dir%dxds_max%fit_2d_pc, c_dir%dxds_max%fit_2d_r, c_dir%dxds_max%c0
      write (d_unit) c_dir%dyds_max%fit_1d_r, c_dir%dyds_max%fit_2d_pc, c_dir%dyds_max%fit_2d_r, c_dir%dyds_max%c0
    enddo
  enddo
endif

! Cartesian map
! If the field is the same as the field of a previous element then just reference the prior element's field

do i = 1, n_cart
  ct_map => ele%cartesian_map(i)

  if (ct_map%ptr%file == '') then
    call out_io (s_error$, r_name, 'CARTESIAN_MAP FILE REFERENCE IS BLANK!!??. PLEASE REPORT THIS.')
    write (ct_map%ptr%file, '(3i0)') ele%ix_branch, ele%ix_ele, i  ! Something unique
  endif

  write (d_unit) ct_map%field_scale, ct_map%master_parameter, ct_map%ele_anchor_pt, ct_map%field_type, ct_map%r0

  call find_matching_fieldmap (ct_map%ptr%file, ele, cartesian_map$, ele2, ix_ptr) 
  if (ix_ptr > 0) then
    write (d_unit) ele2%ix_ele, ele2%ix_branch, ix_ptr, size(ct_map%ptr%term)
  else
    write (d_unit) -1, -1, -1, size(ct_map%ptr%term)
    write (d_unit) ct_map%ptr%file
    do j = 1, size(ct_map%ptr%term)
      write (d_unit) ct_map%ptr%term(j)
    enddo
  endif
enddo

! Cylindrical map
! If the field is the same as the field of a previous element then just reference the prior element's field

do i = 1, n_cyl
  cl_map => ele%cylindrical_map(i)

  if (cl_map%ptr%file == '') then
    call out_io (s_error$, r_name, 'CYLINDRICAL_MAP FILE REFERENCE IS BLANK!!??. PLEASE REPORT THIS.')
    write (cl_map%ptr%file, '(3i0)') ele%ix_branch, ele%ix_ele, i  ! Something unique
  endif

  write (d_unit) cl_map%field_scale, cl_map%master_parameter, cl_map%harmonic, &
                cl_map%phi0_fieldmap, cl_map%theta0_azimuth, cl_map%ele_anchor_pt, cl_map%m, cl_map%dz, cl_map%r0

  call find_matching_fieldmap (cl_map%ptr%file, ele, cylindrical_map$, ele2, ix_ptr) 
  if (ix_ptr > 0) then
    write (d_unit) ele2%ix_ele, ele2%ix_branch, ix_ptr, size(cl_map%ptr%term)
  else
    write (d_unit) -1, -1, -1, size(cl_map%ptr%term)
    write (d_unit) cl_map%ptr%file
    do j = 1, size(cl_map%ptr%term)
      write (d_unit) cl_map%ptr%term(j)
    enddo
  endif
enddo

! Gen_grad_field

do i = 1, n_gen
  gg_map => ele%gen_grad_map(i)

  write (d_unit) gg_map%field_scale, gg_map%master_parameter, gg_map%curved_ref_frame, &
          gg_map%ele_anchor_pt, gg_map%field_type, gg_map%dz, gg_map%r0, size(gg_map%gg), gg_map%iz0, gg_map%iz1

  do j = 1, size(gg_map%gg)
    ggcoef => gg_map%gg(j)
    write (d_unit) ggcoef%m, ggcoef%sincos, ggcoef%n_deriv_max, ubound(ggcoef%deriv,2)
    do k = gg_map%iz0, gg_map%iz1
      write (d_unit) ggcoef%deriv(k,:)
    enddo
  enddo
enddo

! Grid field
! If the field is the same as the field of a previous element then just reference the prior element's field

do i = 1, n_grid
  g_field => ele%grid_field(i)

  if (g_field%ptr%file == '') then
    call out_io (s_error$, r_name, 'GRID_FIELD FILE REFERENCE IS BLANK!!??. PLEASE REPORT THIS.')
    write (g_field%ptr%file, '(3i0)') ele%ix_branch, ele%ix_ele, i  ! Something unique
  endif

  write (d_unit) g_field%field_scale, g_field%master_parameter, &
                g_field%ele_anchor_pt, g_field%phi0_fieldmap, g_field%dr, &
                g_field%r0, g_field%harmonic, g_field%geometry, &
                g_field%curved_ref_frame, g_field%field_type, g_field%interpolation_order

  call find_matching_fieldmap (g_field%ptr%file, ele, grid_field$, ele2, ix_ptr) 
  if (ix_ptr > 0) then
    write (d_unit) ele2%ix_ele, ele2%ix_branch, ix_ptr, lbound(g_field%ptr%pt), ubound(g_field%ptr%pt)
  else
    write (d_unit) -1, -1, -1, lbound(g_field%ptr%pt), ubound(g_field%ptr%pt)
    write (d_unit) g_field%ptr%file
    do j = lbound(g_field%ptr%pt, 3), ubound(g_field%ptr%pt, 3)
      write (d_unit) g_field%ptr%pt(:, :, j)
    enddo
  endif
enddo

!

if (mode3) write (d_unit) ele%mode3

if (associated(ele%r)) then
  write (d_unit) lbound(ele%r), ubound(ele%r)
  do i = lbound(ele%r, 3), ubound(ele%r, 3)
    write (d_unit) ele%r(:,:,i)
  enddo
endif

if (associated(ele%custom)) then
  write (d_unit) ele%custom
endif

if (associated (ele%photon)) then
  ph => ele%photon
  write (d_unit) ph%target, ph%material, ph%curvature, ph%grid%active, ph%grid%type, &
    ph%grid%dr, ph%grid%r0, allocated(ph%grid%pt), ph%pixel%dr, ph%pixel%r0, allocated(ph%pixel%pt), &
    allocated(ph%reflectivity_table_sigma%angle), allocated(ph%reflectivity_table_pi%angle), allocated(ph%init_energy_prob)

  if (allocated(ph%grid%pt)) then
    write (d_unit) lbound(ph%grid%pt), ubound(ph%grid%pt)
    do i = lbound(ph%grid%pt, 1), ubound(ph%grid%pt, 1)
    do j = lbound(ph%grid%pt, 2), ubound(ph%grid%pt, 2)
      write (d_unit) ph%grid%pt(i,j)
    enddo
    enddo
  endif

  if (allocated(ph%pixel%pt)) then
    write (d_unit) lbound(ph%pixel%pt), ubound(ph%pixel%pt)
    ! Detectors do not have any pixel data that needs saving
  endif

  if (allocated(ph%reflectivity_table_sigma%angle)) then
    prt => ph%reflectivity_table_sigma
    n_energy = size(prt%energy)
    n_angle  = size(prt%angle)
    write (d_unit) n_energy, n_angle
    write (d_unit) prt%angle
    write (d_unit) prt%energy
    write (d_unit) prt%bragg_angle
    do j = 1, n_energy
      write (d_unit) prt%p_reflect(:,j)
    enddo
  endif

  if (allocated(ph%reflectivity_table_pi%angle)) then
    prt => ph%reflectivity_table_pi
    n_energy = size(prt%energy)
    n_angle  = size(prt%angle)
    write (d_unit) n_energy, n_angle
    write (d_unit) prt%angle
    write (d_unit) prt%energy
    write (d_unit) prt%bragg_angle
    do j = 1, n_energy
      write (d_unit) prt%p_reflect(:,j)
    enddo
  endif

  if (allocated(ph%init_energy_prob)) then
    write (d_unit) size(ph%init_energy_prob)
    write (d_unit) ph%init_energy_prob
    write (d_unit) ph%integrated_init_energy_prob
  endif
endif

if (associated(ele%descrip))      write (d_unit) ele%descrip
if (associated(ele%a_pole))       write (d_unit) ele%a_pole, ele%b_pole
if (associated(ele%a_pole_elec))  write (d_unit) ele%a_pole_elec, ele%b_pole_elec
    
do j = 1, size(ele%taylor)
  if (.not. associated(ele%taylor(j)%term)) cycle
  write (d_unit) ele%taylor(j)%ref
  do k = 1, size(ele%taylor(j)%term)
    write (d_unit) ele%taylor(j)%term(k)
  enddo
enddo

do i = 0, 3
  if (.not. associated(ele%spin_taylor(i)%term)) cycle
  write (d_unit) ele%spin_taylor(i)%ref
  do k = 1, size(ele%spin_taylor(i)%term)
    write (d_unit) ele%spin_taylor(i)%term(k)
  enddo
enddo

if (associated(wake) .and. write_wake) then
  write (d_unit) wake%sr%z_ref_long, wake%sr%z_ref_trans, wake%sr%z_max, wake%sr%scale_with_length, wake%sr%amp_scale, wake%sr%z_scale
  do i = 1, size(wake%sr%long)
    write (d_unit) wake%sr%long(i)
  enddo
  do i = 1, size(wake%sr%trans)
    write (d_unit) wake%sr%trans(i)
  enddo
  write (d_unit) wake%lr%t_ref, wake%lr%freq_spread, wake%lr%self_wake_on, wake%lr%amp_scale, wake%lr%time_scale
  do i = 1, size(wake%lr%mode)
    write (d_unit) wake%lr%mode(i)
  enddo
endif

call write_this_wall3d (ele%wall3d, (ix_wall3d > 0))

if (associated(ele%rad_map)) then
  write (d_unit) ele%rad_map%rm0
  write (d_unit) ele%rad_map%rm1, ele%rad_map%stale
endif

end subroutine

!-------------------------------------------------------------------------------------
! contains

subroutine write_this_wall3d (wall3d, write_wall)

type (wall3d_struct), pointer :: wall3d(:)
integer i, j, k
logical write_wall

!

if (write_wall) then

  write (d_unit) size(wall3d)
  do i = 1, size(wall3d)
    write (d_unit) size(wall3d(i)%section), wall3d(i)%type, &
                   wall3d(i)%ele_anchor_pt, wall3d(i)%superimpose, &
                   wall3d(i)%thickness, wall3d(i)%clear_material, wall3d(i)%opaque_material

    do j = lbound(wall3d(i)%section, 1), ubound(wall3d(i)%section, 1)
      call write_this_wall3d_section (wall3d(i)%section(j))
    enddo
  enddo

else
  write (d_unit) 0
endif

end subroutine

!-------------------------------------------------------------------------------------
! contains

subroutine write_this_wall3d_section (sec)

type (wall3d_section_struct), target :: sec
type (wall3d_vertex_struct), pointer :: v
integer nv, k

!

if (allocated(sec%v)) then
  nv = size(sec%v)
else
  nv = 0
endif

write (d_unit) sec%name, sec%material, sec%type, sec%n_vertex_input, &
                   sec%ix_ele, sec%ix_branch, sec%patch_in_region, &
                   sec%thickness, sec%s, sec%r0, sec%dx0_ds, sec%dy0_ds, &
                   sec%x0_coef, sec%y0_coef, sec%dr_ds, sec%p1_coef, sec%p2_coef, nv
do k = 1, nv
  write (d_unit) sec%v(k)
enddo

end subroutine write_this_wall3d_section

end subroutine

