!+
! Subroutine write_digested_bmad_file (digested_name, lat, n_files, file_names, extra, err_flag)
!
! Subroutine to write a digested file. The names of the original files used
! to create the LAT structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! Modules Needed:
!   use bmad
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
type (ptc_parameter_struct) ptc_param

real(rp) value(num_ele_attrib$)

integer, intent(in), optional :: n_files
integer d_unit, i, j, k, n, n_file, ix_value(num_ele_attrib$), ierr
integer stat_b(24), stat, n_wake, n_wall_section
integer, allocatable :: ix_wake(:)

character(*) digested_name
character(*), optional :: file_names(:)
character(200) fname, full_digested_name
character(32) :: r_name = 'write_digested_bmad_file'
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

call fullfilename (digested_name, full_digested_name)
inquire (file = full_digested_name, name = full_digested_name)
call simplify_path (full_digested_name, full_digested_name)
open (unit = d_unit, file = full_digested_name, form = 'unformatted', err = 9000)

write (d_unit, err = 9000) n_file+1, bmad_inc_version$

! Write digested file name

stat_b = 0
ierr = stat(full_digested_name, stat_b)

fname = '!DIGESTED:' // full_digested_name
write (d_unit) fname, stat_b(2), stat_b(8), stat_b(10), 0, 0  ! stat_b(10) = Modification date
 
! write other file names.
! file names starting with '!' are not true file names but information to be stored in file.

do j = 1, n_file
  stat_b = 0
  if (file_names(j)(1:1) /= '!') then  
    call simplify_path (file_names(j), fname)
    ierr = stat(fname, stat_b)
  endif
  write (d_unit) fname, stat_b(2), stat_b(8), stat_b(10), 0, 0  ! stat_b(10) = Modification date
enddo

! Write the lat structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.

write (d_unit) &
        lat%use_name, lat%lattice, lat%input_file_name, lat%title, &
        lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track, &
        lat%n_ele_track, lat%n_ele_max, lat%lord_state, &
        lat%n_control_max, lat%n_ic_max, lat%input_taylor_order, &
        lat%absolute_time_tracking, lat%auto_scale_field_phase, lat%auto_scale_field_amp, &
        lat%use_ptc_layout, lat%pre_tracker
write (d_unit) ubound(lat%branch, 1)

! custom attribute names

if (allocated(lat%attribute_alias)) then
  n = size(lat%attribute_alias)
  write (d_unit) n
  do i = 1, n
    write (d_unit)lat%attribute_alias(i)
  enddo
else
  write (d_unit) 0
endif

! Branches

allocate (ix_wake(100))

n_wake = 0  ! number of wakes written to the digested file for this branch.
do i = 0, lat%n_ele_max
  call write_this_ele (lat%ele(i))
enddo

call write_this_wall3d (lat%branch(0)%wall3d, associated(lat%branch(0)%wall3d))

do i = 1, ubound(lat%branch, 1)
  n_wake = 0  ! number of wakes written to the digested file for this branch.
  branch => lat%branch(i)
  write (d_unit) branch%param
  write (d_unit) branch%name, branch%ix_from_branch, &
                 branch%ix_from_ele, branch%n_ele_track, branch%n_ele_max, 0
  do j = 0, branch%n_ele_max
    call write_this_ele(branch%ele(j))
  enddo
  call write_this_wall3d (branch%wall3d, associated(branch%wall3d))
enddo

! write the control info, etc

do i = 1, lat%n_control_max
  write (d_unit) lat%control(i)
enddo

do i = 1, lat%n_ic_max
  write (d_unit) lat%ic(i)
enddo

write (d_unit) lat%beam_start

! Write PTC info

call get_ptc_params (ptc_param)
write (d_unit) ptc_param

! Write random state info

if (present(extra)) then
  write (d_unit) .true.
  write (d_unit) extra
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
type (em_field_mode_struct), pointer :: mode, mode2
type (photon_surface_struct), pointer :: surf
type (surface_grid_pt_struct), pointer :: s_pt

integer ix_wig, ix_wall3d, ix_r, ix_d, ix_m, ix_t(6), ie, ib, ix_wall3d_branch
integer ix_sr_long, ix_sr_trans, ix_lr, ie_max, ix_s
integer i, j, k, n, ng, nf, n_em_field_mode, ix_ele, ix_branch, ix_wig_branch

logical write_wake, mode3

!

ix_wig = 0; ix_d = 0; ix_m = 0; ix_t = 0; ix_r = 0; ix_s = 0
ix_sr_long = 0; ix_sr_trans = 0; ix_lr = 0
mode3 = .false.; ix_wall3d = 0; n_em_field_mode = 0; ix_wig_branch = 0

if (associated(ele%mode3))          mode3 = .true.
if (associated(ele%wig))            ix_wig = size(ele%wig%term)
if (associated(ele%r))              ix_r = 1
if (associated(ele%photon))         ix_s = 1
if (associated(ele%descrip))        ix_d = 1
if (associated(ele%a_pole))         ix_m = 1
if (associated(ele%taylor(1)%term)) ix_t = [(size(ele%taylor(j)%term), j = 1, 6)]
if (associated(ele%wall3d))         ix_wall3d = size(ele%wall3d%section)
if (associated(ele%em_field))       n_em_field_mode = size(ele%em_field%mode)

! Since some large lattices with a large number of wakes can take a lot of time writing the wake info,
! we only write a wake when needed and ix_lr serves as a pointer to a previously written wake.

write_wake = .true.
if (associated(ele%wake)) then
  do j = 1, n_wake
    if (.not. ele%branch%ele(ix_wake(j))%wake == ele%wake) cycle
    write_wake = .false.
    ix_lr = -ix_wake(j)        
  enddo

  if (write_wake) then
    if (allocated(ele%wake%sr_long%mode))  ix_sr_long  = size(ele%wake%sr_long%mode)
    if (allocated(ele%wake%sr_trans%mode)) ix_sr_trans = size(ele%wake%sr_trans%mode)
    if (allocated(ele%wake%lr))       ix_lr       = size(ele%wake%lr)
    n_wake = n_wake + 1
    if (n_wake > size(ix_wake)) call re_allocate(ix_wake, 2*size(ix_wake))
    ix_wake(n_wake) = ele%ix_ele
  endif
endif

! Wiggler

if (associated(ele%wig)) then
  wig_branch_loop: do ib = 0, ele%ix_branch
    ie_max = lat%branch(ib)%n_ele_max
    if (ib == ele%ix_branch) ie_max = ele%ix_ele - 1
    do ie = 1, ie_max
      ele2 => lat%branch(ib)%ele(ie)
      if (.not. associated(ele2%wig)) cycle
      if (size(ele2%wig%term) /= size(ele%wig%term)) cycle
      if (.not. all(ele2%wig%term == ele%wig%term)) cycle
      ix_wig = -ie
      ix_wig_branch = ib
      exit wig_branch_loop
    enddo
  enddo wig_branch_loop
endif

! Wall3d

ix_wall3d_branch = 0

if (associated(ele%wall3d)) then
  wall3d_branch_loop: do ib = 0, ele%ix_branch
    ie_max = lat%branch(ib)%n_ele_max
    if (ib == ele%ix_branch) ie_max = ele%ix_ele - 1
    do ie = 1, ie_max
      ele2 => lat%branch(ib)%ele(ie)
      if (.not. associated(ele2%wall3d)) cycle
      if (size(ele2%wall3d%section) /= size(ele%wall3d%section)) cycle
      if (.not. all(ele2%wall3d%section == ele%wall3d%section)) cycle
      ix_wall3d = -ie
      ix_wall3d_branch = ib
      exit wall3d_branch_loop
    enddo
  enddo wall3d_branch_loop
endif

! Now write the element info. 
! The last zero is for future use.

write (d_unit) mode3, ix_wig, ix_wig_branch, ix_r, ix_s, ix_wall3d_branch, 0, 0, ix_d, ix_m, ix_t, &
          0, ix_sr_long, ix_sr_trans, ix_lr, ix_wall3d, n_em_field_mode, 0

write (d_unit) &
          ele%name, ele%type, ele%alias, ele%component_name, ele%x, ele%y, &
          ele%a, ele%b, ele%z, ele%gen0, ele%vec0, ele%mat6, &
          ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
          ele%is_on, ele%sub_key, ele%lord_status, ele%slave_status, ele%ix_value, &
          ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
          ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
          ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
          ele%spin_tracking_method, ele%symplectify, ele%mode_flip, &
          ele%multipoles_on, ele%taylor_map_includes_offsets, ele%Field_master, &
          ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
          ele%aperture_type, ele%csr_calc_on, ele%orientation, &
          ele%map_ref_orb_in, ele%map_ref_orb_out, ele%offset_moves_aperture, &
          ele%ix_branch, ele%ref_time, ele%scale_multipoles, 0, &
          0, ele%bookkeeping_state, ele%ptc_integration_type

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

! EM field def.

if (n_em_field_mode > 0) then

  write (d_unit) ele%em_field%mode_to_autoscale

  ix_ele = 0
  ix_branch = 0

  ! If the field is the same as the field of a previous element then set ix_ele and ix_branch.

  mode_loop: do i = 1, n_em_field_mode
    mode => ele%em_field%mode(i)
    if (.not. associated(mode%map) .and. .not. associated(mode%grid)) cycle

    do ib = 0, ele%ix_branch
      ie_max = lat%branch(ib)%n_ele_max
      if (ib == ele%ix_branch) ie_max = ele%ix_ele - 1
      do ie = 1, ie_max
        ele2 => lat%branch(ib)%ele(ie)
        if (.not. associated(ele2%em_field)) cycle
        mode2 => ele2%em_field%mode(i)
        if (size(ele2%em_field%mode) /= size(ele%em_field%mode)) cycle
        if (associated(mode%map) .and. .not. associated (mode2%map, mode%map)) cycle
        if (associated(mode%grid) .and. .not. associated (mode2%grid, mode%grid)) cycle
        ix_ele = ie
        ix_branch = ib
        exit mode_loop
      enddo  
    enddo

    exit

  enddo mode_loop

  !

  do i = 1, n_em_field_mode
    mode => ele%em_field%mode(i)

    nf = 0
    if (associated(mode%map)) nf = size(mode%map%term)
    ng = 0
    if (associated(mode%grid)) ng = size(mode%grid%pt)

    write (d_unit) nf, ng, ix_ele, ix_branch, mode%harmonic, mode%f_damp, mode%phi0_ref, &
                   mode%stored_energy, mode%m, mode%phi0_azimuth, mode%field_scale, mode%master_scale

    if (ix_ele == 0 .and. nf > 0) then
      write (d_unit) mode%map%file, mode%map%dz, mode%map%ele_anchor_pt
      write (d_unit) mode%map%term
    endif

    if (ix_ele == 0 .and. ng > 0) then
      write (d_unit) (lbound(mode%grid%pt, j), ubound(mode%grid%pt, j), j = 1, 3), &
                      mode%grid%type, mode%grid%file, mode%grid%dr, mode%grid%r0, mode%grid%ele_anchor_pt
      do j = lbound(mode%grid%pt, 3), ubound(mode%grid%pt, 3)
        write (d_unit) mode%grid%pt(:,:,j)
      enddo
    endif

  enddo 
endif

!

if (mode3) write (d_unit) ele%mode3

if (ix_wig > 0) then
  do j = 1, ix_wig
    write (d_unit) ele%wig%term(j)
  enddo
endif

if (associated(ele%r)) then
  write (d_unit) lbound(ele%r), ubound(ele%r)
  do i = lbound(ele%r, 3), ubound(ele%r, 3)
    write (d_unit) ele%r(:,:,i)
  enddo
endif

if (associated (ele%photon)) then

  surf => ele%photon%surface
  write (d_unit) ele%photon%target, ele%photon%material, &
          surf%curvature_xy, surf%has_curvature, &
          surf%grid%type, surf%grid%dr, surf%grid%r0, surf%segment, allocated(surf%grid%pt)

  if (allocated(surf%grid%pt)) then
    write (d_unit) lbound(surf%grid%pt), ubound(surf%grid%pt)
    ! Detectors do not have any grid data that needs saving
    if (ele%key /= detector$) then
      do i = lbound(surf%grid%pt, 1), ubound(surf%grid%pt, 1)
      do j = lbound(surf%grid%pt, 2), ubound(surf%grid%pt, 2)
        write (d_unit) surf%grid%pt(i,j)
      enddo
      enddo
    endif
  endif

endif

if (associated(ele%descrip))  write (d_unit) ele%descrip
if (associated(ele%a_pole))   write (d_unit) ele%a_pole, ele%b_pole
    
do j = 1, 6
  if (ix_t(j) == 0) cycle
  write (d_unit) ele%taylor(j)%ref
  do k = 1, ix_t(j)
    write (d_unit) ele%taylor(j)%term(k)
  enddo
enddo

if (associated(ele%wake) .and. write_wake) then
  write (d_unit) ele%wake%sr_file
  write (d_unit) ele%wake%sr_long%mode
  write (d_unit) ele%wake%sr_trans%mode
  write (d_unit) ele%wake%lr_file
  write (d_unit) ele%wake%lr
  write (d_unit) ele%wake%z_sr_max
endif

call write_this_wall3d (ele%wall3d, (ix_wall3d > 0))

end subroutine

!-------------------------------------------------------------------------------------
! contains

subroutine write_this_wall3d (wall3d, write_wall)

type (wall3d_struct), pointer :: wall3d
integer j, k
logical write_wall

!

if (write_wall) then

  write (d_unit) size(wall3d%section)
  write (d_unit) wall3d%ele_anchor_pt, wall3d%superimpose, &
      wall3d%thickness, wall3d%clear_material, wall3d%opaque_material

  do j = lbound(wall3d%section, 1), ubound(wall3d%section, 1)
    call write_this_wall3d_section (wall3d%section(j))
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

write (d_unit) sec%type, sec%material, sec%thickness, sec%s, sec%x0, sec%y0, &
                   sec%dx0_ds, sec%dy0_ds, sec%x0_coef, sec%y0_coef, sec%dr_ds, sec%p1_coef, &
                   sec%p2_coef, nv, sec%n_vertex_input, sec%ix_ele, sec%type
do k = 1, nv
  write (d_unit) sec%v(k)
enddo

end subroutine write_this_wall3d_section

end subroutine

