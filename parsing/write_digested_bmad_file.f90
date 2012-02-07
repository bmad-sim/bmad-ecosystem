!+
! Subroutine write_digested_bmad_file (digested_name, lat, n_files, file_names, ran_p, err_flag)
!
! Subroutine to write a digested file. The names of the original files used
! to create the LAT structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! The ran_p argument is used to encode in the digested file information about
! whether a random number generator was used to construct the lattice.
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
!   ran_p         -- ran_parsing_struct, optional: 
!     %ran_function_was_called
!     %deterministic_ran_function_was_called
!     %initial_state
!
! Output:
!   err_flag -- Logical, optional: Set True if there is a problem. EG: No write permission.
!                 Set False if everything is OK.
!-

#include "CESR_platform.inc"

subroutine write_digested_bmad_file (digested_name, lat,  n_files, file_names, ran_p, err_flag)

use bmad_struct
use equality_mod, only: operator(==)

implicit none

type (lat_struct), target, intent(in) :: lat
type (branch_struct), pointer :: branch
type (ran_parsing_struct), optional :: ran_p
type (random_state_struct), save :: dummy_ran_state

real(rp) value(n_attrib_maxx)

integer, intent(in), optional :: n_files
integer d_unit, i, j, k, n_file, ix_value(n_attrib_maxx), ierr
integer stat_b(24), stat, n_wake, n_wall_section
integer :: ix_wake(lat%n_ele_max)

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

! Ran function called?

if (present(ran_p)) then
  if (ran_p%ran_function_was_called) then
    write (d_unit, err = 9000) n_file+2, bmad_inc_version$
    fname = '!RAN FUNCTION WAS CALLED'
    write (d_unit) fname, 0
  elseif (ran_p%deterministic_ran_function_was_called) then
    write (d_unit, err = 9000) n_file+2, bmad_inc_version$
    fname = '!DETERMINISTIC RAN FUNCTION WAS CALLED'
    write (d_unit) fname, 0
  else
    write (d_unit, err = 9000) n_file+1, bmad_inc_version$
  endif
else
  write (d_unit, err = 9000) n_file+1, bmad_inc_version$
endif

! Write digested file name

stat_b = 0

#if defined (CESR_VMS) 
  call get_file_time_stamp (full_digested_name, time_stamp)
  full_digested_name = trim(full_digested_name) // '@' // time_stamp
#else
  ierr = stat(full_digested_name, stat_b)
#endif

fname = '!DIGESTED:' // full_digested_name
write (d_unit) fname, stat_b(2), stat_b(8), stat_b(10), 0, 0  ! stat_b(10) = Modification date
 
! write other file names.

do j = 1, n_file
  call simplify_path (file_names(j), fname)
  stat_b = 0
#if defined (CESR_VMS) 
  call get_file_time_stamp (fname, time_stamp)
  fname = trim(fname) // '@' // time_stamp
#else
  ierr = stat(fname, stat_b)
#endif
  write (d_unit) fname, stat_b(2), stat_b(8), stat_b(10), 0, 0  ! stat_b(10) = Modification date
enddo

! Write the lat structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.

write (d_unit) &
        lat%use_name, lat%lattice, lat%input_file_name, lat%title, &
        lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track, &
        lat%n_ele_track, lat%n_ele_max, lat%lord_status, &
        lat%n_control_max, lat%n_ic_max, lat%input_taylor_order, &
        lat%absolute_time_tracking, lat%rf_auto_scale_phase, lat%rf_auto_scale_amp, &
        lat%use_ptc_layout, lat%pre_tracker

call write_this_wall3d (lat%wall3d)

n_wake = 0  ! number of wakes written to the digested file.

do i = 0, lat%n_ele_max
  call write_this_ele (lat%ele(i))
enddo

! write the control info, etc

do i = 1, lat%n_control_max
  write (d_unit) lat%control(i)
enddo

do i = 1, lat%n_ic_max
  write (d_unit) lat%ic(i)
enddo

write (d_unit) lat%beam_start

! Write the branch line info

write (d_unit) ubound(lat%branch, 1)
do i = 1, ubound(lat%branch, 1)
  branch => lat%branch(i)
  write (d_unit) branch%param
  write (d_unit) branch%name, -1, branch%ix_from_branch, &
                 branch%ix_from_ele, branch%n_ele_track, branch%n_ele_max, branch%wall3d%anchor_pt
  do j = 0, branch%n_ele_max
    call write_this_ele(branch%ele(j))
  enddo
  call write_this_wall3d (branch%wall3d)
enddo

! Write random state info

if (present(ran_p)) then
  write (d_unit) ran_p%initial_state
else
  write (d_unit) dummy_ran_state
endif

close (d_unit)

if (present(err_flag)) err_flag = .false.

return

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
type (taylor_struct), pointer :: tt(:)
type (rf_wake_struct), pointer :: wake
type (em_field_mode_struct), pointer :: mode

integer ix_wig, n_wall_section, ix_const, ix_r, ix_d, ix_m, ix_t(6)
integer ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr
integer i, j, k, n, ng, nf, n_rf_field_mode

logical write_wake, mode3

!

tt => ele%taylor
    
ix_wig = 0; ix_d = 0; ix_m = 0; ix_t = 0; ix_const = 0; ix_r = 0
ix_sr_table = 0; ix_sr_mode_long = 0; ix_sr_mode_trans = 0; ix_lr = 0
mode3 = .false.; n_wall_section = 0; n_rf_field_mode = 0

if (associated(ele%mode3))          mode3 = .true.
if (associated(ele%wig))            ix_wig = size(ele%wig%term)
if (associated(ele%const))          ix_const = size(ele%const)
if (associated(ele%r))              ix_r = 1
if (associated(ele%descrip))        ix_d = 1
if (associated(ele%a_pole))         ix_m = 1
if (associated(tt(1)%term))         ix_t = [(size(tt(j)%term), j = 1, 6)]
if (associated(ele%wall3d%section)) n_wall_section = size(ele%wall3d%section)
if (associated(ele%em_field))       n_rf_field_mode = size(ele%em_field%mode)

! Since some large lattices with a large number of wakes can take a lot of time writing 
! the wake info we only write a wake when needed.
! The idea is that ix_lr serves as a pointer to a previously written wake.

write_wake = .true.
if (associated(ele%rf_wake)) then
  do j = 1, n_wake
    if (.not. lat%ele(ix_wake(j))%rf_wake == ele%rf_wake) cycle
    write_wake = .false.
    ix_lr = -ix_wake(j)        
  enddo

  if (write_wake) then
    if (associated(ele%rf_wake%sr_table))      ix_sr_table      = size(ele%rf_wake%sr_table)
    if (associated(ele%rf_wake%sr_mode_long))  ix_sr_mode_long  = size(ele%rf_wake%sr_mode_long)
    if (associated(ele%rf_wake%sr_mode_trans)) ix_sr_mode_trans = size(ele%rf_wake%sr_mode_trans)
    if (associated(ele%rf_wake%lr))            ix_lr            = size(ele%rf_wake%lr)
    n_wake = n_wake + 1
    ix_wake(n_wake) = ele%ix_ele
  endif
endif

! Now write the element info. 
! The last zero is for future use.

write (d_unit) mode3, ix_wig, ix_const, ix_r, 0, 0, 0, ix_d, ix_m, ix_t, &
          ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, n_wall_section, n_rf_field_mode, 0

write (d_unit) &
          ele%name, ele%type, ele%alias, ele%component_name, ele%x, ele%y, &
          ele%a, ele%b, ele%z, ele%gen0, ele%vec0, ele%mat6, &
          ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
          ele%is_on, ele%sub_key, ele%lord_status, ele%slave_status, ele%ix_value, &
          ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
          ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
          ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, ele%ref_orbit, &
          ele%spin_tracking_method, ele%symplectify, ele%mode_flip, &
          ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
          ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
          ele%aperture_type, ele%on_a_girder, ele%csr_calc_on, ele%reversed, &
          ele%map_ref_orb_in, ele%map_ref_orb_out, ele%offset_moves_aperture, &
          ele%ix_branch, ele%ref_time, ele%scale_multipoles, ele%wall3d%anchor_pt

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

! RF field def

if (n_rf_field_mode > 0) then
  do i = 1, n_rf_field_mode
    mode => ele%em_field%mode(i)
    nf = 0
    if (associated(mode%map)) nf = size(mode%map%term)
    ng = 0
    if (associated(mode%grid)) ng = size(mode%grid%pt)
    write (d_unit) nf, ng, mode%freq, mode%f_damp, mode%dphi0_ref, mode%stored_energy, &
                                    mode%m, mode%phi0_azimuth, mode%field_scale
    if (nf > 0) then
      write (d_unit) mode%map%file, mode%map%dz, mode%map%anchor_pt
      write (d_unit) mode%map%term
    endif

    if (ng > 0) then
      write (d_unit) (lbound(mode%grid%pt, j), ubound(mode%grid%pt, j), j = 1, 3), &
                      mode%grid%type, mode%grid%file, mode%grid%dr, mode%grid%r0, mode%grid%anchor_pt
      do j = lbound(mode%grid%pt, 3), ubound(mode%grid%pt, 3)
        write (d_unit) mode%grid%pt(:,:,j)
      enddo
    endif

  enddo 
endif

!

if (mode3) write (d_unit) ele%mode3

do j = 1, ix_wig
  write (d_unit) ele%wig%term(j)
enddo

if (associated(ele%const))    write (d_unit) ele%const

if (associated(ele%r)) then
  write (d_unit) lbound(ele%r), ubound(ele%r)
  do i = lbound(ele%r, 3), ubound(ele%r, 3)
    write (d_unit) ele%r(:,:,i)
  enddo
endif

if (associated(ele%descrip))  write (d_unit) ele%descrip
if (associated(ele%a_pole))   write (d_unit) ele%a_pole, ele%b_pole
    
do j = 1, 6
  if (ix_t(j) == 0) cycle
  write (d_unit) tt(j)%ref
  do k = 1, ix_t(j)
    write (d_unit) tt(j)%term(k)
  enddo
enddo

if (associated(ele%rf_wake) .and. write_wake) then
  write (d_unit) ele%rf_wake%sr_file
  write (d_unit) ele%rf_wake%sr_table
  write (d_unit) ele%rf_wake%sr_mode_long
  write (d_unit) ele%rf_wake%sr_mode_trans
  write (d_unit) ele%rf_wake%lr_file
  write (d_unit) ele%rf_wake%lr
  write (d_unit) ele%rf_wake%z_sr_mode_max
endif

call write_this_wall3d (ele%wall3d)

end subroutine

!-------------------------------------------------------------------------------------
! contains

subroutine write_this_wall3d (wall3d)

type (wall3d_struct) wall3d
type (wall3d_section_struct), pointer :: sec
type (wall3d_vertex_struct), pointer :: v

integer j, k

!

if (associated(wall3d%section)) then

  write (d_unit) size(wall3d%section)

  do j = lbound(wall3d%section, 1), ubound(wall3d%section, 1)
    sec => wall3d%section(j)
    write (d_unit) sec%type, sec%s, sec%x0, sec%y0, sec%dx0_ds, sec%dy0_ds, sec%x0_coef, sec%y0_coef, &
                   sec%dr_ds, sec%p1_coef, sec%p2_coef, size(sec%v), sec%n_vertex_input
    do k = 1, size(sec%v)
      write (d_unit) sec%v(k)
    enddo
  enddo

else
  write (d_unit) 0
endif

end subroutine

end subroutine

