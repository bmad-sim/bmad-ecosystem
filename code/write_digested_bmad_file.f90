!+
! Subroutine write_digested_bmad_file (digested_name, lat, n_files, file_names)
!
! Subroutine to write a digested file. The names of the original files used
! to create the LAT structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! Modules Needed:
!   use bmad
!
! Input:
!     digested_name -- Character(*): Name for the digested file.
!     lat           -- lat_struct: Input lat structure.
!     n_files       -- Integer, optional: Number of original files
!     file_names(:) -- Character(*), optional: Names of the original 
!                       files used to create the lat structure.
!-

#include "CESR_platform.inc"

subroutine write_digested_bmad_file (digested_name, lat,  &
                                                  n_files, file_names)

use bmad_struct
use equality_mod, only: operator(==)
use bmad_parser_mod, except_dummy => write_digested_bmad_file

implicit none

type (lat_struct), target, intent(in) :: lat
type (branch_struct), pointer :: branch
  
real(rp) value(n_attrib_maxx)

integer, intent(in), optional :: n_files
integer d_unit, i, j, k, n_file, ix_value(n_attrib_maxx)
integer ix_wig, ix_const, ix_r(4), ix_d, ix_m, ix_t(6)
integer stat_b(24), stat, n_wake
integer ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, ierr
integer :: ix_wake(lat%n_ele_max)

character(*) digested_name
character(*), optional :: file_names(:)
character(200) fname, full_digested_name
character(32) :: r_name = 'write_digested_bmad_file'
character(30) time_stamp

logical write_wake, mode3

!! external stat ! Removed because it caused a mac link problem. DCS.

! Write input file names to the digested file
! The idea is that read_digested_bmad_file can look at these files and see
! if a file has been modified since the digested file has been created.
! Additionally record if one of the random number functions was called.

n_file = 0
if (present(n_files)) n_file = n_files

d_unit = lunget()

call fullfilename (digested_name, full_digested_name)
inquire (file = full_digested_name, name = full_digested_name)
open (unit = d_unit, file = full_digested_name, form = 'unformatted', err = 9000)

! Ran function called?

if (bp_com%ran_function_was_called) then
  write (d_unit, err = 9010) n_file+2, bmad_inc_version$
  fname = '!RAN FUNCTION WAS CALLED'
  write (d_unit) fname, 0
else
  write (d_unit, err = 9010) n_file+1, bmad_inc_version$
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
write (d_unit) fname, stat_b(10)  ! stat_b(10) = Modification date
 
! write other file names.

do j = 1, n_file
  fname = file_names(j)
  stat_b = 0
#if defined (CESR_VMS) 
  call get_file_time_stamp (fname, time_stamp)
  fname = trim(fname) // '@' // time_stamp
#else
  ierr = stat(fname, stat_b)
#endif
  write (d_unit) fname, stat_b(10)  ! stat_b(10) = Modification date
enddo

! Write the lat structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.

write (d_unit) &
        lat%name, lat%lattice, lat%input_file_name, lat%title, &
        lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track, &
        lat%n_ele_track, lat%n_ele_max, &
        lat%n_control_max, lat%n_ic_max, lat%input_taylor_order
  
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
  write (d_unit) branch%name, branch%key, branch%ix_from_branch, &
                 branch%ix_from_ele, branch%n_ele_track, branch%n_ele_max
  do j = 0, branch%n_ele_max
    call write_this_ele(branch%ele(j))
  enddo
enddo

close (d_unit)

return

! Errors

9000  print *
print *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT OPEN FILE FOR OUTPUT:'
print *, '    ', trim(digested_name)
print *, '     [This does not affect program operation]'
return

9010  print *
print *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT WRITE TO FILE FOR OUTPUT:'
print *, '    ', trim(digested_name)
print *, '     [This does not affect program operation]'
close (d_unit)
return

!-------------------------------------------------------------------------------------
contains

subroutine write_this_ele (ele)

type (ele_struct), target :: ele
type (taylor_struct), pointer :: tt(:)
type (wake_struct), pointer :: wake

integer j

!

tt => ele%taylor
    
ix_wig = 0; ix_d = 0; ix_m = 0; ix_t = 0; ix_const = 0; ix_r = 0
ix_sr_table = 0; ix_sr_mode_long = 0; ix_sr_mode_trans = 0; ix_lr = 0
mode3 = .false.

if (associated(ele%mode3))    mode3 = .true.
if (associated(ele%wig_term)) ix_wig = size(ele%wig_term)
if (associated(ele%const))    ix_const = size(ele%const)
if (associated(ele%r))        ix_r = (/ lbound(ele%r), ubound(ele%r) /)
if (associated(ele%descrip))  ix_d = 1
if (associated(ele%a_pole))   ix_m = 1
if (associated(tt(1)%term))   ix_t = (/ (size(tt(j)%term), j = 1, 6) /)

! Since some large lattices with a large number of wakes can take a lot of time writing 
! the wake info we only write a wake when needed.
! The idea is that ix_lr serves as a pointer to a previously written wake.

write_wake = .true.
if (associated(ele%wake)) then
  do j = 1, n_wake
    if (.not. lat%ele(ix_wake(j))%wake == ele%wake) cycle
    write_wake = .false.
    ix_lr = -ix_wake(j)        
  enddo

  if (write_wake) then
    if (associated(ele%wake%sr_table))      ix_sr_table      = size(ele%wake%sr_table)
    if (associated(ele%wake%sr_mode_long))  ix_sr_mode_long  = size(ele%wake%sr_mode_long)
    if (associated(ele%wake%sr_mode_trans)) ix_sr_mode_trans = size(ele%wake%sr_mode_trans)
    if (associated(ele%wake%lr))            ix_lr            = size(ele%wake%lr)
    n_wake = n_wake + 1
    ix_wake(n_wake) = i
  endif
endif

! Now write the element info

write (d_unit) mode3, ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
          ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
          ele%name, ele%type, ele%alias, ele%component_name, ele%x, ele%y, &
          ele%a, ele%b, ele%z, ele%gen0, ele%vec0, ele%mat6, &
          ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
          ele%is_on, ele%sub_key, ele%lord_status, ele%slave_status, ele%ix_value, &
          ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
          ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
          ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
          ele%num_steps, ele%integrator_order, ele%ref_orbit, &
          ele%taylor_order, ele%symplectify, ele%mode_flip, &
          ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
          ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
          ele%coupler_at, ele%on_a_girder, ele%csr_calc_on, &
          ele%map_ref_orb_in, ele%map_ref_orb_out, ele%offset_moves_aperture, &
          ele%ix_branch, ele%ref_time

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

! 

if (mode3) write (d_unit) ele%mode3

do j = 1, ix_wig
  write (d_unit) ele%wig_term(j)
enddo

if (associated(ele%const))    write (d_unit) ele%const
if (associated(ele%r))        write (d_unit) ele%r
if (associated(ele%descrip))  write (d_unit) ele%descrip
if (associated(ele%a_pole))   write (d_unit) ele%a_pole, ele%b_pole
    
do j = 1, 6
  if (ix_t(j) == 0) cycle
  write (d_unit) tt(j)%ref
  do k = 1, ix_t(j)
    write (d_unit) tt(j)%term(k)
  enddo
enddo

if (associated(ele%wake) .and. write_wake) then
  write (d_unit) ele%wake%sr_file
  write (d_unit) ele%wake%sr_table
  write (d_unit) ele%wake%sr_mode_long
  write (d_unit) ele%wake%sr_mode_trans
  write (d_unit) ele%wake%lr_file
  write (d_unit) ele%wake%lr
  write (d_unit) ele%wake%z_sr_mode_max
endif

end subroutine

end subroutine

