!+
! Subroutine read_digested_bmad_file (digested_name, lat, version)
!
! Subroutine to read in a digested file. The subroutine will check that
! the version of the digested file is up to date and that the digested file
! is current with respect to the original BMAD files that were used. [See
! write_digested_bmad_file.]
!
! Note: This subroutine also reads in the common structures for BMAD_PARSER2
!
! Modules Needed:
!   use bmad
!
! Input:
!   digested_name -- Character(*): Name of the digested file.
!
! Output:
!   lat         -- lat_struct: Output lattice structure
!   version     -- Integer: bmad_inc_version number stored in the lattice file.
!                   If the file is current this number should be the same 
!                   as the global parameter bmad_inc_version$.
!   bmad_status -- Bmad_status_struct: Common block status structure.
!     %ok           -- Set .false. if read failure.
!-

#include "CESR_platform.inc"

subroutine read_digested_bmad_file (digested_name, lat, version)

use bmad_struct
use bmad_interface, except_dummy => read_digested_bmad_file
use multipole_mod

implicit none

type old_mode_info_struct
  real(rp) tune      ! "fractional" tune in radians: 0 < tune < 2pi
  real(rp) emit      ! Emittance.
  real(rp) chrom     ! Chromaticity.
end type
type (old_mode_info_struct) old_a, old_b, old_z

type (lat_struct), target, intent(inout) :: lat
type (branch_struct), pointer :: branch

real(rp) value(n_attrib_maxx)

integer d_unit, n_files, version, i, j, k, ix, ix_value(n_attrib_maxx)
integer ix_wig, ix_const, ix_r(4), ix_d, ix_m, ix_t(6), ios, k_max
integer ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, ierr, stat
integer stat_b(12), idate_old, n_branch, n, control_type

character(*) digested_name
character(200) fname_read, fname_versionless, fname_full
character(200) input_file_name, full_digested_name
character(200), allocatable :: file_names(:)
character(25) :: r_name = 'read_digested_bmad_file'
character(40) old_time_stamp, new_time_stamp

logical found_it, v86, v87, v88, v89, v_old, mode3, error, is_open

! init all elements in lat

call init_ele (lat%ele_init)  ! init pointers
call deallocate_lat_pointers (lat)

! Read the digested file.
! Some old versions can be read even though they are not the current version.

d_unit = lunget()
bmad_status%ok = .true.
lat%n_ele_track = 0

call fullfilename (digested_name, full_digested_name)
inquire (file = full_digested_name, name = full_digested_name)
call simplify_path (full_digested_name, full_digested_name)
open (unit = d_unit, file = full_digested_name, status = 'old',  &
                     form = 'unformatted', action = 'READ', err = 9000)

read (d_unit, err = 9010) n_files, version

v86 = (version == 86)
v87 = (version == 87)
v88 = (version == 88)
v89 = (version == 89)

v_old = v86 .or. v87 .or. v88

if (version < bmad_inc_version$) then
  if (bmad_status%type_out) call out_io (s_warn$, r_name, &
         (/ 'DIGESTED FILE VERSION OUT OF DATE \i0\ > \i0\ ' /),  &
          i_array = (/ bmad_inc_version$, version /) )
  if (v_old) then 
    allocate (file_names(n_files))
    bmad_status%ok = .false.
  else
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif
endif

if (version > bmad_inc_version$) then
  if (bmad_status%type_out) call out_io (s_warn$, r_name, &
       'DIGESTED FILE HAS VERSION: \i0\ ', &
       'GREATER THAN VERSION OF THIS PROGRAM: \i0\ ', &
       'WILL NOT USE THE DIGESTED FILE. YOU SHOULD RECOMPILE THIS PROGRAM.', &
       i_array = (/ version, bmad_inc_version$ /) )
  close (d_unit)
  bmad_status%ok = .false.
  return
endif

! if the digested file is out of date then we still read in the file since
! we can possibly reuse the taylor series.

do i = 1, n_files

  old_time_stamp = ''
  new_time_stamp = ''
  stat_b = 0

  read (d_unit, err = 9020) fname_read, idate_old

#if defined (CESR_VMS) 
  ix = index(fname_read, '@')
  if (ix /= 0) then
    old_time_stamp = fname_read(ix+1:)
    fname_read = fname_read(:ix-1)
  endif
#endif

  if (fname_read(1:10) == '!DIGESTED:') then
    fname_read = fname_read(11:)
    inquire (file = fname_read, opened = is_open)
    if (.not. is_open) then
      if (bmad_status%type_out .and. bmad_status%ok) call out_io(s_warn$, &
                                          r_name, ' NOTE: MOVED DIGESTED FILE.')
      bmad_status%ok = .false.
    endif
    cycle
  endif

  if (fname_read == '!RAN FUNCTION WAS CALLED') then
    if (bmad_status%type_out) call out_io(s_warn$, r_name, &
                'NOTE: THE RANDOM NUMBER FUNCTION WAS USED IN THE LATTICE FILE SO THIS', &
                '      LATTICE WILL DIFFER FROM OTHER LATTICES GENERATED FROM THE SAME FILE.')
    bmad_status%ok = .false.
    cycle
  endif

  call simplify_path (fname_read, fname_read)
  if (v_old) file_names(i) = fname_read  ! fake out

#if defined (CESR_VMS) 
  ix = index(fname_read, ';')
  fname_versionless = fname_read(:ix-1)
  call get_file_time_stamp (fname_versionless, new_time_stamp)
#else
  ierr = stat(fname_read, stat_b)
  fname_versionless = fname_read
#endif

  inquire (file = fname_versionless, exist = found_it, name = fname_full)
  call simplify_path (fname_full, fname_full)
  if (.not. found_it .or. fname_read /= fname_full .or. stat_b(10) /= idate_old .or. &
                                                    old_time_stamp /= new_time_stamp) then
    if (bmad_status%type_out .and. bmad_status%ok) call out_io(s_warn$, &
                                      r_name, 'NOTE: DIGESTED FILE OUT OF DATE.')
    bmad_status%ok = .false.
  endif

enddo

! we read (and write) the lat in pieces since it is
! too big to write in one piece

if (v86) then
  read (d_unit, err = 9030)  &   
            lat%name, lat%lattice, lat%input_file_name, lat%title, &
            old_a, old_b, old_z, lat%param, lat%version, lat%n_ele_track, &
            lat%n_ele_track, lat%n_ele_max, &
            lat%n_control_max, lat%n_ic_max, lat%input_taylor_order
  lat%a%chrom = old_a%chrom; lat%a%emit = old_a%emit; lat%a%tune = old_a%tune
  lat%b%chrom = old_b%chrom; lat%b%emit = old_b%emit; lat%b%tune = old_b%tune
  lat%z%chrom = old_z%chrom; lat%z%emit = old_z%emit; lat%z%tune = old_z%tune

else
  read (d_unit, err = 9030)  &   
            lat%name, lat%lattice, lat%input_file_name, lat%title, &
            lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track, &
            lat%n_ele_track, lat%n_ele_max, &
            lat%n_control_max, lat%n_ic_max, lat%input_taylor_order
endif

! %control and %ic are allocated to the same length for convenience.

call allocate_lat_ele_array(lat, lat%n_ele_max+10)
allocate (lat%control(lat%n_control_max+10))
allocate (lat%ic(lat%n_control_max+10))

!

do i = 0, lat%n_ele_max
  call read_this_ele(lat%ele(i), i, error)
  if (error) return
enddo

!

do i = 1, lat%n_control_max
  read (d_unit, err = 9040) lat%control(i)
enddo

do i = 1, lat%n_ic_max
  read (d_unit, err = 9050) lat%ic(i)
enddo

read (d_unit, iostat = ios) lat%beam_start

! read branch lines

if (version > 86) then
  read (d_unit, err = 9060) n_branch
  call allocate_branch_array (lat%branch, n_branch, lat)  ! Initial allocation

  do i = 1, n_branch
    branch => lat%branch(i)
    branch%ix_branch = i
    read (d_unit) branch%key, branch%ix_from_branch, branch%ix_from_ele, &
                                      branch%n_ele_track, branch%n_ele_max
    call allocate_ele_array (branch%ele, branch%n_ele_max)
    do j = 0, branch%n_ele_max
      call read_this_ele (branch%ele(j), j, error)
      if (error) return
    enddo
  enddo
endif

! And finish

close (d_unit)

lat%param%stable = .true.  ! Assume this 

return

!------------------

9000  continue
if (bmad_status%type_out) then
   call out_io (s_warn$, r_name, 'DIGESTED FILE DOES NOT EXIST.')
endif
close (d_unit)
bmad_status%ok = .false.
version = -1
return

!--------------------------------------------------------------

9010  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE VERSION.')
endif
close (d_unit)
bmad_status%ok = .false.
return

!--------------------------------------------------------------

9020  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE FILE AND DATE.')
endif
close (d_unit)
bmad_status%ok = .false.
return

!--------------------------------------------------------------

9030  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE LATTICE GLOBALS.')
endif
close (d_unit)
bmad_status%ok = .false.
return

!--------------------------------------------------------------

9040  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE CONTROL.')
endif
close (d_unit)
bmad_status%ok = .false.
return

!--------------------------------------------------------------

9050  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE IC.')
endif
close (d_unit)
bmad_status%ok = .false.
return

!--------------------------------------------------------------

9060  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE N_BRANCH.')
endif
close (d_unit)
bmad_status%ok = .false.
return

!-----------------------------------------------------------------------------------
contains

subroutine read_this_ele (ele, ix_ele, error)

type (ele_struct) ele
integer ix_ele
logical error

!

error = .true.

if (v86) then
  read (d_unit, err = 9100) mode3, ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
          ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
          ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, ele%y, &
          ele%a, ele%b, ele%z, ele%gen0, ele%vec0, ele%mat6, &
          ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
          ele%is_on, ele%sub_key, control_type, ele%ix_value, &
          ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
          ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
          ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
          ele%num_steps, ele%integrator_order, ele%ref_orbit, &
          ele%taylor_order, ele%symplectify, ele%mode_flip, &
          ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
          ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
          ele%coupler_at, ele%on_a_girder, ele%csr_calc_on, &
          ele%map_ref_orb_in, ele%map_ref_orb_out, ele%offset_moves_aperture
 
elseif (v87 .or. v88) then
  read (d_unit, err = 9100) mode3, ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
          ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
          ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, ele%y, &
          ele%a, ele%b, ele%z, ele%gen0, ele%vec0, ele%mat6, &
          ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
          ele%is_on, ele%sub_key, control_type, ele%ix_value, &
          ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
          ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
          ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
          ele%num_steps, ele%integrator_order, ele%ref_orbit, &
          ele%taylor_order, ele%symplectify, ele%mode_flip, &
          ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
          ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
          ele%coupler_at, ele%on_a_girder, ele%csr_calc_on, &
          ele%map_ref_orb_in, ele%map_ref_orb_out, ele%offset_moves_aperture, &
          ele%ix_branch
elseif (v89) then
  read (d_unit, err = 9100) mode3, ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
          ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
          ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, ele%y, &
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
endif

! Decompress value array

read (d_unit, err = 9110) k_max
read (d_unit, err = 9120) ix_value(1:k_max), value(1:k_max)
do k = 1, k_max
  ele%value(ix_value(k)) = value(k)
enddo

!

if (mode3) then
  allocate(ele%mode3)
  read (d_unit, err = 9150) ele%mode3
endif

if (ix_wig /= 0) then
  allocate (ele%wig_term(ix_wig))
  do j = 1, ix_wig
    read (d_unit, err = 9200) ele%wig_term(j)
  enddo
endif

! This is to cover up the change where periodic_type wigglers now have a single wig_term
! where before they did not have any.

if (ele%key == wiggler$ .and. ele%sub_key == periodic_type$ .and. &
                                                       .not. associated(ele%wig_term)) then
  allocate (ele%wig_term(1))
endif

if (ix_const /= 0) then
  allocate (ele%const(ix_const))
  read (d_unit, err = 9300) ele%const
endif

if (any (ix_r /= 0)) then
  allocate (ele%r(ix_r(1):ix_r(3), ix_r(2):ix_r(4)))
  read (d_unit, err = 9400) ele%r
endif

if (ix_d /= 0) then
  allocate (ele%descrip)
  read (d_unit, err = 9500) ele%descrip
endif

if (ix_m /= 0) then
  call multipole_init (ele)
  read (d_unit, err = 9600) ele%a_pole, ele%b_pole
endif
  
do j = 1, 6
  if (ix_t(j) == 0) cycle
  read (d_unit) ele%taylor(j)%ref
  allocate (ele%taylor(j)%term(ix_t(j)))
  do k = 1, ix_t(j)
    read (d_unit, err = 9700) ele%taylor(j)%term(k)
  enddo
enddo

! If ix_lr is negative then it is a pointer to a previously read wake. 
! See write_digested_bmad_file.

if (ix_sr_table /= 0 .or. ix_sr_mode_long /= 0 .or. ix_sr_mode_trans /= 0 .or. ix_lr /= 0) then
  if (ix_lr < 0) then
    call transfer_wake (lat%ele(abs(ix_lr))%wake, ele%wake)

  else
    call init_wake (ele%wake, ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr)
    read (d_unit, err = 9800) ele%wake%sr_file
    read (d_unit, err = 9810) ele%wake%sr_table
    read (d_unit, err = 9840) ele%wake%sr_mode_long
    read (d_unit, err = 9850) ele%wake%sr_mode_trans
    read (d_unit, err = 9820) ele%wake%lr_file
    read (d_unit, err = 9830) ele%wake%lr
    read (d_unit, err = 9860) ele%wake%z_sr_mode_max
  endif
endif

ele%value(check_sum$) = 0
if (associated(ele%a_pole)) ele%value(check_sum$) = sum(ele%a_pole) + sum(ele%b_pole)
 
ele%old_value = ele%value

error = .false.

return

!--------------------------------------------------------------

9100  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING ELEMENT # \i0\ ', &
                                  i_array = (/ ix_ele /) )
endif
close (d_unit)
bmad_status%ok = .false.
return

9110  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING K_MAX OF ELEMENT # \i0\ ', &
                                  i_array = (/ ix_ele /) )
endif
close (d_unit)
bmad_status%ok = .false.
return

9120  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
                                 'ERROR READING VALUES OF ELEMENT # \i0\ ', &
                                  i_array = (/ ix_ele /) )
endif
close (d_unit)
bmad_status%ok = .false.
return

9150  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING MODE3 COMPONENT FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9200  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
        'ERROR READING WIGGLER TERM FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9300  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING IX_CONST TERM FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9400  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING R TERM FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9500  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING DESCRIP TERM FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9600  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING AN,BN TERMS FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9700  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING TAYLOR TERMS FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9800  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%SR_FILE FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9810  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_table FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9820  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%LR_FILE FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9830  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%LR FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9840  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_mode_long FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9850  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_mode_trans FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

9860  continue
if (bmad_status%type_out) then
   call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%Z_CUT_SR FOR ELEMENT: ' // ele%name)
endif
close (d_unit)
bmad_status%ok = .false.
return

end subroutine

end subroutine
