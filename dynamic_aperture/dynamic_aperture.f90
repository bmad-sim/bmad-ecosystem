!+
! Program dynamic_aperture
!
! Program for maping out the dynamic aperture.
! Not to be confused with the CESR centric dynamic_aperture_cesr program.
!-

program dynamic_aperture_program

use da_program_mod

implicit none

type (lat_struct), target :: lat
type (aperture_param_struct) :: da_param
type (aperture_point_struct), pointer :: da_point
type (aperture_scan_struct), allocatable, target :: aperture_scan(:)
type (aperture_scan_struct), pointer :: da
type (coord_struct), allocatable :: closed_orb(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele

real(rp) dpz(20)
real(rp) :: ramping_start_time = 0
integer nargs, ios, i, j, n_dpz

logical :: ramping_on = .false.
logical set_rf_off, err

character(160) in_file, lat_file, dat_file, gnu_command

namelist / params / lat_file, bmad_com, set_rf_off, da_param, dpz, dat_file, &
            ramping_start_time, ramping_on

! Get parameters

in_file = 'dynamic_aperture.init'

nargs = command_argument_count()
if (nargs == 1) then
  call get_command_argument(1, in_file)
elseif (nargs > 1) then
  print *, 'ERROR: EXTRA STUFF ON COMMAND LINE'
  stop
endif

print *, 'Parameter input file: ', trim(in_file)

open (unit = 1, file = in_file, status = 'old', iostat = ios)

if (ios /= 0) then
  print *, 'ERROR: CANNOT OPEN FILE: ', trim(in_file)
  stop
endif

set_rf_off = .false.
bmad_com%auto_bookkeeper = .false.   ! Makes tracking faster
bmad_com%absolute_time_tracking_default = .true.
dpz = real_garbage$

read (1, nml = params)

close (1)

da_param%min_angle = da_param%min_angle * pi / 180
da_param%max_angle = da_param%max_angle * pi / 180

da_com%ramping_on = ramping_on
da_com%ramping_start_time = ramping_start_time

! Read in lattice

call bmad_parser (lat_file, lat)

call twiss_and_track(lat, closed_orb)

n_dpz = count(dpz /= real_garbage$)
print *, 'Note: Number of dpz points: ', n_dpz

if (set_rf_off) then
  print *, 'Note: RF being turned off for tracking.'
  call set_on_off (rfcavity$, lat, off$)
endif

if (.not. lat%absolute_time_tracking) then
  print *, 'Note: absolute time tracking is OFF!'
endif

print *, 'Data file: ', trim(dat_file)

write (gnu_command, '(a, i0, 3a)') 'gnuplot plotting command: plot for [IDX=1:', &
                  n_dpz, '] "', trim(dat_file), '" index (IDX-1) u 1:2 w lines title columnheader(1)'

! Ramper setup

if (ramping_on) then
  call lat_ele_locator ('RAMPER::*', lat, eles, da_com%n_ramper_loc, err)
  if (da_com%n_ramper_loc == 0) then
    print '(2a)', 'Warning! NO RAMPER ELEMENTS FOUND IN LATTICE.'
    stop
  endif

  do i = 1, da_com%n_ramper_loc
    ele => eles(i)%ele
    if (ele%control%var(1)%name /= 'TIME') then
      print *, 'Note! This ramper does not use "time" as the control variable: ' // trim(ele%name)
      print *, '      This ramper will not be directly varied in the simulation.'
    endif
    da_com%ramper(i) = lat_ele_loc_struct(ele%ix_branch, ele%ix_ele)
  enddo
endif

! Scan

open (1, file = dat_file)
write (1, '(2a)')        '# lat_file              = ', trim(lat_file)
write (1, '(a, l1)')     '# set_rf_off            = ', set_rf_off
write (1, '(a, f10.1)')  '# da_param%min_angle    =', da_param%min_angle
write (1, '(a, f10.1)')  '# da_param%max_angle    =', da_param%max_angle
write (1, '(a, es10.2)') '# da_param%rel_accuracy =', da_param%rel_accuracy
write (1, '(a, es10.2)') '# da_param%abs_accuracy =', da_param%abs_accuracy
write (1, '(a, f10.5)')  '# da_param%x_init       =', da_param%x_init
write (1, '(a, f10.5)')  '# da_param%y_init       =', da_param%y_init
write (1, '(a, i0)')     '# da_param%n_turn       = ', da_param%n_turn
write (1, '(a, i0)')     '# da_param%n_angle      = ', da_param%n_angle
write (1, '(2a)')        '# ', trim(gnu_command)

call dynamic_aperture_scan (aperture_scan, da_param, dpz(1:n_dpz), lat)

do i = 1, n_dpz
  da => aperture_scan(i)

  write (1, *)
  write (1, *)
  write (1, '(a, f10.6, a)') '"Dpz =', dpz(i), '"'
  write (1, '(a, f10.6, a)') '"x_ref_orb =', da%ref_orb%vec(1), '"'
  write (1, '(a, f10.6, a)') '"y_ref_orb =', da%ref_orb%vec(1), '"'
  do j = 1, da_param%n_angle
    da_point => da%point(j)
    write (1, '(2f11.6, i7, 6x, a, 3x, a)') da_point%x, da_point%y, da_point%i_turn, &
                              coord_state_name(da_point%plane), trim(lat%ele(da_point%ix_ele)%name)
  enddo
enddo

close (1)

print '(a)', trim(gnu_command)

end program
