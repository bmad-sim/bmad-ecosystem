!+
! Program dynamic_aperture
!
! Program for maping out the dynamic aperture.
! Not to be confused with the CESR centric dynamic_aperture_cesr program.
!-

program dynamic_aperture_program

use bmad
use dynamic_aperture_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orbit(:)
type (coord_struct) orb0
type (branch_struct), pointer :: branch
type (aperture_param_struct) :: da_param
type (aperture_data_struct), pointer :: ddat
type (aperture_scan_struct), target :: da

real(rp) dpz(20)
integer nargs, ios, i, j, n_dpz

logical set_rf_off

character(160) in_file, lat_file, dat_file, gnu_command

namelist / params / lat_file, bmad_com, set_rf_off, da_param, dpz, dat_file

! Get parameters

in_file = 'dynamic_aperture.init'

nargs = cesr_iargc()
if (nargs == 1) then
  call cesr_getarg(1, in_file)
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

da%param = da_param
da%param%min_angle = da%param%min_angle * pi / 180
da%param%max_angle = da%param%max_angle * pi / 180

! Read in lattice

call bmad_parser (lat_file, lat)
branch => lat%branch(0)

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

! Scan

open (1, file = dat_file)
write (1, '(2a)')        '# lat_file           = ', trim(lat_file)
write (1, '(a, l1)')     '# set_rf_off         = ', set_rf_off
write (1, '(a, f10.1)')  '# da_param%min_angle =', da_param%min_angle
write (1, '(a, f10.1)')  '# da_param%max_angle =', da_param%max_angle
write (1, '(a, es10.2)') '# da_param%accuracy  =', da_param%accuracy
write (1, '(a, f10.5)')  '# da_param%x_init    =', da_param%x_init
write (1, '(a, f10.5)')  '# da_param%y_init    =', da_param%y_init
write (1, '(a, i0)')     '# da_param%n_turn    = ', da_param%n_turn
write (1, '(a, i0)')     '# da_param%n_angle   = ', da_param%n_angle
write (1, '(2a)')        '# ', trim(gnu_command)

call reallocate_coord(orbit, lat, branch%ix_branch)
if (rf_is_on(lat%branch(0))) then
  call closed_orbit_calc (lat, orbit, 6)
  orb0 = orbit(0)
endif

do i = 1, size(dpz)
  if (dpz(i) == real_garbage$) cycle
  print *, 'Tracking at dpz:', dpz(i)

  if (rf_is_on(lat%branch(0))) then
    orbit = orb0
    orbit%vec(6) = orbit%vec(6) + dpz(i)
  else
    orbit(0)%vec(6) = dpz(i)
    call twiss_and_track (lat, orbit)
!    call closed_orbit_calc (lat, orbit, 4)
  endif

  da%ref_orb = orbit(0)
  call dynamic_aperture_scan (lat%branch(0), da)

  write (1, *)
  write (1, *)
  write (1, '(a, f10.6)') '"Dpz =', dpz(i), '"'
  do j = 1, da_param%n_angle
    ddat => da%aperture(j)
    write (1, '(2f11.6, i7, 6x, a, 3x, a)') ddat%x, ddat%y, ddat%i_turn, &
                              coord_state_name(ddat%plane), trim(lat%ele(ddat%ix_ele)%name)
  enddo

enddo

close (1)

print '(a)', trim(gnu_command)

end program
