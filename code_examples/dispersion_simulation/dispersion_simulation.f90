!+
! Program dispersion_simulation
!
! Program to simulate a dispersion measurement where the RF frequency is varied and the change
! in orbit measured.
!-

program dispersion_simulation

use bmad

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orb0(:), orb1(:), orb2(:)
type (branch_struct), pointer :: branch

real(rp) rf_freq0, dfreq, pz1, pz2, dpz
real(rp), allocatable :: on_off_vals(:)

integer n, j, ix, n_ele_track, ix_branch
integer, allocatable :: ix_rf(:)

character(100) lat_file

logical err_flag, radiation_damping_on

namelist / params / ix_branch, lat_file, dfreq, radiation_damping_on

! Get lattice and calculate Twiss parameters

open (1, file = 'dispersion_simulation.init', status = 'old')
read (1, nml = params)
close (1)

call bmad_parser (lat_file, lat)
branch => lat%branch(ix_branch)       
n_ele_track = branch%n_ele_track

call set_on_off (rfcavity$, lat, off_and_save$, ix_branch = ix_branch, saved_values = on_off_vals)
call twiss_and_track (lat, orb0)
call set_on_off (rfcavity$, lat, restore_state$, ix_branch = ix_branch, saved_values = on_off_vals)

! Make a list of where the RF cavities are.

allocate (ix_rf(count(branch%ele%key == rfcavity$)))

j = 0
do n = 1, branch%n_ele_max
  if (branch%ele(n)%key /= rfcavity$) cycle
  j = j + 1
  ix_rf(j) = n
enddo

rf_freq0 = branch%ele(ix_rf(1))%value(rf_frequency$)  ! Assume all cavities have the same frequency

! Switch to absolute time tracking since relative time tracking will not give correct results.
! See the Bmad manual discussion on relative verses absolute time tracking.
! With absolute time tracking, use auto phasing to set phi0_ref.

bmad_com%absolute_time_tracking = .true.

do n = 1, size(ix_rf)
  ix = ix_rf(n)
  call autoscale_phase_and_amp (branch%ele(n), branch%param, err_flag, scale_amp = .false., scale_phase = .true.)
enddo

! Turn on radiation dampling for extra realism

bmad_com%radiation_damping_on = radiation_damping_on

! Closed orbit at lower energy (higher frequency)

call reallocate_coord(orb1, lat, ix_branch)

do n = 1, size(ix_rf)
  ix = ix_rf(n)
  branch%ele(ix)%value(rf_frequency$) = rf_freq0 + dfreq / 2
enddo

call closed_orbit_calc (lat, orb1, 6, 1, ix_branch)

pz1 = sum(orb1(1:n_ele_track)%vec(6)) / n_ele_track

! Closed orbit at higher energy (lower frequency)

call reallocate_coord(orb2, lat, ix_branch)

do n = 1, size(ix_rf)
  ix = ix_rf(n)
  branch%ele(ix)%value(rf_frequency$) = rf_freq0 - dfreq / 2
enddo

call closed_orbit_calc (lat, orb2, 6, 1, ix_branch)

pz2 = sum(orb2(1:n_ele_track)%vec(6)) / n_ele_track

! Compute some numbers

dpz = pz2 - pz1

print *, '%Change in energy:            ', dpz
print *, 'Average horizontal orbit diff:', sum(orb2(1:n_ele_track)%vec(1) - orb1(1:n_ele_track)%vec(1)) / n_ele_track

! Write results

open (1, file = 'dispersion_simulation.dat')

do n = 1, n_ele_track
  write (1, '(i4, f10.4, 3x, 2f11.6, 3x, f10.4, 3x, 2f11.6)') n, branch%ele(n)%s, &
              (orb2(n)%vec(1) - orb1(n)%vec(1)) / dpz, branch%ele(n)%x%eta, &
              ((orb2(n)%vec(6) - orb1(n)%vec(6)) - dpz) / dpz, &
              (orb2(n)%vec(3) - orb1(n)%vec(3)) / dpz, branch%ele(n)%y%eta
              
enddo

close (1)
print *, 'Data file: dispersion_simulation.dat'

end program
