!+
! Subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
!                        low_E_lat, high_E_lat, low_E_orb, high_E_orb)
!
! Subroutine to calculate the chromaticities by computing the tune change
! when then energy is changed.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat          -- lat_struct: Lat
!   delta_e       -- Real(rp): Delta energy used for the calculation.
!                      If 0 then default of 1.0e-4 is used.
!
! Output:
!   delta_e      -- Real(rp): Set to 1.0e-4 if on input DELTA_E =< 0.
!   chrom_x      -- Real(rp): Horizontal chromaticity.
!   chrom_y      -- Real(rp): Vertical chromaticity.
!   err_flag     -- Logical, optional: Set true if there is an error. False otherwise.
!   low_E_lat    -- lat_struct, optional: Lattice with RF off and matrices computed at E_lat - delta_e
!   high_E_lat   -- lat_struct, optional: Lattice with RF off and matrices computed at E_lat + delta_e
!   low_E_orb(:) -- coord_struct, allocatable, optional: Orbit computed at E_lat - delta_e.
!   high_E_orb(:)-- coord_struct, allocatable, optional: Orbit computed at E_lat + delta_e.
!-

subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
                       low_E_lat, high_E_lat, low_E_orb, high_E_orb)

use bookkeeper_mod, except_dummy => chrom_calc

implicit none

type (lat_struct)  lat
type (lat_struct), optional, target :: low_E_lat, high_E_lat
type (coord_struct), allocatable, optional, target :: low_E_orb(:), high_E_orb(:)
type (lat_struct), target :: this_lat
type (lat_struct), pointer :: lat2
type (coord_struct), allocatable, target :: this_orb(:)

real(rp) high_tune_x, high_tune_y, low_tune_x, low_tune_y
real(rp) delta_e, chrom_x, chrom_y

integer nt, stat

logical, optional, intent(out) :: err_flag
logical err, used_this_lat

! Init setup

if (present(err_flag)) err_flag = .true.
if (delta_e <= 0) delta_e = 1.0e-4

nt = lat%n_ele_track

! lower energy tune

if (present(low_E_lat)) then
  lat2 => low_E_lat
  used_this_lat = .false.
else
  lat2 => this_lat
  used_this_lat = .true.
endif

lat2 = lat
call set_on_off (rfcavity$, lat2, off$)

if (present(low_E_orb)) then
  call reallocate_coord (low_E_orb, lat%n_ele_max)
  low_E_orb(0)%vec(6) = -delta_e
  call closed_orbit_calc (lat2, low_E_orb, 4, err_flag = err)
  if (err) return
  call lat_make_mat6 (lat2, -1, low_E_orb)
else
  call reallocate_coord (this_orb, lat%n_ele_max)
  this_orb(0)%vec(6) = -delta_e
  call closed_orbit_calc (lat2, this_orb, 4, err_flag = err)
  if (err) return
  call lat_make_mat6 (lat2, -1, this_orb)
endif

call twiss_at_start (lat2, stat)
if (stat /= ok$) return
call twiss_propagate_all (lat2)
low_tune_x = lat2%ele(nt)%a%phi / twopi
low_tune_y = lat2%ele(nt)%b%phi / twopi

! higher energy tune

if (present(high_E_lat)) then
  lat2 => high_E_lat
  used_this_lat = .false.
else
  lat2 => this_lat
endif

if (.not. used_this_lat) then
  call transfer_lat (lat, lat2)
  call set_on_off (rfcavity$, lat2, off$)
endif

if (present(high_E_orb)) then
  call reallocate_coord (high_E_orb, lat%n_ele_max)
  high_E_orb(0)%vec(6) = delta_e
  call closed_orbit_calc (lat2, high_E_orb, 4, err_flag = err)
  if (err) return
  call lat_make_mat6 (lat2, -1, high_E_orb)
else
  call reallocate_coord (this_orb, lat%n_ele_max)
  this_orb(0)%vec(6) = delta_e
  call closed_orbit_calc (lat2, this_orb, 4, err_flag = err)
  if (err) return
  call lat_make_mat6 (lat2, -1, this_orb)
endif

call twiss_at_start (lat2, stat)
if (stat /= ok$) return
call twiss_propagate_all (lat2)
high_tune_x = lat2%ele(nt)%a%phi / twopi
high_tune_y = lat2%ele(nt)%b%phi / twopi

! compute the chrom

chrom_x = (high_tune_x - low_tune_x) / (2 * delta_e)
chrom_y = (high_tune_y - low_tune_y) / (2 * delta_e)

if (present(err_flag)) err_flag = .false.

end subroutine
