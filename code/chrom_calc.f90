!+
! Subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, exit_on_error)
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
!   exit_on_error -- Logical, optional: If True then subroutine will terminate 
!                         programif the orbit does not converge. Default is
!                         determined by bmad_status%exit_on_error
!   bmad_status   -- Bmad status common block
!     %exit_on_error -- Default for exit_on_error argument.
!
! Output:
!   delta_e     -- Real(rp): Set to 1.0e-4 if on input DELTA_E =< 0.
!   chrom_x     -- Real(rp): Horizontal chromaticity.
!   chrom_y     -- Real(rp): Vertical chromaticity.
!   bmad_status -- Status structure in common block
!                  See MAT_SYMP_DECOUPLE for for more info.
!       %ok          -- Set False if orbit does not converge, True otherwise.
!       %type_out    -- Logical: If .true. then will type a message for
!                       non ok$ bmad_status
!-

#include "CESR_platform.inc"

subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, exit_on_error)

  use bookkeeper_mod, except => chrom_calc

  implicit none

  type (lat_struct)  lat
  type (lat_struct), save :: lat2
  type (coord_struct), allocatable, save :: orb(:)

  real(rp) high_tune_x, high_tune_y, low_tune_x, low_tune_y
  real(rp) delta_e, chrom_x, chrom_y
  integer nt
  logical, optional :: exit_on_error

!

  call reallocate_coord (orb, lat%n_ele_max)

  if (delta_e <= 0) delta_e = 1.0e-4

  lat2 = lat
  call set_on_off (rfcavity$, lat2, off$)

  nt = lat%n_ele_track

! lower energy tune

  orb(0)%vec(6) = -delta_e
  call closed_orbit_calc (lat2, orb, 4, exit_on_error = exit_on_error)
  if (.not. bmad_status%ok) return
  call lat_make_mat6 (lat2, -1, orb)

  call twiss_at_start (lat2)
  if (.not. bmad_status%ok) return
  call twiss_propagate_all (lat2)
  low_tune_x = lat2%ele(nt)%a%phi / twopi
  low_tune_y = lat2%ele(nt)%b%phi / twopi

! higher energy tune

  orb(0)%vec(6) = delta_e
  call closed_orbit_calc (lat2, orb, 4, exit_on_error = exit_on_error)
  if (.not. bmad_status%ok) return
  call lat_make_mat6 (lat2, -1, orb)

  call twiss_at_start (lat2)
  if (.not. bmad_status%ok) return
  call twiss_propagate_all (lat2)
  high_tune_x = lat2%ele(nt)%a%phi / twopi
  high_tune_y = lat2%ele(nt)%b%phi / twopi

! compute the chrom

  chrom_x = (high_tune_x - low_tune_x) / (2 * delta_e)
  chrom_y = (high_tune_y - low_tune_y) / (2 * delta_e)

end subroutine
