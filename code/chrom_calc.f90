!+
! Subroutine chrom_calc (ring, delta_e, chrom_x, chrom_y)
!
! Subroutine to calculate the chromaticities by computing the tune change
! when then energy is changed.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring      -- Ring_struct: Ring
!   delta_e   -- Real: Delta energy used for the calculation.
!                    If 0 then default of 1.0e-4 is used.
!
! Output:
!     delta_e     -- Real: Set to 1.0e-4 if on input DELTA_E =< 0.
!     chrom_x     -- Real: Horizontal chromaticity.
!     chrom_y     -- Real: Vertical chromaticity.
!     bmad_status -- Status structure in common block
!                    See MAT_SYMP_DECOUPLE for for more info.
!           %type_out  -- Logical: If .true. then will type a message for
!                         non ok$ bmad_status
!-


subroutine chrom_calc (ring, delta_e, chrom_x, chrom_y)

  use bmad_struct
  implicit none

  type (ring_struct)  ring, ring2
  type (coord_struct)  c0, coord_(0:n_ele_maxx)

  integer i, key

  real high_tune_x, high_tune_y, low_tune_x, low_tune_y
  real delta_e, chrom_x, chrom_y

  if (delta_e <= 0) delta_e = 1.0e-4
  ring2 = ring

! lower energy tune

  coord_(0)%z%vel = -delta_e
  call closed_orbit_at_start (ring2, coord_(0), 4, .true.)
  call track_all (ring2, coord_)
  call ring_make_mat6 (ring2, -1, coord_)

  call twiss_at_start (ring2)
  low_tune_x = ring2%x%tune / twopi
  if (low_tune_x < 0) low_tune_x = 1 + low_tune_x
  low_tune_y = ring2%y%tune / twopi
  if (low_tune_y < 0) low_tune_y = 1 + low_tune_y

! higher energy tune

  coord_(0)%z%vel = delta_e
  call closed_orbit_at_start (ring2, coord_(0), 4, .true.)
  call track_all (ring2, coord_)
  call ring_make_mat6 (ring2, -1, coord_)

  call twiss_at_start (ring2)
  high_tune_x = ring2%x%tune / twopi
  if (high_tune_x < 0) high_tune_x = 1 + high_tune_x
  high_tune_y = ring2%y%tune / twopi
  if (high_tune_y < 0) high_tune_y = 1 + high_tune_y

! compute the chrom

  chrom_x = (high_tune_x - low_tune_x) / (2 * delta_e)
  chrom_y = (high_tune_y - low_tune_y) / (2 * delta_e)

end subroutine
