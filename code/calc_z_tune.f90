!+
! Subroutine calc_z_tune (ring)
!
! Subroutine to calculate the synchrotron tune from the full 6X6 1-turn matrix
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring  -- Ring_struct: Ring
!
! Output:
!   ring -- Ring_struct
!     %z%tune            -- Synchrotron tune (radians)
!     %param%t1_with_RF  -- 6x6 1-turn matrix.
!-

#include "CESR_platform.inc"

subroutine calc_z_tune ( ring)

  use bmad_struct
  use bmad_interface, except => calc_z_tune
  use nrtype
  use nr

  implicit none

  type (ring_struct) ring

  real(rp) a(6,6), wr(6), wi(6), cos_z

  integer i
!

  call one_turn_matrix (ring, .true., a)
  ring%param%t1_with_RF = a

  cos_z = (a(5,5) + a(6,6)) / (2 * (a(5,5)*a(6,6) - a(5,6)*a(6,5)))

  if (cos_z > 0.9999999) then
    ring%z%tune = 0
    return
  endif

  call balanc(a)
  call elmhes(a)
  call hqr(a,wr,wi)

! we need to find which eigen-value is closest to the z_tune

  i = minloc(abs(wr-cos_z), 1)
  ring%z%tune = -abs(atan2(wi(i),wr(i)))

end subroutine
