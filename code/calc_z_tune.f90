!+
! Subroutine calc_z_tune (lat)
!
! Routine to calculate the synchrotron tune from the full 6X6 1-turn matrix.
!
! Note: The tune will be negative above transition which corresponds to
! counter-clockwise rotation in z-pz space.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat  -- lat_struct: Lat
!
! Output:
!   lat -- lat_struct
!     %z%tune            -- Synchrotron tune (radians)
!     %param%t1_with_RF  -- 6x6 1-turn matrix.
!-

subroutine calc_z_tune (lat)

  use bmad_struct
  use bmad_interface, except_dummy => calc_z_tune
  use nrtype
  use nr

  implicit none

  type (lat_struct) lat

  real(rp) a(6,6), wr(6), wi(6), cos_z

  integer i
!

  call transfer_matrix_calc (lat, .true., a)
  lat%param%t1_with_RF = a

  cos_z = (a(5,5) + a(6,6)) / (2 * (a(5,5)*a(6,6) - a(5,6)*a(6,5)))

  if (cos_z > 0.9999999) then
    lat%z%tune = 0
    return
  endif

  call balanc(a)
  call elmhes(a)
  call hqr(a,wr,wi)

! we need to find which eigen-value is closest to the z_tune

  i = minloc(abs(wr-cos_z), 1)
  lat%z%tune = -abs(atan2(wi(i),wr(i)))

end subroutine
