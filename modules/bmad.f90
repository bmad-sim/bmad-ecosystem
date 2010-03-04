!+
! Module bmad
!
! A "use bmad" statement will define the most common bmad structures needed 
! by a program.
!-

module bmad

  use runge_kutta_mod
  use track1_mod
  use boris_mod
  use reverse_mod
  use dynamic_aperture_mod
  use bookkeeper_mod
  use radiation_mod
  use symp_lie_mod
  use lat_ele_loc_mod
  use twiss_and_track_mod

  ! This is to suppress the ranlib "has no symbols" message
  integer, private :: private_dummy

end module
