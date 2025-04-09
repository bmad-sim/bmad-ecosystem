!+
! Subroutine init_custom (lat)
!
! Dummy routine for initializing custom elements or elements that need custom calculations. 
! Custom calculations are done if any one of the following ele_struct components is set to custom$:
!   ele%tracking_method
!   ele%mat6_calc_method
!   ele%field_calc
!   ele%aperture_type
!
! Also procedure pointers can be setup.
!
! If called, this routine will do nothing.
! This routine needs to be replaced for custom initialization.
!
! Note!! The linker for MacOS does not allow this dummy init_custom to be overridden. 
! And there is no guarantee that overriding will work on any other platform.
!
! Input:
!   lat     -- lat_struct: Lattice to be customized.
!
! Output:
!   lat     -- lat_struct: 
!+

subroutine init_custom (lat)

use bmad_struct
use bmad_interface, except_dummy => init_custom

implicit none

type (lat_struct), target :: lat

!

end subroutine
