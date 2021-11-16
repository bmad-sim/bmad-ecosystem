!+
! Subroutine init_bmad ()
!
! Routine to initalize bmad related stuff.
!-

subroutine init_bmad()

use attribute_mod, dummy => init_bmad
use ptc_interface_mod, dummy2 => init_bmad
implicit none

!

call set_ptc_com_pointers()
call init_attribute_name_array()

end subroutine
