!+
! Program ptc_layout_example
!
! Program to show how to setup a PTC layout.
!-

program ptc_layout_example

use ptc_layout_mod

implicit none

type (lat_struct) lat
type (fibre), pointer :: ptc_fibre

!

call bmad_parser('cell1.bmad', lat)
call lat_to_ptc_layout (lat)

ptc_fibre => lat%ele(5)%ptc_fibre
call type_ptc_fibre (ptc_fibre)

end program
