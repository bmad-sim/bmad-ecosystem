!+
! Subroutine type_coord (coord)
!
! Subroutine to type out a coordinate.
!
! Input:
!   coord -- Coord_struct: Coordinate
!-

subroutine type_coord (coord)

use bmad_interface, except_dummy => type_coord

implicit none
type (coord_struct)  coord

character(16), parameter :: r_name = 'type_coord'

!

call out_io (s_blank$, r_name, '(X, X''): \2es16.6\ ', '(Y, Y''): \2es16.6\ ', '(Z, Z''): \2es16.6\ ', &
                        r_array = coord%vec)

end subroutine

