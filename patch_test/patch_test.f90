!+
! Program girder_test
!
! This program is part of the Bmad regression testing suite.
!-

program patch_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: patch, slave, slave2
type (coord_struct) :: start_orb, end_orb

integer ip
character(40) fmt

!

call bmad_parser ('patch_test.bmad', lat)

!

open (1, file = 'output.now')
fmt = '(3a, 3f20.15, 5x, 3f20.15)'

do ip = 1, lat%n_ele_max
  patch => lat%ele(ip)
  if (patch%key == marker$) cycle
!  if (patch%key /= patch$) cycle
  call init_coord (start_orb, lat%beam_start, patch, .false.)
  call track1 (start_orb, patch, lat%param, end_orb)
  write (1, '(3a, 6es14.6)') '"', trim(patch%name), '" ABS 0', end_orb%vec

enddo

!

close (1)

end program
