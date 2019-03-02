!+
! Program patch_test
!
! This program is part of the Bmad regression testing suite.
!-

program patch_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave, slave2
type (coord_struct) :: start_orb, end_orb

integer ip
character(40) fmt

!

call bmad_parser ('patch_test.bmad', lat)

!

open (1, file = 'output.now')
fmt = '(3a, 3f20.15, 5x, 3f20.15)'

do ip = 1, 3
  ele => lat%ele(ip)
  if (ele%key == marker$) cycle
!  if (ele%key /= patch$) cycle
  call init_coord (start_orb, lat%particle_start, ele, upstream_end$)
  call track1 (start_orb, ele, lat%param, end_orb)
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '" ABS 0', end_orb%vec
  if (ele%key == patch$) then
    write (1, '(a, f20.14)') '"L" REL 1E-12 ', ele%value(l$)
  endif
enddo

ele => lat%ele(4)
write (1, '(a, 6es10.2)') '"Flexible" REL 1E-15 ', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi

!

close (1)

end program
