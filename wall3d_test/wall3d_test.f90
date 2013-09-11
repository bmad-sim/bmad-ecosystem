program wall3d_test

use wall3d_mod

implicit none

type (lat_struct), target :: lat
type (wall3d_struct), pointer :: wall

integer i

!

call bmad_parser ('wall3d_test.bmad', lat)
wall => lat%branch(0)%wall3d

open (1, file = 'output.now')


write (1, '(a, i0)') '"n_wall"  ABS 0  ', size(wall%section)

do i = 1, size(wall%section)
  write (1, '(a, i3, a, i3, f10.3)') '"Section', i, '" ABS 0  ', &
                                wall%section(i)%ix_ele, wall%section(i)%s
enddo

end program
