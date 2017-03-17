program aperture_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb

integer :: i, j, ib, nargs, isn

logical print_extra
 
character(100) lat_file

!

global_com%exit_on_error = .false.
lat_file = 'aperture_test.bmad'

print_extra = .false.
nargs = cesr_iargc()
if (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit

elseif (nargs > 0)then
  call cesr_getarg(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
endif

open (1, file = 'output.now')

!

call bmad_parser (lat_file, lat)
ele => lat%ele(2)
start_orb = lat%beam_start

ele%tracking_method = runge_kutta$
call track1(start_orb, lat%ele(2), lat%param, end_orb)
write (1, '(a, 6f10.6, i4)') '"RK"   ABS 0', end_orb%vec, end_orb%state

ele%tracking_method = time_runge_kutta$
call track1(start_orb, lat%ele(2), lat%param, end_orb)
write (1, '(a, 6f10.6, i4)') '"TRK"  ABS 0', end_orb%vec, end_orb%state


!

close (1)

end program
