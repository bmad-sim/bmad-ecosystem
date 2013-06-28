program control_test

use bmad
use mad_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

character(40) :: lat_file  = 'control_test.bmad'

integer :: i, j, k, nargs
logical print_extra

!

print_extra = .false.
nargs = cesr_iargc()
if (nargs == 1)then
   call cesr_getarg(1, lat_file)
   print *, 'Using ', trim(lat_file)
   print_extra = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

call bmad_parser (lat_file, lat, make_mats6 = .false.)

!

open (1, file = 'output.now', recl = 200)

do i = 1, lat%n_ele_max
  ele => lat%ele(i)

  if (ele%name == 'Q1') then
    write (1, '(a, f10.4)') '"Q1:L"      ABS 0', ele%value(l$)
    write (1, '(a, f10.4)') '"Q1:K1"     ABS 0', ele%value(k1$) 
    write (1, '(a, f10.4)') '"Q1:TILT"   ABS 0', ele%value(tilt$) 
    write (1, '(a, f10.4)') '"Q1:HKICK"  ABS 0', ele%value(hkick$) 
  endif

  if (ele%name == 'Q2') then
    write (1, '(a, f10.4)') '"Q2:L"      ABS 0', ele%value(l$)
  endif

enddo

close(1)

end program
