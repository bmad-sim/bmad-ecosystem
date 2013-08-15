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

do i = 1, lat%n_ele_track
  ele => lat%ele(i)

  if (ele%name == 'Q1') then
    write (1, '(a, f10.4)') '"Q1[K1]"     ABS 0', ele%value(k1$) 
    write (1, '(a, f10.4)') '"Q1[TILT]"   ABS 0', ele%value(tilt$) 
    write (1, '(a, f10.4)') '"Q1[HKICK]"  ABS 0', ele%value(hkick$) 
  endif

  write (1, '(3a, f10.4)') '"', trim(ele%name), '[L]"      ABS 0', ele%value(l$)
enddo

do i = lat%n_ele_track+1, lat%n_ele_max

  ele => lat%ele(i)
  if (ele%lord_status == group_lord$) then
    write (1, '(3a, f10.4)') '"', trim(ele%name), '[COMMAND]"      ABS 0', ele%value(command$)
  endif

  if (ele%lord_status == super_lord$) then
    write (1, '(3a, f10.4)') '"', trim(ele%name), '[L]"      ABS 0', ele%value(l$)
  endif

  if (ele%lord_status == overlay_lord$) then
    write (1, '(5a, f10.4)') '"', trim(ele%name), '[', &
                      trim(ele%component_name), ']"      ABS 0', ele%value(ele%ix_value)
  endif

enddo

close(1)

end program
