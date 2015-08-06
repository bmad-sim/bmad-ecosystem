program bookkeeper_test

use bmad
use mad_mod
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat, lat2
type (ele_struct), pointer :: ele, nele

character(40) :: lat_file  = 'bookkeeper_test.bmad'
character(100) str

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

!-------------

call bmad_parser ('bookkeeper_test2.bmad', lat)

do i = 1, lat%n_ele_max
  lat%ele(i)%select = .false.
  if (lat%ele(i)%type == 'A') lat%ele(i)%select = .true.
enddo

call make_hybrid_lat (lat, lat2)
call write_bmad_lattice_file('out.bmad', lat2)

!-------------

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
  if (ele%lord_status == super_lord$) then
    write (1, '(3a, f10.4)') '"', trim(ele%name), '[L]"      ABS 0', ele%value(l$)
  endif

  if (ele%key == overlay$ .or. ele%key == group$) then
    do j = 1, size(ele%control_var)
      write (1, '(5a, f10.4)') '"', trim(ele%name), '[', &
                      trim(ele%control_var(j)%name), ']"      ABS 0', ele%control_var(j)%value
    enddo
  endif

  if (ele%name == 'GRN') then
    j = ele%ix1_slave
    str = expression_stack_to_string(lat%control(j)%stack)
    write (1, '(5a, f10.4)') '"GRN[string]" STR "', trim(str), '"'
  endif

enddo

! pointer_to_next_ele test

nele => pointer_to_next_ele(lat%ele(1), 7)
write (1, '(a, i4)') '"Next-1" ABS 0', nele%ix_ele

nele => pointer_to_next_ele(lat%ele(1), 7, .true.)
write (1, '(a, i4)') '"Next-2" ABS 0', nele%ix_ele

nele => pointer_to_next_ele(lat%ele(1), -7)
write (1, '(a, i4)') '"Next-3" ABS 0', nele%ix_ele

nele => pointer_to_next_ele(lat%ele(1), -7, .true.)
write (1, '(a, i4)') '"Next-4" ABS 0', nele%ix_ele


close(1)

end program
