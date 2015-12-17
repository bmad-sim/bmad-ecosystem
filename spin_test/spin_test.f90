program spin_test

use bmad
use ptc_spin_mod

implicit none

type (lat_struct) lat
type (ele_struct) ele
type (coord_struct) orb0

integer nargs

character(40) :: lat_file = 'spin_test.bmad'
logical print_extra, err_flag
                   
!                  

global_com%exit_on_error = .false.

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

call bmad_parser (lat_file, lat)


!

call init_coord (orb0, lat%beam_start, lat%ele(0), downstream_end$)
call transfer_map_calc_with_spin (lat, ele%taylor, ele%spin_taylor, orb0, err_flag, 0, 1)

call type_taylors (ele%taylor)
print *, '--------------------------------'
call type_spin_taylors (ele%spin_taylor)


end program
