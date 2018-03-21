program field_test

use bmad

implicit none

type (lat_struct) lat


!

call bmad_parser ('field_test.bmad', lat)
print *, mat_symp_error(lat%ele(1)%mat6)

end program
