program ptc_to_bmad

use bmad

implicit none

type (lat_struct) lat
integer i, n
logical err_flag
character(200) flat_file(2)

!

n = command_argument_count()

if (n < 1 .or. n > 2) then
  print *, 'Usage:'
  print *, '  ptc_to_bmad ptc_flat_file1 ptc_flat_file2'
  print *, 'Note: flat_file2 is optional.'
  stop
endif

!

do i = 1, n
  call get_command_argument(i, flat_file(i))
enddo

call ptc_read_flat_file(flat_file(1:n), err_flag, lat)

end program
