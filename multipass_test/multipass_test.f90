program test

use bmad

implicit none

type (lat_struct) lat

!

call bmad_parser ('lat.bmad', lat, use_line = '')

open (1, file = 'output.now')

write (1, '(3a)') '"01"  STR  "', trim(lat%ele(01)%name), '"'
write (1, '(3a)') '"02"  STR  "', trim(lat%ele(02)%name), '"'
write (1, '(3a)') '"03"  STR  "', trim(lat%ele(03)%name), '"'
write (1, '(3a)') '"04"  STR  "', trim(lat%ele(04)%name), '"'
write (1, '(3a)') '"05"  STR  "', trim(lat%ele(05)%name), '"'
write (1, '(3a)') '"06"  STR  "', trim(lat%ele(06)%name), '"'
write (1, '(3a)') '"07"  STR  "', trim(lat%ele(07)%name), '"'
write (1, '(3a)') '"08"  STR  "', trim(lat%ele(08)%name), '"'
write (1, '(3a)') '"09"  STR  "', trim(lat%ele(09)%name), '"'
write (1, '(3a)') '"10"  STR  "', trim(lat%ele(10)%name), '"'
write (1, '(3a)') '"11"  STR  "', trim(lat%ele(11)%name), '"'
write (1, '(3a)') '"12"  STR  "', trim(lat%ele(12)%name), '"'
write (1, '(3a)') '"13"  STR  "', trim(lat%ele(13)%name), '"'

close (1)

end program
