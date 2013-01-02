program test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

! Init

open (1, file = 'output.now')

! Forking with a branch element

call bmad_parser ('branch_fork.bmad', lat)
call write_bmad_lattice_file ('lat.bmad', lat)
call bmad_parser ('lat.bmad', lat)

ele => lat%branch(1)%ele(2)
write (1, '(3a)')       '"BF-01"  STR  "', trim(ele%name), '"'
write (1, '(a, 2i3)')   '"BF-02"  ABS   0 ', nint(ele%value(ix_to_branch$)), nint(ele%value(ix_to_element$))
write (1, '(a, f10.4)') '"BF-03"  REL  1e-12 ', ele%value(p0c$)

! Multipass and superimpose

call bmad_parser ('multipass_and_superimpose.bmad', lat)

write (1, '(3a)') '"MS-01"  STR  "', trim(lat%ele(01)%name), '"'
write (1, '(3a)') '"MS-02"  STR  "', trim(lat%ele(02)%name), '"'
write (1, '(3a)') '"MS-03"  STR  "', trim(lat%ele(03)%name), '"'
write (1, '(3a)') '"MS-04"  STR  "', trim(lat%ele(04)%name), '"'
write (1, '(3a)') '"MS-05"  STR  "', trim(lat%ele(05)%name), '"'
write (1, '(3a)') '"MS-06"  STR  "', trim(lat%ele(06)%name), '"'
write (1, '(3a)') '"MS-07"  STR  "', trim(lat%ele(07)%name), '"'
write (1, '(3a)') '"MS-08"  STR  "', trim(lat%ele(08)%name), '"'
write (1, '(3a)') '"MS-09"  STR  "', trim(lat%ele(09)%name), '"'
write (1, '(3a)') '"MS-10"  STR  "', trim(lat%ele(10)%name), '"'
write (1, '(3a)') '"MS-11"  STR  "', trim(lat%ele(11)%name), '"'
write (1, '(3a)') '"MS-12"  STR  "', trim(lat%ele(12)%name), '"'
write (1, '(3a)') '"MS-13"  STR  "', trim(lat%ele(13)%name), '"'

!


! And close

close (1)

end program
