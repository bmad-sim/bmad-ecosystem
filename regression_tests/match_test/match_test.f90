program match_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer ib, ie
character(100) lat_file

!

open (1, file = 'output.now')

lat_file = 'match_test.bmad'
call bmad_parser(lat_file, lat)

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele => branch%ele(ie)

enddo


end program
