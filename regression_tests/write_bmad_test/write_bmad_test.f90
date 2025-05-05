!+
! Program write_bmad_test
!
! This program is part of the Bmad regression testing suite.
!-

program write_bmad_test

use bmad

implicit none

type (lat_struct), target :: lat, lat2
type (ele_struct), pointer :: ele, ele2
type (branch_struct), pointer :: branch, branch2
integer ib, ie

!

open (1, file = 'output.now')

call bmad_parser ('write_bmad_test.bmad', lat, make_mats6 = .false.)
call write_bmad_lattice_file ('z.bmad', lat)
call bmad_parser ('z.bmad', lat2, make_mats6 = .false.)
call bmad_parser ('z.bmad', lat2, make_mats6 = .false.)   ! Make sure to read digested file

do ib = 0, ubound(lat%branch,1)
  branch => lat%branch(ib)
  branch2 => lat2%branch(ib)
  do ie = 1, branch%n_ele_max
    ele  => branch%ele(ie)
    ele2 => branch2%ele(ie)
    call compare_eles(ele, ele2)
  enddo
enddo

close(1)

!--------------------------------------------------------------------------------------------
contains

subroutine compare_eles(ele, ele2)

type (ele_struct) ele, ele2

real(rp) val(num_ele_attrib$), val2(num_ele_attrib$)
integer i
character(40) attrib_name
logical good

!

val = ele%value
val2 = ele2%value
good = .true.

do i = 1, num_ele_attrib$
  attrib_name = attribute_name(ele, i)
  if (attrib_name(1:1) == '!') cycle  ! Note a valid attribute
  if (is_real_matched(val(i), val2(i), ele, attrib_name)) cycle
  write (1, '(3a)') quote(trim(ele%name)), ' STR ', quote(attrib_name)
  good = .false.
enddo

if (good) write (1, '(3a)') quote(trim(ele%name)), ' STR  "GOOD"'

end subroutine compare_eles

!--------------------------------------------------------------------------------------------
! contains

function is_real_matched(val1, val2, ele, name) result (is_matched)

type (ele_struct) ele
real(rp) val1, val2
character(*) name
logical is_matched

!

is_matched = .true.
if (val1 == 0 .and. val2 == 0) return
if (abs(val1 - val2) < 1d-14 * (abs(val1) + abs(val2))) return
is_matched = .false.

end function is_real_matched

end program
