!+
! Subroutine tao_help (cmd_line)
!
! Online help for TAO commmands. 
! Interfaces with the documentation.
!
! Input:
!   help_what   -- Character(*): command to query
!
! Output:
!   none
!
!-

subroutine tao_help (help_what)

use tao_struct
use tao_interface
use cesr_utils

implicit none

integer nl, iu, ios, n, ix, ix2

character(*) :: help_what
character(16) :: r_name = "TAO_HELP"
character(40) start_tag
character(200) line, file_name

logical blank_line_before

! Help depends upon if we are in single mode or not.
! Determine what file to open and starting tag.

if (s%global%single_mode) then
  call fullfilename ('TAO_DIR:doc/single_mode.tex', file_name)
  start_tag = '%% keys'
else
  call fullfilename ('TAO_DIR:doc/command_list.tex', file_name)
  start_tag = '%% ' // help_what
endif

! Open the file 

iu = lunget()
open (iu, file = file_name, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // file_name)
  return
endif

! Skip all lines before the start tag.

n = len_trim(start_tag)
do 
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'CANNOT FIND TAG: ' // start_tag, &
                                                     'IN FILE: ' // file_name)
    return
  endif
  if (line(1:n) == start_tag(1:n)) exit
enddo

! Print all lines to the next tag or the end of the file.

blank_line_before = .true.
do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) return
  if (line(1:2) == '%%') return

  if (line(1:8) == '\section') cycle
  if (line(1:6) == '\label') cycle
  if (line(1:6) == '\begin') cycle
  if (line(1:4) == '\end') cycle

  call eliminate ("``", '"')
  call eliminate ("''", '"')
  call eliminate ("$")
  call eliminate ("\{", "{")
  call eliminate ("\}", "}")
  call eliminate ("\_", "_")
  call eliminate ("\tao", "Tao")
  call eliminate2 ('\item[', ']')
  call eliminate2 ('\vn{', '}')

  if (line == ' ') then
    if (blank_line_before) cycle
    blank_line_before = .true.
  else
    blank_line_before = .false.
  endif

  call out_io (s_blank$, r_name, line)

enddo

!-------------------------------------------------------------------------------------
contains

subroutine eliminate (str1, sub)

character(*) str1
character(*), optional :: sub
integer n1

!

n1 = len(str1)

do
  ix = index(line, str1)
  if (ix == 0) exit
  if (present(sub)) then
    line = line(1:ix-1) // sub // line(ix+n1:)
  else
    line = line(1:ix-1) // line(ix+n1:)
  endif
enddo

end subroutine


!-------------------------------------------------------------------------------------
! contains

subroutine eliminate2 (str1, str2)

character(*) str1, str2
integer n1, n2

n1 = len(str1)
n2 = len(str2)

do
  ix = index (line, str1)
  if (ix == 0) return
  ix2 = index (line(ix+1:), str2) + ix
  if (ix2 == 0) return
  line = line(1:ix-1) // line(ix+n1:ix2-1) // line(ix2+n2:)
enddo

end subroutine

end subroutine tao_help
