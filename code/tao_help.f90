!+
! Subroutine tao_help (what1, what2)
!
! Online help for TAO commmands. 
! Interfaces with the documentation.
!
! Input:
!   what1   -- Character(*): command to query
!
! Output:
!   none
!
!-

subroutine tao_help (what1, what2)

use tao_struct
use tao_interface
use sim_utils

implicit none

integer i, nl, iu, ios, n, ix, ix2

character(*) :: what1, what2
character(16) :: r_name = "TAO_HELP"
character(40) start_tag, left_over_eliminate, left_over_sub
character(200) line, file_name, full_file_name

logical blank_line_before, in_example

! Help depends upon if we are in single mode or not.
! Determine what file to open and starting tag.

if (tao_com%single_mode) then
  file_name = 'TAO_DIR:doc/single-mode.tex'
else
  file_name = 'TAO_DIR:doc/command-list.tex'
endif

call fullfilename (file_name, full_file_name)

if (what1 == '') then
  start_tag = '%% command_table'
else
  start_tag = '%% ' // what1
endif

! Open the file 

iu = lunget()
open (iu, file = full_file_name, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // file_name)
  return
endif

! Skip all lines before the start tag.

n = len_trim(start_tag)
do 
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) then
    call out_io (s_error$, r_name, &
           'CANNOT FIND ANY INFO FOR: ' // trim(what1) // ' ' // trim(what2), &
           'IN FILE: ' // file_name)
    return
  endif
  ! If a match for what1 then check for a match for what2.
  if (line(1:n) == start_tag(1:n)) then
    if (what2 == '') exit
    if (line(n+1:n+1) == ' ') then
      call string_trim(line(n+1:), line, ix)
    else
      call string_trim(line(n+1:), line, ix)
      call string_trim(line(ix+1:), line, ix)
    endif
    if (index(line, trim(what2)) == 1) exit
  endif
enddo

! Print all lines to the next tag or the end of the file.

if (what1 == '') then
  call out_io (s_blank$, r_name, &
                 "Type 'help <command>' for help on an individual command", &
                 "Available commands:")
endif

blank_line_before = .true.
left_over_eliminate = ''
left_over_sub = ''
in_example = .false.

do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) return
  if (line(1:2) == '%%') return

  if (index(line, '\begin{example}') /= 0) in_example = .true.
  if (index(line, '\end{example}') /= 0) in_example = .true.

  if (line(1:8)  == '\section')   cycle
  if (line(1:6)  == '\label')     cycle
  if (line(1:6)  == '\begin')     cycle
  if (line(1:4)  == '\end')       cycle
  if (line(1:10) == '\centering') cycle
  if (line(1:8)  == '\caption') cycle
  
  if (line(1:6)  == '\vskip') then
    call string_trim (line(7:), line, ix)
    call string_trim (line(ix+1:), line, ix)
  endif

  if (left_over_eliminate /= '') then
    ix = index(line, trim(left_over_eliminate))
    if (ix /= 0) then
      call substitute (left_over_eliminate, left_over_sub)
      left_over_eliminate = ''
      left_over_sub = ''
    endif
  endif

  call substitute ("``", '"')
  call substitute ("''", '"')
  call substitute ("$")
  call substitute ("\protect")
  call substitute ("\_", "_")
  call substitute ("\#", "#")
  call substitute ("\tao", "Tao")
  call eliminate2 ('\item[', ']', '     ', '')
  call eliminate2 ('\vn{', '}', '"', '"')
  call eliminate_inbetween ('& \sref{', '}', .true.)
  call eliminate_inbetween ('\hspace*{', '}', .true.)
  call eliminate_inbetween ('\sref{', '}', .false.)
  call eliminate_inbetween ('{\it ', '}', .false.)
  call eliminate_inbetween ('\parbox{', '}', .false.)
  call substitute ("] \Newline")
  call substitute ("\Newline")
  if (.not. in_example) call substitute (" &")
  call substitute ('\\ \hline')
  call substitute ('\\')
  call substitute ('\W ', '^')
  call substitute ('"\W"', '"^"')

  do
    if (line(1:1) /= '{' .and. line(1:1) /= '}') exit
    line = line(2:)
  enddo
  do i = 2, len(line)
    if (line(i-1:i-1) == '\') cycle !'
    if (line(i:i) == '{' .or. line(i:i) == '}') line = line(:i-1) // line(i+1:)
  enddo

  call substitute ("\{", "{")
  call substitute ("\}", "}")
  
  if (line == ' ') then
    if (blank_line_before) cycle
    blank_line_before = .true.
  else
    blank_line_before = .false.
  endif

  call out_io (s_blank$, r_name, line)

enddo

close (iu)

!-----------------------------------------------------------------------------
contains
!
! substitutes a string and optionally replaces it with another

subroutine substitute (str1, sub)

character(*) str1
character(*), optional :: sub
integer n1

!

n1 = len(str1)

do
  ix = index(line, str1)
  if (ix == 0) exit
  if (present(sub)) then
    line = line(1:ix-1) // trim(sub) // line(ix+n1:)
  else
    line = line(1:ix-1) // line(ix+n1:)
  endif
enddo

end subroutine


!-----------------------------------------------------------------------------
! contains
!
! eliminates two strings, but only if they both exist on the same line

subroutine eliminate2 (str1, str2, sub1, sub2)

character(*) str1, str2
character(*), optional :: sub1, sub2
integer n1, n2, ix1, ix2

n1 = len(str1)
n2 = len(str2)

do
  ix1 = index (line, str1)
  if (ix1 == 0) return
  ix2 = index (line(ix1+1:), str2) + ix1

  ! If ending string is not found then must be on a later line.
  ! If so, mark for future deletion

  if (ix2 == ix1) then
    left_over_eliminate = str2
    if (present(sub2)) left_over_sub = sub2 
    if (present(sub1)) then
      line = line(1:ix1-1) // sub1 // line(ix1+n1:)
    else
      line = line(1:ix1-1) // line(ix1+n1:)
    endif
    return
  endif

  !

  if (present(sub1)) then
    line = line(1:ix1-1) // sub1 // line(ix1+n1:ix2-1) // sub2 // line(ix2+n2:)    
  else
    line = line(1:ix1-1) // line(ix1+n1:ix2-1) // line(ix2+n2:)
  endif

enddo

end subroutine

!-----------------------------------------------------------------------------
! contains
!
! eliminates everything between strings, including the strings

subroutine eliminate_inbetween (str1, str2, pad_with_blanks)

character(*) str1, str2
character(100) :: blank = ''

integer n1, n2, ix1, ix2

logical pad_with_blanks

!

n1 = len(str1)
n2 = len(str2)

do
  ix1 = index (line, str1)
  if (ix1 == 0) return

  ix2 = index (line(ix1+1:), str2) + ix1
  if (ix2 == ix1) return

  if (pad_with_blanks) then
    line = line(1:ix1-1) // blank(:ix2+n2-ix1) // line(ix2+n2:)
  else
    line = line(1:ix1-1) // line(ix2+n2:)
  endif
enddo

end subroutine

end subroutine tao_help
