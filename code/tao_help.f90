!+
! Subroutine tao_help (what1, what2, lines, n_lines)
!
! Online help for TAO commmands. 
!
! Input:
!   what1   -- Character(*): command to query. EG: "show".
!   what2   -- Character(*): subcommand to query. EG: "element".
!
! Output:
!   lines(:) -- character(200), optional, allocatable: If present then the output will 
!                 be put in this string array instead of printing to the terminal.
!   n_lines  -- integer, optional: Must be present if lines is present.
!                 Number of lines used in the lines(:) array.
!-

subroutine tao_help (what1, what2, lines, n_lines)

use tao_struct
use tao_interface, dummy => tao_help

implicit none

integer, optional :: n_lines
integer i, iu, ios, n, ix, ix2, nl

character(*) :: what1, what2
character(*), parameter :: r_name = "tao_help"
character(40) start_tag, left_over_eliminate, left_over_sub
character(200) line, file_name, full_file_name
character(*), optional, allocatable :: lines(:)

logical blank_line_before, in_example, has_subbed, python_search

! This help system depends upon parsing one of three files:
!       tao/doc/single-mode.tex
!       tao/doc/command-list.tex
!       tao/code/tao_python_cmd.f90
! The code here will look for the appropriate string (start_tag) that signals that
! the wanted documentation has been found.

! Help depends upon if we are in single mode or not.
! Determine what file to open and starting tag.

python_search = .false.

if (index('python', trim(what1)) == 1 .and. what2 /= '') then
  file_name = '$TAO_DIR/code/tao_python_cmd.f90'
  python_search = .true.
elseif (s%com%single_mode) then
  file_name = '$TAO_DIR/doc/single-mode.tex'
else
  file_name = '$TAO_DIR/doc/command-list.tex'
endif

call fullfilename (file_name, full_file_name)

if (python_search) then
  start_tag = '!%% ' // what2
elseif (what1 == '' .or. what1 == 'help-list') then
  start_tag = '%% command_table'
else
  start_tag = '%% ' // what1
endif

! Open the file 

iu = lunget()
open (iu, file = full_file_name, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // file_name, &
                                 'TYPICALLY THIS IS DUE TO THE ENVIRONMENT VARIABLE TAO_DIR NOT BEING SET.', &
                                 'TYPICALLY TAO_DIR LOOKS LIKE: "<path-to-dist>/bmad_dist_2023_0509-0/tao"', &
                                 'WHERE "<path-to-dist>" IS THE PATH TO THE BMAD DISTRIBUTION.')
  return
endif

! Python search

if (python_search) then
  n = len_trim(start_tag)
  ! Find start of desired comment block
  do
    if (.not. read_this_line(line)) return
    if (line(1:n) == start_tag(1:n)) exit
  enddo

  ! Print comment block
  do
    read (iu, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line(1:1) /= '!') exit
    call this_line_out (line)
  enddo

  close (iu)
  if (present(n_lines)) n_lines = nl

  return
endif

! Skip all lines before the start tag.

n = len_trim(start_tag)
do
  if (.not. read_this_line(line)) return

  ! If a match for what1 then check for a match for what2.
  if (line(1:n) == start_tag(1:n)) then
    if (what2 == '') exit
    if (line(n+1:n+1) == ' ') then  ! If what1 exact match.
      call string_trim(line(n+1:), line, ix)
    else  ! Else not an exact match then remove rest of what1 word from line.
      call string_trim(line(n+1:), line, ix)
      call string_trim(line(ix+1:), line, ix)
    endif
    if (index(line, trim(what2)) == 1) exit
    if (line(1:1) == '*') exit    ! "*" means match to anything.
  endif
enddo

! Print all lines to the next tag or the end of the file.

if (what1 == '') then
  call out_io (s_blank$, r_name, &
                 "Type 'help <command>' for help on an individual command", &
                 "Available commands:")
endif

blank_line_before = .true.
left_over_eliminate = ''  ! To handle multiline constructs
left_over_sub = ''        ! To handle multiline constructs
in_example = .false.
nl = 0

do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  if (line(1:2) == '%%') exit

  if (index(line, '\begin{example}') /= 0) in_example = .true.
  if (index(line, '\end{example}') /= 0) in_example = .false.

  if (line(1:8)  == '\section')    cycle
  if (line(1:11) == '\subsection') cycle
  if (line(1:6)  == '\label')      cycle
  if (line(1:6)  == '\begin')      cycle
  if (line(1:4)  == '\end')        cycle
  if (line(1:10) == '\centering')  cycle
  if (line(1:8)  == '\caption')    cycle
  if (line(1:6)  == '\vskip')      cycle
  if (index(line, '\index{') /= 0) cycle

  if (left_over_eliminate /= '') then
    ix = index(line, trim(left_over_eliminate))
    if (ix /= 0) then
      call substitute (line, left_over_eliminate, left_over_sub)
      left_over_eliminate = ''
      left_over_sub = ''
    endif
  endif

  has_subbed = .false.
  call substitute (line, "{}", "")
  call eliminate2 (line, '\vn{', '}', '"', '"')
  call eliminate2 (line, '\item[\vn{\{', '\}}]', '   Argument: ', '')
  call eliminate2 (line, '\item[', ']', '     ', '')
  call substitute (line, '\item', '*')
  call substitute (line, "``\vn", '"')
  call substitute (line, "``", '"')
  call substitute (line, "\`", '`')
  call substitute (line, "}''", '"')
  call substitute (line, "''", '"')
  call substitute (line, "\bf", '')
  call substitute (line, "\arrowbf", '')
  call substitute (line, "\$", '$', has_subbed)
  if (.not. has_subbed .and. .not. in_example) call substitute (line, "$")  ! Do not remove "$" if "\$" -> "$" has been done
  call substitute (line, "\protect")
  call substitute (line, "\_", "_")
  call substitute (line, "\#", "#")
  call substitute (line, "\%", "%")
  call substitute (line, "\tao", "Tao")
  call substitute (line, "\bmad", "Bmad")
  call eliminate_inbetween (line, '& \sref{', '}', .true.)
  call eliminate_inbetween (line, '\hspace*{', '}', .true.)
  call eliminate_inbetween (line, '(\sref{', '})', .false.)
  call eliminate_inbetween (line, ' \sref{', '}', .false.)
  call eliminate_inbetween (line, '\sref{', '}', .false.)
  call eliminate_inbetween (line, '{\it ', '}', .false.)
  call eliminate_inbetween (line, '\parbox{', '}', .false.)
  call substitute (line, "] \Newline")
  call substitute (line, "\Newline")
  if (.not. in_example) call substitute (line, " &")
  call substitute (line, '\vfill')
  call substitute (line, '\vfil')
  call substitute (line, '\hfill')
  call substitute (line, '\hfil')
  call substitute (line, '\break')
  call substitute (line, '\midrule')
  call substitute (line, '\toprule')
  call substitute (line, '\bottomrule')
  call substitute (line, '\\ \hline')
  call substitute (line, '\\')
  call substitute (line, '\W ', '^')
  call substitute (line, '"\W"', '"^"')
  call substitute (line, "\Bf ", "")
  call substitute (line, "\B", "\")       ! "

  if (line(1:2) == '% ') line = line(3:)
  if (line(1:1) == '%')  line = line(2:)

  do
    if (line(1:1) /= '{' .and. line(1:1) /= '}') exit
    line = line(2:)
  enddo

  i = 2
  do while (i /= len(line))
    if ((line(i:i) == '{' .or. line(i:i) == '}') .and. line(i-1:i-1) /= '\') then   ! '
      line = line(:i-1) // line(i+1:)
    else
      i = i + 1
    endif
  enddo

  call substitute (line, "\{", "{")
  call substitute (line, "\}", "}")
  call substitute (line, "\(", "")
  call substitute (line, "\)", "")

  n = max(1, len_trim(line))
  if (line(n:n) == '!') line(n:n) =  ' '

  if (line == ' ') then
    if (blank_line_before) cycle
    blank_line_before = .true.
  else
    blank_line_before = .false.
  endif

  call this_line_out(line)
enddo

close (iu)
if (present(n_lines)) n_lines = nl

!-----------------------------------------------------------------------------
contains

subroutine this_line_out(line)

character(*) line

!

if (present(lines)) then
  if (.not. allocated(lines)) allocate(lines(100))
  if (nl >= size(lines)) call re_allocate (lines, nl+100)
  nl = nl+1; lines(nl) = line
else
  call out_io (s_blank$, r_name, line)
endif

end subroutine this_line_out

!-----------------------------------------------------------------------------
! contains

function read_this_line(line) result (ok)

character(*) line
logical ok

!

read (iu, '(a)', iostat = ios) line
ok = .true.

if (ios /= 0) then
  call out_io (s_error$, r_name, &
         'CANNOT FIND ANY INFO FOR: ' // trim(what1) // ' ' // trim(what2), &
         'IN FILE: ' // file_name)
  close(iu)
  ok = .false.
endif

end function read_this_line

!-----------------------------------------------------------------------------
! contains
!
! Removes a string and optionally replaces it with another.

subroutine substitute (line, str1, sub, has_subbed)

character(*) line, str1
character(*), optional :: sub
integer n1
logical, optional :: has_subbed

!

n1 = len(str1)
if (present(has_subbed)) has_subbed = .false.

do
  ix = index(line, str1)
  if (ix == 0) exit
  if (present(has_subbed))has_subbed = .true.
  if (present(sub)) then
    line = line(1:ix-1) // trim(sub) // line(ix+n1:)
  else
    line = line(1:ix-1) // line(ix+n1:)
  endif
enddo

end subroutine substitute

!-----------------------------------------------------------------------------
! contains
!
! eliminates two strings, but only if they both exist on the same line

subroutine eliminate2 (line, str1, str2, sub1, sub2)

character(*) line, str1, str2
character(*), optional :: sub1, sub2
integer n1, n2, ix1, ix2

n1 = len(str1)
n2 = len(str2)
ix1 = 0

main: do

  ! Find str1

  ix1 = ix1 + 1
  if (ix1+n1-1 > len(line)) return
  if (line(ix1:ix1+n1-1) /= str1) cycle
  if (ix1 > 1) then
    if (line(ix1-1:ix1-1) == '\') cycle   ! '
  endif

  ! Find str2

  ix2 = ix1 + n1 - 1
  do
    ix2 = ix2 + 1

    ! If ending string is not found then must be on a later line.
    ! If so, mark for future deletion

    if (ix2+n2-1 > len(line)) then
      left_over_eliminate = str2
      if (present(sub2)) left_over_sub = sub2 
      if (present(sub1)) then
        line = line(1:ix1-1) // sub1 // line(ix1+n1:)
      else
        line = line(1:ix1-1) // line(ix1+n1:)
      endif
      return
    endif

    if (line(ix2:ix2+n2-1) /= str2) cycle
    if (line(ix2-1:ix2-1) == '\') cycle   ! '
    exit
  enddo

  ! substitute

  if (present(sub1)) then
    line = line(1:ix1-1) // sub1 // line(ix1+n1:ix2-1) // sub2 // line(ix2+n2:)    
    ix1 = ix1 + len(sub1) - 1
  else
    line = line(1:ix1-1) // line(ix1+n1:ix2-1) // line(ix2+n2:)
    ix1 = ix1 - 1
  endif

enddo main

end subroutine eliminate2

!-----------------------------------------------------------------------------
! contains
!
! eliminates everything between strings, including the strings

subroutine eliminate_inbetween (line, str1, str2, pad_with_blanks)

character(*) line, str1, str2
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

end subroutine eliminate_inbetween

end subroutine tao_help
