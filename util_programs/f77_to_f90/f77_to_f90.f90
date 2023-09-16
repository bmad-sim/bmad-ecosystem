!+
! Program F77_TO_F90.F90
!
! Program to convert .F77 files to .F90 format.
! Additionally the program will "clean-up" .F90 files by
! replacing "structure" witH "type" etc.
!-

program f77_to_f90

use sim_utils

implicit none

character(80) file_name

integer i, n

!-------------------------------------------------------------------------
! get the input file and open input and output files


n = command_argument_count()

if (n == 0) then
  print '(a, $)', ' Input f77 (or f90) file name: '
  read (*, '(a)') file_name
  call convert1 (file_name)
  stop
endif

do i = 1, n
  call get_command_argument (1, file_name)
  call convert1 (file_name)
enddo

end program

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine convert1 (file_name)

implicit none

integer ix, ii, ios, ix2

character(*) file_name
character(132) line1, line2, line_f90, line4
character(100) f90_file
character(1) tab

integer :: f_type, f77$ = 1, f90$ = 2
integer istat, nret

parameter (tab = char(9))

logical :: is_there, doit, c2, c3, debug_out = .false.

! if no suffix supplied then look for existing file with appropriate suffix

inquire (file = file_name, exist = is_there)

if (.not. is_there) then
  print *, 'ERROR: CANNOT FIND: ', trim(file_name)
  stop
endif

! find input file type 

ix = len_trim(file_name)
if (file_name(ix-3:ix) == '.f77') then
  f_type = f77$
elseif (file_name(ix-1:ix) == '.f') then
  f_type = f77$
elseif (file_name(ix-3:ix) == '.f90') then
  f_type = f90$
else
  print *, 'I do not know what type of file this is: ' // trim(file_name)
  stop
endif

! if f90 input file then make sure we are doing the right thing

if (f_type == f90$) then
!    doit = .true.
!    call logic_get ('Y', 'N', 'F90 to F90 Conversion?', doit)
!    if (.not. doit) then
    print *, 'f90 to f90 conversion disabled.'
    stop
!    endif
endif

call file_suffixer (file_name, f90_file, '.f90', .true.)

print *, 'Output file: ', trim(f90_file)

open (unit = 1, file = file_name, status = 'old')
open (unit = 2, file = f90_file)

!-------------------------------------------------------------------------
! f90_to_f90:

if (f_type == f90$) then
  do
    read (1, '(a)', iostat = ios) line1  
    if (ios /= 0) stop
    do
      ix = index(line1, tab)
      if (ix == 0) exit
      line1 = line1(:ix-1) // '  ' // line1(ix+1:)
    enddo
    call write_line (line1, f_type, .false.)
  enddo
endif

!-------------------------------------------------------------------------
! f77_to_f90
! read a line from the f77 file and transfer it to the f90 file.
! Continuations make things a little complicated:
! We need to check the next line before we write the present line 
!
! LINE1 = f77 input line
! LINE2 = input line with 6 leading blanks trimmed off
! LINE_F90 = f90 output line

ii = 0
read (1, '(a)') line1    ! initialize line1 = input line

! main loop for f77_to_f90

do

  ii = ii + 1

! convert tab form to non-tab form

  if (line1(1:1) == tab) then
    if (index('0123456789', line1(2:2)) == 0) then   ! not a continuation
      line1 = '      ' // line1(2:)
    else   ! a continuation line
      line1 = '     ' // line1(2:)
    endif
  endif

! convert other tabs

  do
    ix = index(line1, tab)
    if (ix == 0) exit
    if (ix < 7 .and. line1(:ix-1) == ' ') then
      line1 = '      ' // line1(ix+1:)
    else
      line1 = line1(:ix-1) // '  ' // line1(ix+1:)
    endif
  enddo

  c2 = .false.

  line2 = line1

! ' type *' -> ' print *', etc.

  ix =index(line2, ' type *')
  if (ix > 0 .and. ix < 10) line2 = line2(1:ix) // 'print *' // line2(ix+7:)

  ix = index(line2, ' type "')
  if (ix > 0 .and. ix < 10) line2 = line2(1:ix) // 'print "' // line2(ix+7:)

  ix = index(line2, " type '")
  if (ix > 0 .and. ix < 10) line2 = line2(1:ix) // "print '" // line2(ix+7:)

! ' character*nn' -> ' character(nn)'

  ix = index(line2, ' character*(*)')
  if (ix > 0 .and. ix < 10) line2 = line2(1:ix) // 'character(*)' // line2(ix+13:)

  ix = index(line2, ' character*')
  if (ix > 0 .and. ix < 10) then
    ix2 = index(line2(ix+11:), ' ')
    line2 = line2(1:ix) // 'character(' // line2(ix+11:ix+9+ix2) // ')' // line2(ix+10+ix2:)
  endif

! if a normal line then just trim off the leading blanks

  if (line1(1:6) == '      ') then     
    line2 = trim(line2(7:))              ! line2 = trimed input line
    
! if a continuation found then must add "&"
       
  elseif (line1(:5) == '     ' .and. line1(6:6) /= '!') then 
    call add_ampersand (line_f90)
    line2 = '        ' // line2(7:)        ! trim off continuation character
    c2 = .true.

! comment lines with a "C" get changed to "!"

  elseif (line1(1:1) == 'c' .or. line1(1:1) == 'C') then
    line2 = '!' // line2(2:)

  endif

! write to output
! but do not output OPTIONS line

  call string_trim (line_f90, line4, ix)
  call str_upcase (line4, line4)
  if (line4(1:7) /= 'OPTIONS') then  
    if (ii .ge. 2) call write_line (line_f90, f_type, c3)
  endif

  line_f90 = line2
  c3 = c2

  read (1, '(a)', end = 1000) line1  

enddo

! Finish. Write final pending f90 line.

1000 continue

call write_line (line_f90, f_type, c3)

end 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine write_line (line, f_type, continue_line)

implicit none
          
integer i, j, ix, i_quote
integer :: f_type, f77$ = 1, f90$ = 2
integer :: ix1, ix2, ix_end, i_stk = 0, i_stk2 = 0

character*80 stack(100), stack2(100), line2, name

character(*) line
character strng*120, str*8
character*52 :: a_to_z = &
      'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'
character*53 :: a_to_p = &
      ')abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'
character*10 :: num = '0123456789'
character*4 :: old_eq(6) = &
      (/ '.eq.', '.ne.', '.gt.', '.lt.', '.ge.', '.le.' /)
character*4 :: old2_eq(6) = &
      (/ '.EQ.', '.NE.', '.GT.', '.LT.', '.GE.', '.LE.' /)
character*2 :: new_eq(6) = &
      (/ '==',   '/=',   '> ',    '< ',    '>=',   '<=' /)

logical :: is_there, initial_if_here, first_to = .true.
logical continue_line

!

if (line(3:11) == 'options /') return

if (line(1:12) == '  subroutine' .or. line(1:9) == '  program') then
  write (2, '(a)') trim(line(3:))
  return
endif

! change ".eq.", ".ne.", ".lt.", ".gt.", ".le.", ".ge."

do i = 1, 6
  do
    ix = max(index(line, old_eq(i)), index(line, old2_eq(i)))
    if (ix == 0) exit
    line = line(:ix-1) // ' ' // trim(new_eq(i)) // ' ' // line(ix+4:)
  enddo
enddo

! change "." to "%" for structures

i_quote = 0
i_loop: do i = 1, len_trim(line)
  if (line(i:i) == "'") then
    i_quote = i_quote + 1
  elseif (line(i:i) == "!") then
    if (mod(i_quote, 2) == 0) exit
  elseif (line(i:i) == ".") then
    if (mod(i_quote, 2) /= 0) cycle
    call str_upcase (str, line(i:i+4))
    if (str(1:5) == '.AND.') cycle 
    if (str(1:4) == '.OR.') cycle 
    if (str(1:5) == '.XOR.') cycle 
    if (str(1:5) == '.NOT.') cycle 
    call str_upcase (str, line(i-4:i))
    if (str(1:5) == '.AND.') cycle 
    if (str(2:5) == '.OR.') cycle 
    if (str(1:5) == '.XOR.') cycle 
    if (str(1:5) == '.NOT.') cycle 
    if (index(a_to_z, line(i+1:i+1)) == 0) cycle
    if (index(a_to_p, line(i-1:i-1)) == 0) then
      do j = i-1, 2, -1
        if (index(num, line(j:j)) == 0) cycle i_loop
        if (index(a_to_z, line(j-1:j-1)) /= 0) exit
      enddo
      if (j == 1) cycle 
    endif
    line(i:i) = "%"
  endif
enddo i_loop

! change "structure" to "type"

ix = max(index(line, 'end structure'), index(line, 'endstructure'))
if (ix /= 0) then
  line = '  end type'
else
  ix = index (line, 'structure')
  if (ix /= 0 .and. line(1:ix-1) == ' ') then
    ix1 = index(line, '/') + 1
    strng = line(ix1:)
    ix1 = index(strng, '/') - 1
    strng = strng(:ix1)   
    call string_trim (strng, strng, ix)
    line = '  type ' // trim(strng)
  endif
endif

! change "record" to "type"

ix = index(line, 'record')
if (ix /= 0 .and. line(1:ix-1) == ' ') then
  ix1 = index(line, '/') + 1
  strng = line(ix1:)
  call string_trim (strng, strng, ix)
  ix1 = index(strng, '/') 
  line = '  type (' // trim(strng(:ix1-1)) // ') ' // trim(strng(ix1+1:))
endif

! for non-flx input we have finished

write (2, '(a)') trim(line)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine push_stack (stack, i_stk, line)

implicit none

character(*) stack(*), line
integer i_stk

i_stk = i_stk + 1
stack(i_stk) = line

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine add_ampersand (line)

character(*) line
integer i
logical in_char

in_char = .false.
i = len_trim(line)

do
  if (line(i:i) == "'") then
    if (line(i-1:i-1) == "'") then
      i = i - 2
      cycle
    else
      in_char = .not. in_char
    endif
  endif
  if (line(i:i) == '!' .and. .not. in_char) then
    line = line(1:i-1) // '& ' // line(i:)
    return
  endif
  i = i - 1
  if (i == 1) then
    line = trim(line) // ' &'
    return
  endif
enddo

end subroutine
