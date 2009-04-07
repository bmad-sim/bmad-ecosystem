!+
! Program F77_TO_F90.F90
!
! Program to convert .F77 .FLX files to .F90 format.
! Additionally the program will "clean-up" .F90 files by
! replacing "structure" witH "type" etc.
!
! Known FLEX problems:
!
! 1) If a line has blanks in cols 1-8 and a number in column 9 then
!    FLEX will take this as a continuation line. Some FLEX code use this
!    feature! 
!
! 2) "D ..." debug lines are translated to "if (debug) ..."
!-

program f77_to_f90

use cesr_utils

implicit none

character(80) file_name

integer i, n

!-------------------------------------------------------------------------
! get the input file and open input and output files


n = cesr_iargc()

if (n == 0) then
  print '(a, $)', ' Input f77 (or flx or f90) file name: '
  read (*, '(a)') file_name
  call convert1 (file_name)
  stop
endif

do i = 1, n
  call cesr_getarg (1, file_name)
  call convert1 (file_name)
enddo

end program

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine convert1 (file_name)

implicit none

integer ix, ii, ios

character(*) file_name
character*132 line1, line2, line_f90, line4
character*100 f90_file
character*1 tab

integer :: f_type, f77$ = 1, f90$ = 2, flx$ = 3
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
if (file_name(ix-3:ix) == '.flx') then
  f_type = flx$
elseif (file_name(ix-3:ix) == '.f77') then
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
open (unit = 2, file = f90_file, status = 'new')

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
    call write_line (line1, f_type)
  enddo
endif

!-------------------------------------------------------------------------
! f77_to_f90, flx_to_f90:
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

! flex #9 continuation?

  if (f_type == flx$ .and. line1(1:8) == ' ' .and. &
                                index('123456789', line1(9:9)) /= 0) then
    call add_ampersand (line_f90)
    line2 = line1(10:)        ! trim off continuation character
    c2 = .true.

! if a normal line then just trim off the leading blanks

  elseif (line1(1:6) == '      ') then     
    line2 = trim(line1(7:))              ! line2 = trimed input line
    
! if a continuation found then must add "&"
       
  elseif (line1(:5) == '     ' .and. line1(6:6) /= '!') then 
    call add_ampersand (line_f90)
    line2 = '        ' // line1(7:)        ! trim off continuation character
    c2 = .true.

! comment lines with a "C" get changed to "!"

  elseif (line1(1:1) == 'c' .or. line1(1:1) == 'C') then
    line2 = '!' // line1(2:)

! flex "D" lines converted to "if (debug)"

  elseif (line1(1:1) == 'D' .and. f_type == flx$) then
    call string_trim (line1(2:), line1, ix)
    line2 = '  if (debug) ' // line1
    if (.not. debug_out) then
      print *, 'WARNING! THIS FILE USES THE FLEX "D ..." DEBUG FEATURE.'
      print *, '         THIS HAS BEEN CONVERTED TO "IF (DEBUG) ...".'
      print *, '         YOU NEED TO ADD A "LOGICAL DEBUG" TO THE CODE.'
      debug_out = .true.
    endif

! only possibility left is a comment or blank line

  else
    line2 = line1
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
integer :: f_type, f77$ = 1, f90$ = 2, flx$ = 3
integer :: ix1, ix2, ix_end, i_stk = 0, i_stk2 = 0

character*80 stack(100), stack2(100), line2, name

character*(*) line
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

if (f_type /= flx$) then
  write (2, '(a)') trim(line)
  return
endif

!--------------------------------------------------------------
! here when we have a flx input file

call string_trim (line, line, ix)

! continue_line

if (continue_line) then
  call write_flx (line, i_stk-i_stk2)
  return
endif

! comment line

if (line(1:1) == '!') then
  write (2, '(a)') trim(line)
  return
endif

! if there is a line number make this a continue statement
! except if this is a format statement

if (index('0123456789', line(1:1)) /= 0) then
  ix = index(line, ' ')
  call string_trim (line(ix:), strng, ix)
  call str_upcase (strng, strng)
  if (strng(1:6) == 'FORMAT') then
    write (2, '(a)') trim(line)
    return
  else
    write (2, '(a)') line(1:5) // ' continue'
    call string_trim (line(6:), line, ix)
  endif
endif

! IF construct

call flx_word_test (line, 'IF', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0 .and. ix_end == 0) then
  strng = trim(line(:ix2)) // ' then ' // line(ix2+1:)
  call write_flx (strng, i_stk-i_stk2)
  call push_stack (stack, i_stk, 'end if')
  return
endif

! UNLESS

call flx_word_test (line, 'UNLESS', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  if (ix_end == 0) then
    line = 'if (.not. ' // line(ix1+1:ix2) // ' then ' // line(ix2+1:)
    call write_flx (line, i_stk-i_stk2)
    call push_stack (stack, i_stk, 'end if')
  else
    line = 'if (.not. ' // line(ix1+1:)
    call write_flx (line, i_stk-i_stk2)
  endif
  return
endif

! WHEN ... 

call flx_word_test (line, 'WHEN', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  if (ix_end == 0) then
    strng = 'if ' // line(ix1:ix2) // ' then ' // line(ix2+1:)
    call write_flx (strng, i_stk-i_stk2)
    call push_stack (stack, i_stk, 'NULL')
  else
    strng = 'if ' // line(ix1:ix2) // ' then '
    call write_flx (strng, i_stk-i_stk2)
    call write_flx (line(ix_end:), i_stk+1-i_stk2)
  endif
  return
endif

! ... ELSE

call flx_word_test (line, 'ELSE', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 == 0) then
  if (ix_end == 0) then
    call write_flx (line, i_stk-i_stk2)
    call push_stack (stack, i_stk, 'end if')
  else
    call write_flx ('else', i_stk-i_stk2)
    call write_flx (line(ix_end:), i_stk+1-i_stk2)
    call write_flx ('end if', i_stk-i_stk2)
  endif
  return
endif

! CONDITIONAL ...

call flx_word_test (line, 'CONDITIONAL', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 == 0) then
  call push_stack (stack, i_stk, 'CONDITIONAL')
  call push_stack (stack2, i_stk2, 'CONDITIONAL')
  initial_if_here = .true.
  return
endif

! SELECT ... 

call flx_word_test (line, 'SELECT', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  call push_stack (stack, i_stk, 'SELECT')
  call push_stack (stack2, i_stk2, 'SELECT')
  call write_flx ('select case ' // line(ix1:), i_stk-i_stk2)
  return
endif

! ... ()

call flx_word_test (line, '()', is_there, ix1, ix2, ix_end)

if (is_there) then

  if (stack2(i_stk2) == 'CONDITIONAL') then
    if (line(ix1:ix2) == '(OTHERWISE)') then
      line2 = 'else'
    elseif (initial_if_here) then
      line2 = 'if ' // line(ix1:ix2) // ' then'
      initial_if_here = .false.
    else
      line2 = 'elseif ' // line(ix1:ix2) // ' then'
    endif
  else
    if (line(ix1:ix2) == '(OTHERWISE)') then
      line2 = 'case default'
    else
      line2 = 'case ' // line(ix1:ix2)
    endif
  endif

  call write_flx (line2, i_stk-i_stk2)

  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'NULL')
  else
    call write_flx (line(ix_end:), i_stk+1-i_stk2)
  endif

  return

endif

! DO

call flx_word_test (line, 'DO', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  call write_flx ('do ' // line(ix1+1:ix2-1), i_stk-i_stk2)
  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'end do')
  else
    call write_flx (line2(ix_end:), i_stk+1-i_stk2)
    call write_flx ('end do', i_stk-i_stk2)
  endif
  return
endif

! WHILE

call flx_word_test (line, 'WHILE', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  call write_flx ('do while' // line(ix1:ix2), i_stk-i_stk2)
  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'end do')
  else
    call write_flx (line2(ix_end:), i_stk+1-i_stk2)
    call write_flx ('end do', i_stk-i_stk2)
  endif
  return
endif

! REPEAT WHILE

call flx_word_test (line, 'REPEAT WHILE', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  call write_flx ('do', i_stk-i_stk2)
  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'REPEAT ' // line(ix1:ix2))
  else
    call write_flx (line(ix_end:), i_stk+1-i_stk2)
    call write_flx ('if ' // line(ix1:ix2) // ' exit', i_stk+1-i_stk2)
    call write_flx ('end do', i_stk-i_stk2)
  endif
  return
endif

! UNTIL

call flx_word_test (line, 'UNTIL', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  call write_flx ('do while (.not. ' // line(ix1+1:ix2), i_stk-i_stk2)
  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'end do')
  else
    call write_flx (line(ix_end:), i_stk+1-i_stk2)
    call write_flx ('end do', i_stk-i_stk2)
  endif
  return
endif

! REPEAT UNTIL

call flx_word_test (line, 'REPEAT UNTIL', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 /= 0) then
  call write_flx ('do', i_stk-i_stk2)
  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'REPEAT (.not. ' // line(ix1+1:ix2))
  else
    call write_flx (line(ix_end:), i_stk+1-i_stk2)
    call write_flx ('if (.not. ' // line(ix1+1:ix2) // ' exit', &
                                                          i_stk+1-i_stk2)
    call write_flx ('end do', i_stk-i_stk2)
  endif
  return
endif
 
! Call Internal Procedure

call flx_procedure_test (line, name, is_there, ix_end)

if (is_there .and. ix_end == 0) then
  call write_flx ('call ' // name, i_stk-i_stk2)
  return
endif

! TO Internal Procedure

call flx_word_test (line, 'TO', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 == 0 .and. line(ix_end:ix_end) /= '=') then

  if (first_to) then
    write (2, *)
    write (2, '(a)') '  contains'
    first_to = .false.
  endif

  write (2, *)
  write (2, '(a)') ' !-------------------------------------------------'
  write (2, *)

  call flx_procedure_test (line(ix_end:), name, is_there, ix_end)
  if (.not. is_there) then
    print *, 'ERROR: NO PROCEDURE NAME AFTER "TO" IN:'
    print *, '    ', line
    call err_exit
  endif

  call write_flx ('subroutine ' // name, i_stk-i_stk2)

  if (ix_end == 0) then
    call push_stack (stack, i_stk, 'end subroutine')
  else
    call write_flx (line(ix_end:), i_stk-i_stk2)
    write (2, '(a)') ' end subroutine'
  endif

  return

endif

! FIN

call flx_word_test (line, 'FIN', is_there, ix1, ix2, ix_end)

if (is_there .and. ix2 == 0 .and. ix_end == 0) then
  strng = stack(i_stk)
  i_stk = i_stk - 1
  if (strng == 'CONDITIONAL') then
    i_stk2 = i_stk2 - 1
    call write_flx ('end if', i_stk-i_stk2)
  elseif (strng == 'SELECT') then
    i_stk2 = i_stk2 - 1
    call write_flx ('end select', i_stk-i_stk2)
  elseif (strng(1:6) == 'REPEAT') then
    call write_flx ('if ' // trim(strng(8:)) // ' exit', i_stk+1-i_stk2)
    call write_flx ('end do', i_stk-i_stk2)
  elseif (strng /= 'NULL') then
    call write_flx (strng, i_stk-i_stk2)
  endif
  return
endif

! write the line

call write_flx (line, i_stk-i_stk2)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine flx_word_test (line, name, is_there, ix1, ix2, ix_end)

implicit none

character*(*) line, name
character*80 nam

integer ix1, ix2, ix_end                                
integer i, j, ix_p, ix0, ix

logical is_there, there

! init default

is_there = .false.
ix1 = 0; ix2 = 0; ix_end = 0

! look for name
                                    
if (name == '()') then
  ix0 = 1
else
  ix = index (line, name)
  if (ix /= 1) return
  ix0 = len(name) + 1
  if (line(ix0:ix0) /= ' ' .and. line(ix0:ix0) /= '(') return
  is_there = .true.
endif

! look for opening "("

if (ix0 > len_trim(line)) return

do i = ix0, len_trim(line)
  if (line(i:i) == '(') then
    if (name == '()') is_there = .true.
    ix1 = i
    exit
  elseif (line(i:i) /= ' ') then
    if (line(i:i) /= '!') then
      ix_end = i
      if (name /= 'TO') then
        call flx_procedure_test (line(ix_end:), nam, there, ix)
        if (there) line(ix_end:) = 'call ' // line(ix_end:)
      endif
    endif
    return
  endif
enddo

! look for ending ")"

ix_p = 1
do j = i+1, len_trim(line)
  if (line(j:j) == '(') ix_p = ix_p + 1
  if (line(j:j) == ')') ix_p = ix_p - 1
  if (ix_p == 0) then
    ix2 = j
    exit
  endif
enddo
 
if (ix2 == 0) then
  print *, 'ERROR: CANNOT FIND ENDING ")" IN:'
  print *, trim(line)
  call err_exit
endif

do j = ix2+1, len_trim(line)
  if (line(j:j) /= ' ') then
    if (line(j:j) /= '!') then
      ix_end = j
      call flx_procedure_test (line(ix_end:), nam, there, ix)
      if (there) line(ix_end:) = 'call ' // line(ix_end:)
    endif
    return
  endif
enddo

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine flx_procedure_test (line, name, is_there, ix_end)

implicit none

character*(*) line, name

integer ix_end, i, j, ix

logical is_there, dash_found

!

is_there = .false.
ix_end = 0

dash_found = .false.

name = line

do j = 1, len_trim(line)
  if (line(j:j) == ' ') then
    if (dash_found) then
      name(j:j) = ' '
      exit
    else
      return
    endif
  elseif (line(j:j) == '-') then
    dash_found = .true.
    name(j:j) = '_'
  elseif (index('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', line(j:j)) == 0) then
    return
  endif
enddo

if (dash_found) then
  is_there = .true.
  ix = len_trim(name)
  line(1:ix) = name(1:ix)
endif

do i = j, len_trim(line)
  if (line(i:i) /= ' ') then
    if (line(i:i) /= '!') ix_end = i
    return
  endif
enddo

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine write_flx (line, i_stk)

implicit none

character*(*) line
integer i_stk, n

n = 2*i_stk + 2
write (2, '(<n>x, a)') trim(line)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine push_stack (stack, i_stk, line)

implicit none

character*(*) stack(*), line
integer i_stk

i_stk = i_stk + 1
stack(i_stk) = line

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine add_ampersand (line)

character*(*) line
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
