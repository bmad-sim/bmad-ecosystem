!+
! Program make_f90_interface.f90
!
! Program to make an f90 interface block for a subroutine or a 
! function.
!
! Known "restrictions"
!
! 1) If there is a "contains" then contained subroutine and functions must
! end with "end subroutine" or "end function". Not just with "end".
!
! 2) All lines containing arguments must not be continuation lines.
!-

program make_f90_interface

use dcslib_interface
use directory_mod

implicit none                    

character(132) f90_file, rline(10), dline(100), aline, line, line2
character(16) rtype 
character(80) d0line*80
character(32) word, function_def, var, full_var, arg(100)
character delim

integer i, i_arg, id, ix, ixa, ixp, ixc, ix_end, ix1, n_rline, ios
integer ixa2, ix_word, ix2, lun, n_contains

logical found(100), delim_found, have_found_1, has_contains, valid

! open the f90 file

valid = dir_open ('.')
if (.not. valid) stop

file_loop: do

  valid = dir_read (f90_file)
  if (.not. valid) stop
  if (index(f90_file, '.f90') /= len_trim(f90_file) - 3) cycle

  open (1, file = f90_file, status = 'old')
  print *, 'Opening: ', trim(f90_file)

  ! init init 

  have_found_1 = .false.  ! set True when the 1st routine block is found

  ! init

  do

    n_rline = 1                      
    i_arg = 0
    found = .false.
    has_contains = .false.
    n_contains = 0

    ! find the subroutine or funtion line 

    do 

      read (1, '(a)', iostat = ios) rline(1)
      if (ios /= 0) then
        if (.not. have_found_1) &
                    type *, 'ERROR: CANNOT FIND "SUBROUTINE" OR "FUNCTION" NAME'
        cycle file_loop        
      endif

      call string_trim (rline(1), rline(1), ix)
      call str_upcase (line, rline(1))
      i = index(line, '!')
      if (i /= 0) line = line(1:ix-1)

      if (line(:ix) == 'CONTAINS') then
        has_contains = .true.
        cycle
      endif

      if (line(:ix) == 'SUBROUTINE') then
        if (has_contains) then
          n_contains = n_contains + 1
          cycle
        else
          rtype = 'subroutine'
          exit
        endif
      endif

      function_def = line(:ix)
      call string_trim(line(ix+1:), line, ix)

      if (function_def == 'END' .and. &
            (line(:ix) == 'FUNCTION' .or. line(:ix) == 'SUBROUTINE')) then
        if (has_contains) then
          n_contains = n_contains - 1
          if (n_contains == -1) has_contains = .false.
        endif
        cycle
      endif

      if (line(:ix) == 'FUNCTION' .or. function_def == 'FUNCTION') then
        if (has_contains) then
          n_contains = n_contains + 1
          cycle
        endif
        rtype = 'function'
        i_arg = i_arg + 1
        if (line(:ix) == 'FUNCTION') then
          found(i_arg) = .true.
          call string_trim (line(ix+1:), line, ix)
        endif
        ix = index(line, '(')
        if (ix == 0) then
          arg(i_arg) = line
        else        
          arg(i_arg) = line(:ix-1)
        endif
        exit
      endif

    enddo

    ! get the argument list

    have_found_1 = .true.

    ix = index(line, '(')
    if (ix == 0) goto 2000

    call  string_trim (line(ix+1:), line, ix)
    if (line(1:1) == ')') goto 2000

    do
      ix = index(line, ',')
      if (ix /= 0) then
        call get_arg
      else
        ix = index(line, ')')
        if (ix /= 0) then
          call get_arg
          exit      
        else
          ix = index(line, '&')                                  
          if (ix /= 0) then
            n_rline = n_rline+1
            read (1, '(a)') rline(n_rline)
            line = line(:ix-1) // rline(n_rline)
            call str_upcase (line, line)
          else
            type *, 'ERROR: ROUTINE LINE DOES NOT END WITH ")" NOR "&"'
            call err_exit
          endif
        endif
      endif
    enddo

    2000 continue

    ! if a function then look for a "result" tag

    if (rtype == 'function') then
      ix = index (line, 'RESULT')
      if (ix /= 0) then
        call string_trim(line(ix+6:), line, ix)
        if (line(1:1) == '&') then
          n_rline = n_rline + 1
          read (1, '(a)', iostat = ios) rline(n_rline)
          call str_upcase (line, rline)
          call string_trim (line, line, ix)
        endif
        if (line(1:1) /= '(') then
          print *, 'ERROR: CANNOT FINE "(" AFTER "RESULT" FOR A FUNCTION!'
          call err_exit
        endif
        ix = index(line, ')')
        arg(1) = line(2:ix-1)
      endif
    endif

    !      type *, 'Args:', i_arg
    !      do i = 1, i_arg   
    !        type '(i4, 2x, a)', i, arg(i)
    !      enddo

    ! read lines to match the argument list with the defs
                                
    id = 0
                         
    do
      read (1, '(a)', iostat = ios) aline
      if (ios /= 0) then
        type *, 'ERROR: CANNOT FIND ALL THE SUBROUTINE ARGUMENTS. CANNOT FIND:'
        do i = 1, i_arg
          if (.not. found(i)) type '(i4, 2x, a)', i, arg(i)
        enddo
        call exit
      endif
      ix = index (aline, '!')
      if (ix /= 0) aline = aline(:ix-1)

      call string_trim (aline, aline, ix)
      if (ix == 0) cycle
      call str_upcase (line, aline)

      if (line(1:3) == 'USE' .or. line(1:7) == 'INCLUDE' &
                                    .or. line(1:8) == 'IMPLICIT') then
        id = id + 1
        dline(id) = aline

        if (line(1:3) == 'USE' .and. index(line, '_INTERFACE') /= 0) then
          print *, 'NOTE: *NOT* using "USE" statement in file: ', trim(aline)
          id = id - 1
        endif

      elseif (line(1:) == 'END SUBROUTINE') then
        type *, 'ERROR: CANNOT FIND ALL THE SUBROUTINE ARGUMENTS. CANNOT FIND:'
        do i = 1, i_arg
          if (.not. found(i)) type '(i4, 2x, a)', i, arg(i)
        enddo
        call exit

      else

        ix1 = index(line, '::')

        if (ix1 == 0) then
          if (line(1:4) == 'TYPE') then
            ix1 = index(line, ')')
          elseif (line(1:6) == 'RECORD') then
            ix1 = index(line, '/')
            ix1 = ix1 + index(line(ix1+1:), '/')
          else
            call string_trim (line, line, ix1)
          endif
        else
          ix1 = ix1 + 1
        endif

        d0line = aline(:ix1)
        aline = aline(ix1+1:)

        do
          call word_read (aline, ',(/*', word, ix_word, delim, delim_found, aline)
          if (delim == ',') then
            full_var = word
          elseif (delim == '(') then
            ix = index(aline, ')')
            full_var = trim(word) // '(' // aline(:ix)
            aline = aline(ix+1:)
          elseif (delim == '/') then
            full_var = word
            ix = index(aline, '/')
            aline = aline(ix+1:)
          elseif (delim == '*') then
            ix = index(aline, ',')
            if (ix == 0) then
              full_var = trim(word) // '*' // aline
              aline = ""
            else
              full_var = trim(word) // '*' // aline(:ix-1)
              aline = aline(ix+1:)
            endif
          else
            full_var = word
          endif
          call str_upcase(var, word)
          do i = 1, i_arg
            if (var == arg(i)) then
              id = id + 1
              dline(id) = trim(d0line) // ' ' // full_var
              found(i) = .true.
            endif
          enddo
          if (.not. delim_found) exit
        enddo
      endif
      if (all(found(1:i_arg))) exit
    enddo

    !  type *
    !  type *, 'Defs:', id
    !  do i = 1, id
    !    type *, trim(dline(i))
    !  enddo

    open (unit = 2, file = 'f90.interface', status = 'unknown', &
                             carriagecontrol = 'list', position = 'append')
    write (2, '(a)') 'interface'
    do i = 1, n_rline
      write (2, '(2x, a)') trim(rline(i))
    enddo
    do i = 1, id
      write (2, '(4x, a)') trim(dline(i))
    enddo
    write (2, '(2x, 2a)') 'end ', trim(rtype)
    write (2, '(a)') 'end interface'
    write (2, *)


    type *
    type '(a)', 'interface'
    do i = 1, n_rline
      type '(2x, a)', trim(rline(i))
    enddo
    do i = 1, id
      type '(4x, a)', trim(dline(i))
    enddo
    type '(2x, 2a)', 'end ', trim(rtype)
    type '(a)', 'end interface'
    type *
    type *, 'Interface appended to file: f90.interface'
    type *

  enddo

  close (1)

enddo file_loop

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

contains

subroutine get_arg ()
  i_arg = i_arg + 1
  call string_trim(line(:ix-1), arg(i_arg), ix2)
  line = line(ix+1:)
end subroutine

end program
