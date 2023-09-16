!+
! Subroutine string_trim (in_string, out_string, word_len)
!
! Subroutine to trim a string of leading blanks and/or tabs and also to 
! return the length of the first word.
!
! Input:
!  in_string - Character(*): String constant to be parsed.
!
! Output:
!   out_string - Character(*): String constant: Left shifted to trim leading blanks
!                      and/or tabs.
!   word_len   - Integer: Set to the length of the first word. 
!                    If string has no non-blank characters
!                    the value of word_len is set to 0.
!
! Example:
!     line = '  tobe or   nuts '
!     call string_trim(line, line, ix)
!
! Result:
!     line => 'tobe or  nuts '
!     ix   => 4
!
! The first word is then: line(1:ix)
!
! Additional words may be obtained by repeated calls to string_trim using 
! out_string(IX+1:) as the in_string for the next call
! 
! Example: 
!     line = '  tobe or   nuts '
!     call string_trim(line, line, ix)
!     call string_trim(line(ix+1:), line, ix)
!
! Result:
!     line => 'or  nuts '
!     ix   => 2
!-

subroutine string_trim (in_string, out_string, word_len)

implicit none

character(*) in_string, out_string
character(len(in_string)) temp_str   ! In case in_str and out_str are the same actual arg
character(1) tab

integer i, j, word_len, len_in, len_out

parameter (tab = char(9))

! find number of leading blenks

len_in = len(in_string)
len_out = len(out_string)

do i = 1, len_in
  if (in_string(i:i) .ne. ' ' .and. in_string(i:i) .ne. tab) exit
enddo            

! If IN_STRING is entirely blanks and/or tab characters then ...

if (i == len_in+1) then
  word_len = 0
  out_string = ' '
  return
endif

! Left shift in_string and put in out_string. 

j = min(len_in, len_out + i - 1)
temp_str = in_string(i:j)
out_string = temp_str

! Count characters in first word

do i = 1, len_out
  if (out_string(i:i) .eq. ' ' .or. out_string(i:i) .eq. tab) then
    word_len = i - 1
    return
  endif
enddo

! here if no blanks or tabs

word_len = len_out 

end subroutine




