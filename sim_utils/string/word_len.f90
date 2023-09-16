!+
! Function WORD_LEN (WORDING)
!
! Function to return the "length" of the first word in a string
! If string(1:1) = ' ' than WORD_LEN is set to 0
!
! Example:
!     ix = word_len ('to be or not    ')
! Result:
!     ix = 2
!-

integer function word_len (wording)

use precision_def

implicit none

integer i1, i
character(*) wording

!

i1 = len(wording)
do i = 1, i1
  if (wording(i:i) == ' ' .or. wording(i:i) == char(9)) then
    word_len = i - 1
    return
  endif
enddo

word_len = i1

end function
