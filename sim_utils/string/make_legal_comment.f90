!+
! Subroutine MAKE_LEGAL_COMMENT (COMMENT_IN, COMMENT_OUT)
!
! Subroutine to double up on "'" characters in a comment line and
! add "'" characters on either side.
! Useful for namelist read.
!
! Input:
!   COMMENT_IN -- Character*(*): Input comment line
!
! Output:
!   COMMENT_OUT -- Character*(*): output comment line
!-

subroutine make_legal_comment (comment_in, comment_out)

use precision_def

implicit none

character*(*) comment_in, comment_out
character*130 line

integer ixc, ixo, ixq, ix

logical quote_pending

!

ixc = 1
ixo = 1
line = comment_in

quote_pending = .true.
do while (quote_pending)
  ixq = index(line(ixc:), "'")
  if (ixq == 0) then
    comment_out(ixo:) = line(ixc:)
    quote_pending = .false.
  else
    comment_out(ixo:) = line(ixc:ixc+ixq-1) // "'"
    ixc = ixc + ixq
    ixo = ixo + ixq + 1
  endif
enddo

ix = len_trim(comment_out)
if (ix == 0) then
  comment_out = "' '"
elseif (ix .le. len(comment_out)) then
  comment_out = "'" // trim(comment_out) // "'"
else
  comment_out = "'" // comment_out(:len(comment_out)-2) // "'"
endif

end subroutine
