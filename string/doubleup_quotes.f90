!+
! Subroutine doubleup_quotes (str_in, str_out, quote)
!
! Subroutine to doubleup on quotation marks in a string
!
! Input:
!   str_in -- Character(*): Input string
!   quote  -- Character(1): quotation mark to double up (" or ')
!
! Output:
!   str_out -- Character(*): Output string
!
! Example:
!   call doubleup_quotes ('This:":', str_out, '"')
!   str_out -> 'This:"":'
!-                               

subroutine doubleup_quotes (str_in, str_out, quote)

implicit none

character(*) str_in, str_out
character(1) quote
character(len(str_out)) str

integer ix_in, ix_out, ixq

!

ix_in = 1
ix_out = 1

do
  ixq = index(str_in(ix_in:), quote)
  if (ixq == 0) then
    str(ix_out:) = str_in(ix_in:)
    if (len(str_out) .lt. len_trim(str)) then
      print *, 'WARNING FROM DOUBLEUP_QUOTES: OUTPUT STRING IS TOO SHORT!'
    endif
    str_out = str
    return
  else
    str(ix_out:) = str_in(ix_in:ix_in+ixq-1) // quote
    ix_in = ix_in + ixq
    ix_out = ix_out + ixq + 1
  endif
enddo

end subroutine
