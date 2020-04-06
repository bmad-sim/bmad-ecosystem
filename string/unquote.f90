!+
! Function unquote (str_in) result (str_out)
!
! Routine to remove quotation marks at the ends of a string.
! Quotation marks will only be removed if they match at both ends.
! If no removal, str_out = str_in
!
! Input:
!   str_in -- Character(*): Input string
!
! Output:
!   str_out -- Character(*): Output string
!
! Example:
!   unquote ('"This"') -> 'This'
!   unquote ('Interior"not"removed') -> 'Interior"not"removed'
!-                               

function unquote (str_in) result (str_out)

implicit none

character(*) str_in
character(len(str_in)) str_out
character(1) quote

integer i, i1

!

do i = 1, len(str_in)
  select case (str_in(i:i))
  case (' ') 
    cycle
  case ('"', "'")
    quote = str_in(i:i)
    i1 = i
    exit
  case default
    str_out = str_in
    return
  end select
enddo

if (i == len(str_in) + 1) then
  str_out = str_in
  return
endif

do i = len(str_in), 1, -1
  select case (str_in(i:i))
  case (' ') 
    cycle
  case ('"', "'")
    if (str_in(i:i) == quote) then
      str_out = str_in(1:i1-1) // str_in(i1+1:i-1)
    else
      str_out = str_in
    endif
    return
  case default
    str_out = str_in
    return
  end select
enddo

str_out = str_in

end function
