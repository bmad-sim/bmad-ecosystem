!+
! Subroutine str_substitute (string, str_in, str_out)
!
! Routine to substitute all instances of one sub-string for another in a string
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   string  -- Character(*): Character string.
!   str_in  -- Character(*), optional: Sub-string to be replaced. 
!                 Default is the tab char.
!   str_out -- Character(*): optional: Sub-string to substitute in. 
!                 Default is the space char.
!
! Output:
!   string  -- Character(*): String with all instances of str_in replaced for str_out.
!-

subroutine str_substitute (string, str_in, str_out)

implicit none

character(*) string
character(*), optional :: str_in, str_out

integer ixs, ix, n_in

!

if (present(str_in)) n_in = len(str_in)
ixs = 1

do 

  if (ixs > len(str_in)) return

  if (present(str_in)) then
    ix = index(string(ixs:), str_in)
  else
    ix = index(string(ixs:), char(9))  ! tab
  endif

  if (ix == 0) return

  if (present(str_out)) then
    string = string(:ixs+ix-2) // str_out // trim(string(ixs+ix+n_in-1:))
    ixs = ixs + ix + n_in - 1
  else
    string = string(:ixs+ix-2) // ' ' // trim(string(ixs+ix:))
    ixs = ixs + ix 
  endif

enddo

end subroutine
