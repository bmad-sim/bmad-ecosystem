!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine convert_blanks_to_underscore (string_in, string_out)

  use sr_struct
  use sr_interface

  implicit none

  integer n
  character(*) string_in, string_out

!

  string_out = string_in
  do n = 1, len_trim(string_out)
    if (string_out(n:n) == ' ') string_out(n:n) = '_'
  enddo

end subroutine
