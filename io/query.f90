!+
! Subroutine : query_real( parameter, default, fmt)
!-

subroutine query_real (parameter, result, fmt)

  use precision_def

  implicit none
  character(*) parameter, fmt
  character(10) string
  real(rp) result
  
  !

  write (*, '(1x, 2a,'// fmt // ', a)', advance = 'NO') parameter, ' (default = ', result, ') ?'
  read (*, '(a)') string
  if (string == '' ) return
  read (string, *) result

end subroutine

!---------------------------------------------------------------------

subroutine query_int( parameter, result, fmt)

  implicit none
  character(*) parameter, fmt
  character(10) string
  integer result
  
  !

  write (*, '(1x, 2a,'// fmt // ', a)', advance = 'NO') parameter, ' (default = ', result, ') ?'
  read (*, '(a)') string
  if (string == '' ) return
  read (string, *) result

end subroutine

!---------------------------------------------------------------------

subroutine query_character( parameter, result, fmt)

  character(*) parameter, fmt, result
  character(72) string
  
  !

  write (*, '(1x, 2a,'// fmt // ', a)', advance = 'NO') parameter, ' (default = ', result, ') ?'
  read (*, '(a)') string
  if (string == '' ) return
  result = string

end subroutine
