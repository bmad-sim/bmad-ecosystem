!  REAL_READ   SUBROUTINE  PROGRAMING  C.DCS.LIB   DCS         96.7
!+
! function REAL_READ(error_message)
!
! function equivalent to TYPSCN function RELTYP with added automatic error
! message generated
!
! example:
!
!     real(rp) = real_read('ERROR: ')
!-

function real_read(error_message) result (r_read)


  use precision_def

  implicit none

  integer idelim

  real(rp) reltyp, r_read

  logical lwait

  character(*) error_message

!

  lwait = .false.

  r_read = reltyp(lwait, idelim)

  if (idelim > 0) then
    print *
    print *, error_message
    call zertyp(' ')
    stop
  endif

  return

end function
