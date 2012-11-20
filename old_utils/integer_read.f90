!+
! function INTEGER_READ(error_message)
!
! function equivalent to TYPSCN function RELTYP with added automatic error
! message generated
!
! example:
!
!     integer = integer_read('ERROR: ')
!-

integer function integer_read(error_message)

use precision_def

implicit none

integer idelim, inttyp

logical lwait

character*(*) error_message

!

lwait = .false.

integer_read = inttyp(lwait, idelim)

if (idelim > 0) then
  print *
  print *, error_message
  call zertyp(' ')
  stop
endif

end function
