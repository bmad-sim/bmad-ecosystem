!+
! Subroutine transfer_lat (lat1, lat2)
!
! Subroutine to set lat2 = lat1. 
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, lat2's pointers point to the same memory as lat1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Input:
!   lat1 -- lat_struct:
!
! Output:
!   lat2 -- lat_struct:
!-

subroutine transfer_lat (lat1, lat2)

use bmad_struct

implicit none

type (lat_struct), intent(in) :: lat1
type (lat_struct), intent(out) :: lat2

lat2 = lat1

end subroutine transfer_lat

