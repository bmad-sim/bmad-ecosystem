!+
! Function tao_universe_number (i_uni) result (i_this_uni)
!
! Fnction to return the universe number.
! i_this_uni = i_uni except when i_uni is -1. 
! In this case i_this_uni = s%com%default_universe.
!
! Input:
!   i_uni -- Integer: Nominal universe number.
!
! Output:
!   i_this_uni -- Integer: Universe number. 
!-

function tao_universe_number (i_uni) result (i_this_uni)

use tao_struct

implicit none

integer i_uni, i_this_uni

i_this_uni = i_uni
if (i_uni == -1) i_this_uni = s%com%default_universe

end function tao_universe_number

