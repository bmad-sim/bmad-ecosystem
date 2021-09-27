!+
! Function time_direction() result (time_sign)
!
! Routine to return +1 if bmad_com%backwards_time_tracking_on = False and -1 otherwise.
!
! Output:
!   time_sign   -- integer: +1 or -1.
!-

function time_direction() result (time_sign)

use bmad_struct, only: bmad_com
implicit none
integer time_sign

!

select case (bmad_com%backwards_time_tracking_on)
case (.true.)
  time_sign = -1
case default
  time_sign = +1
end select

end function time_direction
