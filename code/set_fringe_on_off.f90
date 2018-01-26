!+
! Subroutine set_fringe_on_off (fringe_at, ele_end, on_or_off) 
!
! Routine to modify a ele%value(fringe_at$) setting to either turn on or turn off a fringe
! at either the entrance or exit end of the element. 
!
! Input:
!   fringe_at       -- real(rp): Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
!   ele_end         -- integer: Element edge: entrance_end$ or exit_end$
!   on_or_off       -- integer: Turn on$ or off$
!
! Output:
!   fringe_at       -- real(rp): Modified fringe setting.
!- 

subroutine set_fringe_on_off (fringe_at, ele_end, on_or_off) 

use bmad_struct

implicit none

integer ele_end, on_or_off
real(rp) fringe_at

!

select case (on_or_off)

! On

case (on$)

  select case (ele_end)
  case (entrance_end$)
    select case (nint(fringe_at))
    case (no_end$);   fringe_at = entrance_end$
    case (exit_end$); fringe_at = both_ends$
    end select

  case (exit_end$)
    select case (nint(fringe_at))
    case (no_end$);       fringe_at = exit_end$
    case (entrance_end$); fringe_at = both_ends$
    end select

  case default
    call err_exit  ! Should not be here
  end select

! Off

case (off$)

  select case (ele_end)
  case (entrance_end$)
    select case (nint(fringe_at))
    case (both_ends$);    fringe_at = exit_end$
    case (entrance_end$); fringe_at = no_end$
    end select

  case (exit_end$)
    select case (nint(fringe_at))
    case (both_ends$); fringe_at = entrance_end$
    case (exit_end$);  fringe_at = no_end$
    end select

  case default
    call err_exit  ! Should not be here
  end select

case default
  call err_exit  ! Should not be here
end select


end subroutine set_fringe_on_off

