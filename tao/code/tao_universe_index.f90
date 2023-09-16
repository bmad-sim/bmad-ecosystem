!+
! Function tao_universe_index (i_uni, neg2_to_default) result (i_this_uni)
!
! Fnction to return the universe number.
! Generally i_this_uni = i_uni except:
!   i_this_uni = -1  -> i_this_uni = s%global%default_universe.
!   i_this_uni = -2  -> i_this_uni = -2                      (neg2_to_default = F)
!                    -> i_this_uni = s%global%default_universe  (neg2_to_default = T)
!
! Input:
!   i_uni           -- integer: Nominal universe number.
!   neg2_to_default -- logical, optional: i_uni = -2 (all universes) maps to the default uni?
!                         Default if False.
! Output:
!   i_this_uni      -- integer: Universe number. 
!-

function tao_universe_index (i_uni, neg2_to_default) result (i_this_uni)

use tao_struct

implicit none

integer i_uni, i_this_uni
logical, optional :: neg2_to_default

!

select case (i_uni)
case (-2)
  if (logic_option(.false., neg2_to_default)) then
    i_this_uni = s%global%default_universe
  else
    i_this_uni = i_uni
  endif

case (-1)
  i_this_uni = s%global%default_universe

case default
  i_this_uni = i_uni
end select

end function tao_universe_index

