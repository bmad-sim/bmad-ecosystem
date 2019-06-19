!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_d2_d1_name (d1, show_universe) result (d2_d1_name)
!
! Function to return the datum name in the form:
!   d2_name.d1_name
! If there is only one d1_data array associated with the d2_data 
! array then the name is shortened to:
!   d2_name
! Additionally, if show_universe is True and there is more than one universe
! then "universe@" is prepended to the name.
! For example:
!   2@orbit.x
!   
!
! Input:
!   d1            -- Tao_d1_data_struct: Data array.
!   show_universe -- Logical, optional: Show the datum's universe.
!                       Default is True.
!
! Output:
!   d2_d1_name -- Character(60): Appropriate name.
!-

function tao_d2_d1_name(d1, show_universe) result (d2_d1_name)

use tao_struct

implicit none

type (tao_d1_data_struct) d1
character(60) d2_d1_name, temp_str
logical, optional :: show_universe

! If there is only one d1 array associated with the d2_data array then
! drop the d1 name.

if (size(d1%d2%d1) == 1) then
  write (d2_d1_name, '(a)') trim(d1%d2%name)
else
  write (d2_d1_name, '(3a)') trim(d1%d2%name), '.', trim(d1%name)
endif

!

if (size(s%u) > 1 .and. logic_option(.true., show_universe)) then
  ! Stupid gfortran compiler requires a temp string for this
  temp_str = d2_d1_name
  write (d2_d1_name, '(i0, 2a)') d1%d2%ix_universe, '@', trim(temp_str)
endif

end function tao_d2_d1_name
