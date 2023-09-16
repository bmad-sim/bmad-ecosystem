!+
! Function at_this_ele_end (now_at, where_at) result (is_at_this_end)
!
! Routine to determine if an aperture or fringe field is present.
! When tracking use this routine conjunction with physical_ele_end.
!
! Input:
!   now_at      -- Integer: Which end is under consideration: entrance_end$, exit_end$, surface$, or in_between$.
!   where_at    -- Integer: Which ends have the aperture or fringe field: entrance_end$, exit_end$, 
!                     continuous$, both_ends$, no_aperture$, surface$, wall_transition$.
!
! Output:
!   is_at_this_end   -- Logical: True if at this end. False otherwise.
!- 

function at_this_ele_end (now_at, where_at) result (is_at_this_end)

use bmad_interface, dummy => at_this_ele_end

implicit none

integer now_at, where_at
logical is_at_this_end

!

if (where_at == no_aperture$) then
  is_at_this_end = .false.
  return
endif

if (now_at == surface$ .or. where_at == surface$) then
  is_at_this_end = (now_at == where_at)
  return
endif

if (where_at == continuous$ .or. where_at == wall_transition$) then
  is_at_this_end = .true.
  return
endif

!

select case (now_at)
case (entrance_end$)
  select case (where_at)
  case (entrance_end$, both_ends$); is_at_this_end = .true.
  case default;                     is_at_this_end = .false.
  end select

case (exit_end$)
  select case (where_at)
  case (exit_end$, both_ends$); is_at_this_end = .true.
  case default;                 is_at_this_end = .false.
  end select

case (in_between$)
  is_at_this_end = .false.
end select

end function at_this_ele_end

