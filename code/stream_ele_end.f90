!+
! Function stream_ele_end (physical_end, ele_orientation) result (stream_end)
!
! Rotine to determine which stream end of an element a particle is at given 
! the position in terms of the physical end and the element's orientation
!
! Input:
!   physical_end     -- Integer: entrance_end$, exit_end$, surface$, etc.
!   ele_orientation  -- Integer: Either 1 = Normal or -1 = element reversed.
!
! Output:
!   stream_end       -- Integer: upstream_end$, downstream_end$, or set equal
!                         to physical_end if physical_end is neither entrance_end$
!                         nor exit_end$
!-

function stream_ele_end (physical_end, ele_orientation) result (stream_end)

use bmad_struct

implicit none

integer stream_end, ele_orientation, physical_end
character(*), parameter :: r_name  = 'stream_ele_end'

!

if (physical_end /= entrance_end$ .and. physical_end /= exit_end$) then
  stream_end = physical_end
  return
endif

!

select case (ele_orientation)

case (1) 
  select case (physical_end)
  case (entrance_end$);   stream_end = upstream_end$
  case (exit_end$);       stream_end = downstream_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD PHYSICAL_END: \i0\ ', i_array = [physical_end])
    if (global_com%exit_on_error) call err_exit
  end select

case (-1) 
  select case (physical_end)
  case (entrance_end$);   stream_end = downstream_end$
  case (exit_end$);       stream_end = upstream_end$
  case default;
    call out_io (s_fatal$, r_name, 'BAD PHYSICAL_END: \i0\ ', i_array = [physical_end])
    if (global_com%exit_on_error) call err_exit
  end select

case default
  call out_io (s_fatal$, r_name, 'BAD ELEMENT ORIENTATION: \i0\ ', i_array = [ele_orientation])
  if (global_com%exit_on_error) call err_exit

end select

end function stream_ele_end 

