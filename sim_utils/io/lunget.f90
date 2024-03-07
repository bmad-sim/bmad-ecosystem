!+
! Function lunget()
!
! Function to return a free file unit number to be used with an open statement.
!
! Output:
!   lunget -- Integer: Free file number between 30 and 60
!
! Example:
!   iu = lunget()
!   open (iu, file = ... )
!-

integer function lunget()

use output_mod, dummy => lunget

implicit none

logical op
integer i
character(*), parameter :: r_name = 'lunget'


! Note: If you change the range [30,60] also err_exit needs to be modified.

lunget=0                ! happy compiler

do i = 30, 60
  inquire(unit = i, opened = op)
  if(.not. op) then
    lunget = i
    return
  endif
enddo

!

call out_io (s_fatal$, r_name, 'NO LOGICAL UNIT OPENABLE BETWEEN 30 AND 60!')
if (global_com%exit_on_error) call err_exit
return

end function



