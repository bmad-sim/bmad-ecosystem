!+
! Subroutine check_controller_controls (contrl, name, err)
!
! Routine to check for problems when setting up group or overlay controllers.
!
! Input:
!   contrl(:)   -- Control_struct: control info. 1 element for each slave.
!   name        -- Character(*): Lord name. Used for error reporting.
!
! Output:
!   err         -- Logical: Set true if there is a problem. False otherwise.
!-

subroutine check_controller_controls (contrl, name, err)

use bmad_parser_mod, dummy => check_controller_controls

implicit none

type (control_struct), target :: contrl(:)
type (control_struct), pointer :: c1, c2
integer i, j
logical err
character(*) name

!

err = .true.

do i = 1, size(contrl)
  c1 => contrl(i)
  do j = i+1, size(contrl)
    c2 => contrl(j)
    if (c1%slave == c2%slave .and. c1%attribute == c2%attribute) then
      call parser_error ('DUPLICATE SLAVE CONTROL FOR LORD: ' // name)
      return
    endif
  enddo
enddo

err = .false.

end subroutine check_controller_controls

