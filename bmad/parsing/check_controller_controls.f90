!+
! Subroutine check_controller_controls (ele_key, contrl, name, err)
!
! Routine to check for problems when setting up group, overlay, or ramper controllers.
!
! Input:
!   ele_key     -- integer: Element type. overlay$, etc.
!   contrl(:)   -- control_struct: control info. 1 element for each slave.
!   name        -- character(*): Lord name. Used for error reporting.
!
! Output:
!   err         -- logical: Set true if there is a problem. False otherwise.
!-

subroutine check_controller_controls (ele_key, contrl, name, err)

use bmad_parser_mod, dummy => check_controller_controls

implicit none

type (control_struct), target :: contrl(:)
type (control_struct), pointer :: c1, c2
integer ele_key, i, j
logical err
character(*) name

! Ramper elements do not have control%slave set so check control%slave_name.
! Ignore control%slave_name for other controller types since non-ramper elements may
! control multiple slaves with the same name.

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

