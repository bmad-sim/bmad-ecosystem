!+
! Subroutine type_layout (lay)
!
! Subroutine to print the global information in a layout
!
! Modules Needed:
!   use accelerator
!
! Input:
!   lay - layout: layout to use.
!+

subroutine type_layout (lay)

  use accelerator

  implicit none

  type (layout) lay

!

  if (.not. associated(lay%start)) then
    type *, 'Warning from TYPE_LAYOUT: Layout NOT Associated'
    return
  endif

  type *, 'Name:         ', lay%name
  type *, 'N:            ', lay%N,        '  ! Number of Elements'
  type *, 'LatPos:       ', lay%lastpos,  '  ! Last position'

end subroutine
