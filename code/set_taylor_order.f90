!+
! Subroutine set_taylor_order (order, override_flag)
!
! Subroutine to set the taylor order for the Taylor maps.
!
! Note: override_flag = .false. is generally only used by bmad_parser so that
! if the taylor order has been previously set then the setting in the 
! lattice file will not override it.
!
! Note: Calling this routine after calling bmad_parser will not reset any
! taylor maps made by bmad_parser. Thus when in doubt, call this routine
! before calling bmad_parser.
!
! Note: This routine does not call any of Etienne's PTC routines since this
! routine may be called before PTC has been initialized. This routine
! just sets a global variable and returns.
!
! Modules needed:
!   use bmad
!
! Input:
!   order         -- Integer: Taylor order.
!                     If order = 0. then nothing is done.
!   override_flag -- Logical, optional: If False then if the taylor order 
!                     has been previously set do not reset.
!-

subroutine set_taylor_order (order, override_flag)

  use bmad_struct

  implicit none

  integer, intent(in) :: order
  logical, optional, intent(in) :: override_flag
  logical override

! do nothing if order = 0

  if (order == 0) return

  if (order < 0 .or. order > 100) then
    print *, 'ERROR IN SET_TAYLOR_ORDER: ORDER OUT OF BOUNDS:', order
    call err_exit
  endif

! check for override_flag and do nothing if the taylor order has been set

  override = .true.
  if (present(override_flag)) override = override_flag
  if (.not. override .and. bmad_com%taylor_order_set) return

! set the taylor order.

  bmad_com%taylor_order = order
  bmad_com%taylor_order_set = .true.    

end subroutine
