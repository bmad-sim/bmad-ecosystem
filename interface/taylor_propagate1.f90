!+
! Subroutine taylor_propagate1 (tlr, ele, param)
!
! Subroutine to track a real_8 taylor map through an element.
! The alternative routine if ele has a taylor series is concat_taylor.
!
! Modules needed:
!   use accelerator
!
! Input:
!   tlr(6) -- Taylor_struct: Map to be tracked
!   ele    -- Ele_struct: Element to track through
!   param  -- Param_struct: 
!
! Output:
!   tlr(6)  -- Taylor_struct: Map through element
!-

subroutine taylor_propagate1 (tlr, ele, param)

  use accelerator

  implicit none

  type (taylor_struct) tlr(:)
  type (real_8), save :: y(6)
  type (ele_struct) ele
  type (param_struct) param
  type (fibre), pointer, save :: a_fibre

  logical :: init_needed = .true.

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! init

  if (init_needed) then
    allocate (a_fibre)
    call real_8_init (y)
    init_needed = .false.
  endif

!

  y = tlr

  call alloc_fibre (a_fibre)
  call ele_to_fibre (ele, a_fibre, param)
  call track (a_fibre, y, default, +1)
  call kill (a_fibre)

  tlr = y

end subroutine
