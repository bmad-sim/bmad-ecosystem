!+
! Subroutine NEW_CONTROL (RING, IX_ELE)
!
! Subroutine to create a new control element.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING -- Ring_struct: Ring used
!
! Output
!     IX_ELE -- Integer: Index of the new control element
!-


subroutine new_control (ring, ix_ele)

  use bmad_struct
  implicit none
  type (ring_struct)  ring
  integer ix_ele

!

  ring%n_ele_max = ring%n_ele_max + 1
  ix_ele = ring%n_ele_max

  if (ix_ele > n_ele_maxx) then
    type *, 'ERROR IN NEW_CONTROL: NOT ENOUGH RING ELEMENTS!!!'
    type *, '      YOU NEED TO INCREASE N_ELE_MAXX IN BMAD_STRUCT!!!'
    call exit
  endif

  call init_ele (ring%ele_(ix_ele))

  return
  end
