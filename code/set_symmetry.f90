!+
! Subroutine set_symmetry (symmetry, ring)
!
! Subroutine to set the symmetry of a ring.
! This is not used much now.
!
! Modules Needed:
!   use bmad
!
! Input:
!   symmetry -- Integer: symmetry to set to.
!
! Output:
!   ring     -- Ring_struct: ring%param%symmetry and ring%n_ele_use get
!                  set accordingly
!-

#include "CESR_platform.inc"

subroutine set_symmetry (symmetry, ring)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring

  integer symmetry

!

  ring%param%symmetry = symmetry

  if (symmetry == no_symmetry$) then
    ring%n_ele_use = ring%n_ele_ring
  elseif (symmetry == ew_antisymmetry$) then
    ring%n_ele_use = ring%n_ele_symm
  elseif (symmetry == mobius_symmetry$) then
    ring%n_ele_use = ring%n_ele_ring
  else
    print *, 'ERROR IN SET_SYMMETRY: UNKNOWN SYMMETRY:', symmetry
    call err_exit
  endif

  return
  end
