!+
! Subroutine SET_SYMMETRY (SYMMETRY, RING)
!
! Subroutine to set the symmetry of a ring.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     SYMMETRY  -- Integer: symmetry to set to
!
! Output:
!     RING      -- Ring_struct: RING.PARAM.SYMMETRY and RING.N_ELE_USE get
!                  set accordingly
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine set_symmetry (symmetry, ring)

  use bmad_struct
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
    type *, 'ERROR IN SET_SYMMETRY: UNKNOWN SYMMETRY:', symmetry
    call err_exit
  endif

  return
  end
