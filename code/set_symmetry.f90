!+
! Subroutine SET_SYMMETRY (SYMMETRY, RING)
!
! Subroutine to set the symmetry of a ring.
!
! Modules Needed:
!   use bmad
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
!Revision 1.5  2003/05/02 15:44:02  dcs
!F90 standard conforming changes.
!
!Revision 1.4  2003/01/27 14:40:43  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:24  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

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
