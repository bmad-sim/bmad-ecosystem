!+
! Subroutine NEW_CONTROL (RING, IX_ELE)
!
! Subroutine to create a new control element.
!
! Modules Needed:
!   use bmad
!
! Input:
!     RING -- Ring_struct: Ring used
!
! Output
!     IX_ELE -- Integer: Index of the new control element
!-

!$Id$
!$Log$
!Revision 1.4  2003/01/27 14:40:41  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:21  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:55  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine new_control (ring, ix_ele)

  use bmad_struct
  use bmad_interface

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
