!+
! Subroutine set_on (key, ring, on_switch, orb_)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a ring. An element that is turned off acts like a drift.
! RING_MAKE_MAT6 will be called to remake ring%ele_()%mat6
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   key       -- Integer: Key name of elements to be turned on or off.
!                  [Key = quadrupole$, etc.]
!   ring      -- Ring_struct: Ring structure holding the elements
!   on_switch -- Logical: True  => turn elements on.
!                         False => turn elements off.
!   orb_(0:n_ele_maxx) -- [Optional] Coord_struct: Needed for ring_make_mat6
!                         if ring%param%matrix_order = 1
!
! Output:
!   ring -- Ring_struct: Modified ring.
!-

!$Id$
!$Log$
!Revision 1.3  2002/01/08 21:44:43  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

                                    
subroutine set_on (key, ring, on_switch, orb_)

  use bmad_struct
  use bmad_interface, only: ring_make_mat6

  implicit none

  type (ring_struct) ring
  type (coord_struct), optional :: orb_(0:*)

  integer i, key               

  logical on_switch

!

  do i = 1, ring%n_ele_max

    if (ring%ele_(i)%key == key) then
      ring%ele_(i)%is_on = on_switch
      call ring_make_mat6(ring, i, orb_)
    endif

  enddo

end subroutine
