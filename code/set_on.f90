!+
! ***********************************************************************
! ***********************************************************************
! ***************                                         ***************
! ***************    OBSOLETE: USE SET_ON_OFF INSTEAD     ***************
! ***************                                         ***************
! ***********************************************************************
! ***********************************************************************
!                        
! Subroutine set_on (key, ring, on_switch, orb_)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a ring. An element that is turned off acts like a drift.
! RING_MAKE_MAT6 will be called to remake ring%ele_()%mat6
!
! Modules needed:
!   use bmad
!
! Input:
!   key       -- Integer: Key name of elements to be turned on or off.
!                  [Key = quadrupole$, etc.]
!   ring      -- Ring_struct: Ring structure holding the elements
!   on_switch -- Logical: True  => turn elements on.
!                         False => turn elements off.
!   orb_(0:)  -- Coord_struct, optional: Needed for ring_make_mat6
!
! Output:
!   ring -- Ring_struct: Modified ring.
!-

#include "CESR_platform.inc"
                                    
subroutine set_on (key, ring, on_switch, orb_)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct), optional :: orb_(0:)

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
