!+
! ***********************************************************************
! ***********************************************************************
! ***************                                         ***************
! ***************    OBSOLETE: USE SET_ON_OFF INSTEAD     ***************
! ***************                                         ***************
! ***********************************************************************
! ***********************************************************************
!                        
! Subroutine set_on (key, lat, on_switch, orb)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a lat. An element that is turned off acts like a drift.
! lat_make_mat6 will be called to remake lat%ele()%mat6
!
! Modules needed:
!   use bmad
!
! Input:
!   key       -- Integer: Key name of elements to be turned on or off.
!                  [Key = quadrupole$, etc.]
!   lat      -- lat_struct: Lat structure holding the elements
!   on_switch -- Logical: True  => turn elements on.
!                         False => turn elements off.
!   orb(0:)  -- Coord_struct, optional: Needed for lat_make_mat6
!
! Output:
!   lat -- lat_struct: Modified lat.
!-

#include "CESR_platform.inc"
                                    
subroutine set_on (key, lat, on_switch, orb)

  use bmad_struct
  use bmad_interface, except_dummy => set_on

  implicit none

  type (lat_struct) lat
  type (coord_struct), optional :: orb(0:)

  integer i, key               

  logical on_switch

!

  do i = 1, lat%n_ele_max

    if (lat%ele(i)%key == key) then
      lat%ele(i)%is_on = on_switch
      call lat_make_mat6(lat, i, orb)
    endif

  enddo

end subroutine
