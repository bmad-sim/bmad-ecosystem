!+
! Subroutine s_calc (lat)
!
! Subroutine to calculate the longitudinal distance S for the elements
! in a lat.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat -- lat_struct:
!
! Output:
!   lat -- lat_struct:
!-

#include "CESR_platform.inc"

subroutine s_calc (lat)

  use bmad_struct
  use bmad_interface, except => s_calc
  
  implicit none

  type (lat_struct)  lat

  integer n, ix2
  real*8 ss

! Just go through all the elements and add up the lengths.

  ss = lat%ele(0)%s

  do n = 1, lat%n_ele_track
    ss = lat%ele(n-1)%s + lat%ele(n)%value(l$)
    lat%ele(n)%s = ss
  enddo

  lat%param%total_length = ss - lat%ele(0)%s

! now get fill in the positions of the super_lords

  do n = lat%n_ele_track+1, lat%n_ele_max
    if (lat%ele(n)%control_type == super_lord$) then
      ix2 = lat%control(lat%ele(n)%ix2_slave)%ix_slave
      lat%ele(n)%s = lat%ele(ix2)%s
    else
      lat%ele(n)%s = 0
    endif
  enddo

end subroutine
