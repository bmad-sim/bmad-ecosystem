!+         
! Subroutine transfer_map_calc (lat, t_map, ix1, ix2)
!
! Subroutine to calculate the transfer map between two elements.
!
! The transfer map is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn map is calculated.
! If ix2 < ix1 then the calculation will "wrap around" the lattice end.
! For example if ix1 = 900 and ix2 = 10 then the t_mat is the map from
! element 900 to the lattice end plus from 0 through 10.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- Ring_struct: Lattice used in the calculation.
!   ix1   -- Integer, optional: Element start index for the calculation.
!              Default is 0.
!   ix2   -- Integer, optional: Element end index for the calculation.
!              Default is lat%n_ele_use.
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

#include "CESR_platform.inc"

subroutine transfer_map_calc (lat, t_map, ix1, ix2)

  use bmad_struct
  use bmad_interface, except => transfer_map_calc
  use cesr_utils, only: integer_option

  implicit none

  type (ring_struct) lat

  type (taylor_struct) :: t_map(:), taylor2(6)

  integer, intent(in), optional :: ix1, ix2
  integer i, i1, i2

!

  i1 = integer_option(ix1, 0) 
  i2 = integer_option(ix2, lat%n_ele_use) 

  if (associated(lat%ele_(i1)%taylor(1)%term)) then
    t_map = lat%ele_(i1)%taylor
  else
    call mat6_to_taylor (lat%ele_(i1)%mat6, lat%ele_(i1)%vec0, t_map)
  endif

  if (i2 < i1) then
    do i = i1+1, lat%n_ele_use
      call add_on_to_t_map
    enddo
    do i = 1, i2
      call add_on_to_t_map
    enddo

  else
    do i = i1+1, i2
      call add_on_to_t_map
    enddo
  endif

!--------------------------------------------------------
contains

subroutine add_on_to_t_map

  if (associated(lat%ele_(i)%taylor(1)%term)) then
    call concat_taylor (t_map, lat%ele_(i)%taylor, t_map)
  else
    call mat6_to_taylor (lat%ele_(i)%mat6, lat%ele_(i)%vec0, taylor2)
    call concat_taylor (t_map, taylor2, t_map)
  endif

end subroutine

end subroutine
                                          
