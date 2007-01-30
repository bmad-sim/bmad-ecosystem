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
! Note: If a taylor map does not exist for an element this routine will 
! make one and store it in the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lattice used in the calculation.
!   ix1   -- Integer, optional: Element start index for the calculation.
!              Default is 0.
!   ix2   -- Integer, optional: Element end index for the calculation.
!              Default is lat%n_ele_track.
!
! Output:
!    t_map(6) -- Taylor_struct: Transfer map.
!-

#include "CESR_platform.inc"

subroutine transfer_map_calc (lat, t_map, ix1, ix2)

  use bmad_struct
  use bmad_interface, except => transfer_map_calc
  use cesr_utils, only: integer_option
  use ptc_interface_mod, only: concat_taylor, ele_to_taylor

  implicit none

  type (lat_struct) lat

  type (taylor_struct) :: t_map(:), taylor2(6)

  integer, intent(in), optional :: ix1, ix2
  integer i, i1, i2

!

  i1 = integer_option(0, ix1) 
  i2 = integer_option(lat%n_ele_track, ix2) 

  if (associated(lat%ele(i1)%taylor(1)%term)) then
    t_map = lat%ele(i1)%taylor
  else
    call mat6_to_taylor (lat%ele(i1)%vec0, lat%ele(i1)%mat6, t_map)
  endif

  if (i2 < i1) then
    do i = i1+1, lat%n_ele_track
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

  if (.not. associated(lat%ele(i)%taylor(1)%term)) then
    call ele_to_taylor (lat%ele(i), lat%param)
  endif

  call concat_taylor (t_map, lat%ele(i)%taylor, t_map)

end subroutine

end subroutine
                                          
