!+
! Subroutine find_element_ends (lat, ele, ele1, ele2)
!
! Subroutine to find the end points of an element in the tracking part of the 
! lattice.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat  -- lat_struct: Lat holding the lattice
!   ele  -- Ele_struct: Element to find the ends for.
!
! Output:
!   ele1 -- Ele_struct, pointer:  Pointer to element just before ele.
!   ele2 -- Ele_struct, pointer:  Pointer to ele itself or the last sub-element within ele.
!
! Note: For an element in the regular part of the lat:
!       ele1%ix_ele = ele%ix_ele - 1
!       ele2        => ele
!-

#include "CESR_platform.inc"

subroutine find_element_ends (lat, ele, ele1, ele2)

  use bmad_struct
  use bmad_interface, except_dummy => find_element_ends
  use nr, only: indexx

  implicit none
                                                         
  type (lat_struct) lat
  type (ele_struct), target :: ele
  type (ele_struct), pointer :: ele1, ele2

  integer ix_ele, ix_start, ix_end, ix_start_line, ix_end_line
  integer ix1, ix2, n, n_end, n_slave, ix_slave, ix_line
  integer, allocatable, save :: ix_slave_array(:), ix_line_array(:)

!

  ix_ele = ele%ix_ele
  ix1 = ele%ix1_slave
  ix2 = ele%ix2_slave

  if (ele%n_slave == 0) then
    call pointer_to_element (lat, ele%ix_photon_line, ix_ele-1, ele1)
    ele2 => ele

  elseif (ele%control_type == super_lord$) then
    call pointer_to_element (lat, ele%ix_photon_line, lat%control(ix1)%ix_slave - 1, ele1)
    call pointer_to_element (lat, lat%control(ix2)%ix_slave, ele2)

  else  ! overlay_lord$ or group_lord$
    ix_start = lat%n_ele_track + 1
    ix_end = 0
    n = 0
    n_slave = ele%n_slave
    call re_allocate(ix_slave_array, n_slave)
    call re_allocate(ix_line_array, n_slave)
    ix_slave_array(1:n_slave) = lat%control(ix1:ix2)%ix_slave
    ix_line_array(1:n_slave) = lat%control(ix1:ix2)%ix_photon_line
    n_end = n_slave
    do 
      n = n + 1
      if (n > n_end) return
      ix_slave = ix_slave_array(n)
      ix_line = ix_line_array(n)
      if (ix_slave > lat%n_ele_track .and. ix_line == 0) then
        n_slave = lat%ele(ix_slave)%n_slave
        ix1 = lat%ele(ix_slave)%ix1_slave
        ix2 = lat%ele(ix_slave)%ix2_slave
        call re_allocate(ix_slave_array, n_slave+n_end)
        call re_allocate(ix_line_array, n_slave+n_end)
        ix_slave_array(n_end+1:n_end+n_slave) = lat%control(ix1:ix2)%ix_slave
        ix_line_array(n_end+1:n_end+n_slave) = lat%control(ix1:ix2)%ix_photon_line
        n_end = n_end + n_slave
      else
        if (ix_slave - 1 < ix_start) then
          ix_start = ix_slave - 1
          ix_start_line = ix_line
        endif
        if (ix_start > ix_end) then
          ix_end = ix_slave 
          ix_end_line = ix_line
        endif
      endif
    enddo

    call pointer_to_element (lat, ix_start, ix_start_line, ele1)
    call pointer_to_element (lat, ix_end, ix_end_line, ele2)

  endif

end subroutine
