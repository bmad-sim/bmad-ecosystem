!+
! Subroutine find_element_ends (lat, ix_ele, ix_start, ix_end)
!
! Subroutine to find the end points of an element in the regular part of the 
! lat.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lat holding the lattice
!   ix_ele -- Integer: Index of element to find the ends for.
!
! Output:
!   ix_start -- Integer: Index of the element just before the element.
!   ix_end   -- Integer: Index of element itself or index of the
!                   last sub-element within the element.
!
! Note: For an element in the regular part of the lat:
!       ix_start = ix_ele - 1
!       ix_end = ix_ele
!-

#include "CESR_platform.inc"

subroutine find_element_ends (lat, ix_ele, ix_start, ix_end)

  use bmad_struct
  use bmad_interface, except => find_element_ends
  use nr, only: indexx

  implicit none
                                                         
  type (lat_struct) lat

  integer ix_ele, ix_start, ix_end, ix_slave
  integer ix1, ix2, n, n_end, n_slave
  integer, allocatable, save :: ixs(:)

!

  ix1 = lat%ele(ix_ele)%ix1_slave
  ix2 = lat%ele(ix_ele)%ix2_slave

  if (lat%ele(ix_ele)%n_slave == 0) then
    ix_start = ix_ele - 1
    ix_end = ix_ele

  elseif (lat%ele(ix_ele)%control_type == super_lord$) then
    ix_start = lat%control(ix1)%ix_slave - 1
    ix_end = lat%control(ix2)%ix_slave

  else  ! overlay_lord$ or group_lord$
    ix_start = lat%n_ele_track + 1
    ix_end   = 0
    n = 0
    n_slave = lat%ele(ix_ele)%n_slave
    call re_allocate(ixs, n_slave)
    ixs(1:n_slave) = lat%control(ix1:ix2)%ix_slave
    n_end = n_slave
    do 
      n = n + 1
      if (n > n_end) return
      ix_slave = ixs(n)
      if (ix_slave > lat%n_ele_track) then
        n_slave = lat%ele(ix_slave)%n_slave
        ix1 = lat%ele(ix_slave)%ix1_slave
        ix2 = lat%ele(ix_slave)%ix2_slave
        call re_allocate(ixs, n_slave+n_end)
        ixs(n_end+1:n_end+n_slave) = lat%control(ix1:ix2)%ix_slave
        n_end = n_end + n_slave
      else
        ix_start = min(ix_start, ix_slave-1)
        ix_end   = max(ix_end, ix_slave)
      endif
    enddo

  endif

end subroutine
