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

  integer ix_ele, ix_start, ix_end, ix_start_branch, ix_end_branch
  integer ix1, ix2, n, n_end, n_slave, ix_slave, ix_branch
  integer, allocatable, save :: ix_slave_array(:), ix_branch_array(:)

!

  ix_ele = ele%ix_ele
  ix1 = ele%ix1_slave
  ix2 = ele%ix2_slave

  if (ele%n_slave == 0) then
    call pointer_to_ele (lat, ele%ix_branch, ix_ele-1, ele1)
    ele2 => ele

  elseif (ele%lord_status == super_lord$) then
    call pointer_to_ele (lat, ele%ix_branch, lat%control(ix1)%ix_slave - 1, ele1)
    call pointer_to_ele (lat, ele%ix_branch, lat%control(ix2)%ix_slave, ele2)

  ! For overlays and groups: The idea is to look at all the slave elements in the tracking 
  ! part of the lattice and find the minimum and maximum element indexes.
  ! The small complication is that overlays or groups lords can control other overlays or 
  ! groups, etc.
  ! So we must "recursively" follow the slave tree.
  ! ix_slave_array/ix_branch_array holds the list of slaves we need to look at.

  else  ! overlay_lord$ or group_lord$
    ix_start = lat%n_ele_track + 1
    ix_end = 0
    n = 0       ! Index in ix_slave_array
    n_slave = ele%n_slave
    call re_allocate(ix_slave_array, n_slave)
    call re_allocate(ix_branch_array, n_slave)
    ix_slave_array(1:n_slave) = lat%control(ix1:ix2)%ix_slave
    ix_branch_array(1:n_slave) = lat%control(ix1:ix2)%ix_branch
    n_end = n_slave
    do 
      n = n + 1
      if (n > n_end) exit
      ix_slave = ix_slave_array(n)
      ix_branch = ix_branch_array(n)
      ! If the slave itself has slaves then add the sub-slaves to the list
      if (ix_slave > lat%n_ele_track .and. ix_branch == 0) then
        n_slave = lat%ele(ix_slave)%n_slave
        ix1 = lat%ele(ix_slave)%ix1_slave
        ix2 = lat%ele(ix_slave)%ix2_slave
        call re_allocate(ix_slave_array, n_slave+n_end)
        call re_allocate(ix_branch_array, n_slave+n_end)
        ix_slave_array(n_end+1:n_end+n_slave) = lat%control(ix1:ix2)%ix_slave
        ix_branch_array(n_end+1:n_end+n_slave) = lat%control(ix1:ix2)%ix_branch
        n_end = n_end + n_slave
      ! Else this slave is in the tracking part of the lattice...
      else
        if (ix_slave - 1 < ix_start) then
          ix_start = ix_slave - 1
          ix_start_branch = ix_branch
        endif
        if (ix_slave > ix_end) then
          ix_end = ix_slave 
          ix_end_branch = ix_branch
        endif
      endif
    enddo

    call pointer_to_ele (lat, ix_start_branch, ix_start, ele1)
    call pointer_to_ele (lat, ix_end_branch, ix_end, ele2)

  endif

end subroutine
