!+
! Subroutine elements_locator_by_key (key, lat, indx)
!
! NOTE: THIS ROUTINE IS DEPRECATED SINCE IT DOES NOT WORK WITH BRANCHES.
! USE LAT_ELE_LOCATOR INSTEAD.
!
! Subroutine to locate all the elements of a certain kind in a lat. 
! Note: super_slave elements are not included in the list since super_slaves
! cannot have their attributes changed directly.
!
! Modules Needed:
!   use bmad
!
! Input:
!   key  -- Integer: Element key. For example quadrupole$, etc.
!   lat -- lat_struct: Lat to search through.
!
! Output:
!   indx(:) -- Integer, pointer: Array of indexes of the elements found.
!              Note: This routine does not try to deallocate indx.
!               It is up to you to deallocate indx if needed.
!-

subroutine elements_locator_by_key (key, lat, indx)

  use bmad_struct
  use bmad_interface, except_dummy => elements_locator_by_key
  
  implicit none

  type (lat_struct) lat

  integer key
  integer, pointer :: indx(:)
  integer i, ix

!

  ix = count(lat%ele(:)%key == key .and. &
                                  lat%ele(:)%slave_status /= super_slave$)

  allocate(indx(ix))

  ix = 0
  do i = 0, lat%n_ele_max
    if (lat%ele(i)%key == key .and. &
                              lat%ele(i)%slave_status /= super_slave$) then
      ix = ix + 1
      indx(ix) = i
    endif
  enddo

end subroutine
