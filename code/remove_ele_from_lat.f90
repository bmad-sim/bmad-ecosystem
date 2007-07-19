!+
! Subroutine remove_ele_from_lat (lat, ix_ele)
!
! Subroutine to remove an element from the lattice.
! This subroutine assumes that the element has no controllers.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- Lat_struct: Lattice with element to be removed.
!   ix_ele -- Integer: Index of element in lat%ele(:) array.
!
! Output:
!   lat -- Lat_struct: Lattice with element removed.
!-

subroutine remove_ele_from_lat (lat, ix_ele)

  use bmad_struct
  use bmad_interface, except_dummy => remove_ele_from_lat

  implicit none

  type (lat_struct) lat
  integer ix_ele
  integer i, j
  real(rp) length
  character(20) :: r_name = 'remove_ele_from_lat'

! init

  length = lat%ele(ix_ele)%value(l$)

! Rectify the lord/slave bookkeeping.

  do i = 1, lat%n_control_max
    if (lat%control(i)%ix_lord == ix_ele .or. lat%control(i)%ix_slave == ix_ele) then
      call out_io (s_fatal$, r_name, &
                        'ELEMENT TO BE REMOVED FROM LATTICE IS CONTROLLED:' // lat%ele(ix_ele)%name)
      call err_exit
    endif
    if (lat%control(i)%ix_lord > ix_ele) lat%control(i)%ix_lord = lat%control(i)%ix_lord - 1
    if (lat%control(i)%ix_slave > ix_ele) lat%control(i)%ix_slave = lat%control(i)%ix_slave - 1
  enddo

! Remove the element.

  do i = ix_ele+1, lat%n_ele_max
    lat%ele(i-1) = lat%ele(i)
  enddo

  call init_ele (lat%ele(lat%n_ele_max))

! More bookkeeping adjustment

  if (ix_ele <= lat%n_ele_track) lat%n_ele_track = lat%n_ele_track - 1
  lat%n_ele_max = lat%n_ele_max - 1

  if (length /= 0) then
    call s_calc(lat)
    call lat_geometry (lat)
  endif

end subroutine
