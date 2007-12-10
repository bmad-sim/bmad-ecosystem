!+
! Subroutine remove_ele_from_lat (lat, ix_ele)
!
! Subroutine to remove an element from the tracking part of the lattice.
! The element to be removed must have no controllers.
! If a controller is detected this routine will stop the program.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- Lat_struct: Lattice with element(s) to be removed.
!   ix_ele -- Integer, optional: Index of element in lat%ele(:) array to be removed.
!               If not present than all ele%key = null_ele$ elements will be removed.
!
! Output:
!   lat -- Lat_struct: Lattice with element(s) removed.
!-

subroutine remove_ele_from_lat (lat, ix_ele)

  use bmad_struct
  use bmad_interface, except_dummy => remove_ele_from_lat

  implicit none

  type (lat_struct) lat
  integer, optional :: ix_ele
  integer i, j, ix, n_remove
  real(rp) length
  character(20) :: r_name = 'remove_ele_from_lat'

! Mark elements to be removed

  if (present(ix_ele)) then
    lat%ele(ix_ele)%ix_ele = 0
    lat%ele(ix_ele+1:lat%n_ele_max)%ix_ele = lat%ele(ix_ele+1:lat%n_ele_max)%ix_ele - 1
    n_remove = 1
    
  else
    n_remove = 0
    do i = 1, lat%n_ele_track
      if (lat%ele(i)%key == null_ele$) then
        lat%ele(i)%ix_ele = 0
        n_remove = n_remove + 1
        cycle
      endif
      lat%ele(i)%ix_ele = lat%ele(i)%ix_ele - n_remove
    enddo
  endif

! Rectify the lord/slave bookkeeping.

  do i = 1, lat%n_control_max

    ix = lat%control(i)%ix_lord
    if (lat%ele(ix)%ix_ele == 0) then
      call out_io (s_fatal$, r_name, &
                        'ELEMENT TO BE REMOVED FROM LATTICE IS CONTROLLER:' // lat%ele(ix)%name)
      call err_exit
    endif
    lat%control(i)%ix_lord = lat%ele(ix)%ix_ele

    ix = lat%control(i)%ix_slave
    if (lat%ele(ix)%ix_ele == 0) then
      call out_io (s_fatal$, r_name, &
                        'ELEMENT TO BE REMOVED FROM LATTICE IS CONTROLLED:' // lat%ele(ix)%name)
      call err_exit
    endif
    lat%control(i)%ix_slave = lat%ele(ix)%ix_ele

  enddo

! Remove the elements

  do i = 1, lat%n_ele_max
    ix = lat%ele(i)%ix_ele 
    if (ix == i .or. ix == 0) cycle
    lat%ele(ix) = lat%ele(i)
  enddo

  do i = lat%n_ele_max-n_remove+1, lat%n_ele_max
    call init_ele (lat%ele(i))
  enddo

! More bookkeeping adjustment

  lat%n_ele_track = lat%n_ele_track - n_remove
  lat%n_ele_max = lat%n_ele_max - n_remove

  call s_calc(lat)
  call lat_geometry (lat)

end subroutine
