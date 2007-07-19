!+
! Subroutine compress_lat (LAT, OK)
!
! Subroutine to compress the ele(*) and control(*) arrays to remove
! elements no longer used. Note: to mark an element for removal use:
!     lat%ele(i)%key = -1
!
! Modules Needed:
!   use bmad
!
! Input:
!     LAT -- lat_struct: Lat to compress.
!
! Output:
!     LAT -- lat_struct: Compressed lat.
!     OK   -- Logical: Lat compressed OK.
!-

#include "CESR_platform.inc"

subroutine compress_lat (lat, ok)

  use bmad_struct
  use bmad_interface, except_dummy => compress_lat

  implicit none
                           
  type (lat_struct), target :: lat
  type (ele_struct), pointer :: ele

  type (control_struct), allocatable :: control(:)

  integer i, j, ix, i2
  integer n_ic, n_con, n_lord, n
  integer, allocatable :: ixa(:), ic(:)

  logical ok

! allocate

  n = lat%n_control_max
  allocate (control(n))
  allocate (ixa(n))
  allocate (ic(lat%n_ic_max))

! remove unwanted lat%ele() elements

  ok = .true.

  i2 = 0
  do i = 1, lat%n_ele_max
    if (lat%ele(i)%key == -1) then
      ixa(i) = int_garbage$
    else
      i2 = i2 + 1
      ixa(i) = i2
      lat%ele(i2) = lat%ele(i)
    endif
    if (i == lat%n_ele_track) then
       lat%n_ele_track = i2
       lat%n_ele_track = i2
    endif
  enddo

  lat%n_ele_max = i2

! renumber lat%control()%ix_ele  

  forall (i = 1:lat%n_control_max) 
    lat%control(i)%ix_lord = ixa(lat%control(i)%ix_lord)
    lat%control(i)%ix_slave = ixa(lat%control(i)%ix_slave)
  end forall

! compress lat%control() array

  n_con = 0
  ix = 0

  do i = 1, lat%n_ele_max
    ele => lat%ele(i)
    do j = 1, ele%n_slave
      control(n_con+j) = lat%control(ele%ix1_slave+j-1)
      ixa(ele%ix1_slave+j-1) = n_con+j
    enddo
    if (ele%n_slave > 0) then
      ele%ix1_slave = n_con + 1
      ele%ix2_slave = n_con + ele%n_slave
      n_con = n_con + ele%n_slave
    endif
  enddo

! compress lat%ic() array

  n_ic = 0

  do i = 1, lat%n_ele_max
    ele => lat%ele(i)
    n_lord = 0
    do i2 = ele%ic1_lord, ele%ic2_lord
      ix = lat%ic(i2)
      if (ixa(ix) == 0) then
        if (ele%control_type == super_slave$) then
          print *, 'ERROR IN compress_lat: SUPERPOSITION LORD HAS BEEN REMOVED!'
          ok = .false.
        endif
      else
        n_lord = n_lord + 1                                           
        n_ic = n_ic + 1
        ic(n_ic) = ixa(ix)
      endif
    enddo

    if (n_lord == 0) then
      ele%ic1_lord = 0
      ele%ic2_lord = -1
      ele%n_lord = 0
      if (i <= lat%n_ele_track) ele%control_type = free$
    else
      ele%ic1_lord = n_ic - n_lord + 1
      ele%ic2_lord = n_ic
      ele%n_lord = n_lord
    endif

  enddo

  lat%control(1:n_con) = control(1:n_con)                     
  lat%n_control_max = n_con
  lat%ic(1:n_ic) = ic(1:n_ic)
  lat%n_ic_max = n_ic

! deallocate and do a check

  deallocate (control)
  deallocate (ixa)
  deallocate (ic)

  call check_lat_controls (lat, .true.)

end subroutine
          
