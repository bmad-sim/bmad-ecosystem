!+
! Subroutine create_girder (lat, ix_girder, ix_slave, ele_init)
!
! Subroutine to add the controller information to slave elements of
! an girder_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- lat_struct: Lat to modify.
!   ix_girder    -- Integer: Index of girder element.
!   ix_slave(:)  -- Index of element to control
!   ele_init     -- Element containing attributes to be transfered
!                   to the Girder element:
!                       ele_init%name        
!                       ele_init%alias
!                       ele_init%descrip
!                       ele_init%value(:)
!
! Output:
!   lat    -- lat_struct: Modified lat.
!
! Note: Use NEW_CONTROL to get an index for the girder element
!
! Example: Create the Girder supporting elements 
! lat%ele(10) and lat%ele(12)
!
!   call new_control (lat, ix_ele)        ! get IX_ELE index
!   call create_girder (lat, ix_ele, (/ 10, 12 /))  ! create the girder
!-

subroutine create_girder (lat, ix_girder, ix_slave, ele_init)

  use bmad_struct
  use bmad_interface, except_dummy => create_girder

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), optional :: ele_init
  type (ele_struct), pointer ::  slave, girder

  integer, intent(in) :: ix_girder, ix_slave(:)
  integer i, j, ix, ix2, ixc, n_con2
  integer ixs, idel, n_slave

  real(rp) s_max, s_min

! Mark element as an girder lord

  girder => lat%ele(ix_girder)

  n_slave = size (ix_slave)
  ix = lat%n_control_max
  n_con2 = ix + n_slave

  if (n_con2 > size(lat%control)) call reallocate_control (lat, n_con2+500)

  do j = 1, n_slave
    lat%control(ix+j)%ix_slave  = ix_slave(j)
    lat%control(ix+j)%ix_lord   = ix_girder
    lat%control(ix+j)%coef      = 0
    lat%control(ix+j)%ix_attrib = 0
  enddo

  girder%n_slave = n_slave
  girder%ix1_slave = ix + 1
  girder%ix2_slave = ix + n_slave
  girder%lord_status = girder_lord$
  girder%key = girder$
  lat%n_control_max = n_con2

! Loop over all slaves
! Free elements convert to overlay slaves.

  s_max = -1e30  ! something large and negative
  s_min =  1e30  ! something large and positive

  do i = girder%ix1_slave, girder%ix2_slave

    ixs = lat%control(i)%ix_slave
    if (ixs <= 0) then
      print *, 'ERROR IN CREATE_GIRDER: INDEX OUT OF BOUNDS.', ixs
      call err_exit
    endif

    slave => lat%ele(ixs)

    if (slave%slave_status == free$ .or. slave%slave_status == group_slave$) &
                                                  slave%slave_status = overlay_slave$

! You cannot control super_slaves, group_lords or overlay_lords

    if (slave%slave_status == super_slave$ .or. slave%lord_status == group_lord$ .or. &
                                            slave%lord_status == overlay_lord$) then
      print *, 'ERROR IN CREATE_GIRDER: ILLEGAL GIRDER ON ', slave%name
      print *, '      BY: ', girder%name
      call err_exit
    endif

! update controller info for the slave ele

    slave%n_lord = slave%n_lord + 1
    call add_lattice_control_structs (lat, slave)
    ixc = slave%ic2_lord
    lat%ic(ixc) = i

! compute min/max

    s_max = max(s_max, slave%s)
    s_min = min(s_min, slave%s-slave%value(l$))

  enddo

! center of girder

  girder%value(s_center$) = (s_max + s_min) / 2

! ele_init stuff

  if (present(ele_init)) then
    girder%name    = ele_init%name
    girder%alias   = ele_init%alias
    girder%value   = ele_init%value
    if (associated(ele_init%descrip)) then
      allocate (girder%descrip)
      girder%descrip = ele_init%descrip
    endif
  endif

end subroutine
