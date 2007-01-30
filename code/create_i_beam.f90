!+
! Subroutine create_i_beam (lat, ix_i_beam, ix_slave, ele_init)
!
! Subroutine to add the controller information to slave elements of
! an i_beam_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- lat_struct: Lat to modify.
!   ix_i_beam    -- Integer: Index of i_beam element.
!   ix_slave(:)  -- Index of element to control
!   ele_init     -- Element containing attributes to be transfered
!                   to the I_Beam element:
!                       ele_init%name        
!                       ele_init%alias
!                       ele_init%descrip
!                       ele_init%value(:)
!
! Output:
!   lat    -- lat_struct: Modified lat.
!
! Note: Use NEW_CONTROL to get an index for the i_beam element
!
! Example: Create an I_Beam supporting elements 
! lat%ele(10) and lat%ele(12)
!
!   call new_control (lat, ix_ele)        ! get IX_ELE index
!   call create_i_beam (lat, ix_ele, (/ 10, 12 /))  ! create the i_beam
!-

#include "CESR_platform.inc"

subroutine create_i_beam (lat, ix_i_beam, ix_slave, ele_init)

  use bmad_struct
  use bmad_interface, except => create_i_beam

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), optional :: ele_init
  type (ele_struct), pointer ::  slave, i_beam

  integer, intent(in) :: ix_i_beam, ix_slave(:)
  integer i, j, ix, ix2, ixc, n_con2
  integer ixs, slave_type, idel, n_slave

  real(rdef) s_max, s_min

! Mark element as an i_beam lord

  i_beam => lat%ele(ix_i_beam)

  n_slave = size (ix_slave)
  ix = lat%n_control_max
  n_con2 = ix + n_slave

  if (n_con2 > size(lat%control)) &
                      lat%control => reallocate (lat%control, n_con2+500)

  do j = 1, n_slave
    lat%control(ix+j)%ix_slave  = ix_slave(j)
    lat%control(ix+j)%ix_lord   = ix_i_beam
    lat%control(ix+j)%coef      = 0
    lat%control(ix+j)%ix_attrib = 0
  enddo

  i_beam%n_slave = n_slave
  i_beam%ix1_slave = ix + 1
  i_beam%ix2_slave = ix + n_slave
  i_beam%control_type = i_beam_lord$
  i_beam%key = i_beam$
  lat%n_control_max = n_con2

! Loop over all slaves
! Free elements convert to overlay slaves.

  s_max = -1e30  ! something large and negative
  s_min =  1e30  ! something large and positive

  do i = i_beam%ix1_slave, i_beam%ix2_slave

    ixs = lat%control(i)%ix_slave
    if (ixs <= 0) then
      print *, 'ERROR IN CREATE_I_BEAM: INDEX OUT OF BOUNDS.', ixs
      call err_exit
    endif

    slave => lat%ele(ixs)
    slave_type = slave%control_type

    if (slave_type == free$) slave%control_type = overlay_slave$

! You cannot control super_slaves, group_lords or overlay_lords

    if (slave_type == super_slave$ .or. slave_type == group_lord$ .or. &
                                            slave_type == overlay_lord$) then
      print *, 'ERROR IN CREATE_I_BEAM: ILLEGAL I_BEAM ON ', slave%name
      print *, '      BY: ', i_beam%name
      call err_exit
    endif

! update controller info for the slave ele

    slave%n_lord = slave%n_lord + 1
    call add_lattice_control_structs (lat, ixs)
    ixc = slave%ic2_lord
    lat%ic(ixc) = i

! compute min/max

    s_max = max(s_max, slave%s)
    s_min = min(s_min, slave%s-slave%value(l$))

  enddo

! center of i_beam

  i_beam%value(s_center$) = (s_max + s_min) / 2

! ele_init stuff

  if (present(ele_init)) then
    i_beam%name    = ele_init%name
    i_beam%alias   = ele_init%alias
    i_beam%value   = ele_init%value
    if (associated(ele_init%descrip)) then
      allocate (i_beam%descrip)
      i_beam%descrip = ele_init%descrip
    endif
  endif

end subroutine
