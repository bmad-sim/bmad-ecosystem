!+
! Subroutine create_overlay (lat, ix_overlay, attrib_name, contl, err, err_print_flag)
!
! Subroutine to add the controller information to slave elements of
! an overlay_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: Lat to modify.
!   ix_overlay     -- Integer: Index of overlay element.
!   attrib_name    -- Character(40): Name of attribute in the overlay that is
!                     to be varied.
!   contl(:)       -- Control_struct: control info. 1 element for each slave.
!     %ix_slave      -- Index of element to control
!     %ix_attrib     -- Index of attribute controlled
!     %coef          -- Coefficient
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!   err_print_flag -- Logical, optional: If present and False then supress                                
!                       printing of an error message if attribute is not free.  
!
! Output:
!   lat    -- lat_struct: Modified lat.
!
! Note: Use NEW_CONTROL to get an index for the overlay element
!
! Example:
!   call new_control (lat, ix_ovr)     ! get index of overlay in lat%ele
!   lat%ele(ix_ovr)%name = 'OVERLAY1'  ! overlay name
!   lat%ele(ix_ovr)%value(k1$) = 0.1   ! starting value
!
!   contl(1)%ix_slave = 10   ! LAT%ele(10) is, say, a quadrupole.
!   contl(1)%ix_attrib = k1$ ! The overlay controls the quadrupole strength.
!   contl(1)%coef = 0.1      ! A change in the overlay value of 1 produces
!                            !    a change of 0.1 in k1 of element 10.
!
!   contl(2)%ix_slave = 790  ! LAT%ele(790) is, say, a sextupole.
!   contl(2)%ix_attrib = k2$ ! The overlay controls the sextupole strength.
!   contl(2)%coef = -0.1     ! make changes antisymmetric.
!
!   call create_overlay (lat, ix_ovr, 'K1', contl(1:2))  ! create the overlay
!-

subroutine create_overlay (lat, ix_overlay, attrib_name, contl, err, err_print_flag)

use bmad_struct
use bmad_interface, except_dummy => create_overlay
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: slave, lord
type (control_struct)  contl(:)

integer i, j, nc0, ix_overlay, nc2
integer ix_slave, n_slave, ix_attrib, ix_branch

character(*) attrib_name
character(40) at_name

logical err, free
logical, optional :: err_print_flag

! Mark element as an overlay lord

lord => lat%ele(ix_overlay)
call check_controller_controls (contl, lord%name, err)
if (err) return

n_slave = size (contl)

nc0 = lat%n_control_max
nc2 = nc0

do j = 1, n_slave
  ix_attrib = contl(j)%ix_attrib
  ix_slave = contl(j)%ix_slave
  slave => lat%branch(contl(j)%ix_branch)%ele(ix_slave)

  if (nc2+4 > size(lat%control)) call reallocate_control (lat, nc2+100)

  if (ix_attrib == x_limit$) then
    lat%control(nc2+1:nc2+2) = contl(j)
    lat%control(nc2+1:nc2+2)%ix_lord = ix_overlay
    lat%control(nc2+1)%ix_attrib = x1_limit$
    lat%control(nc2+2)%ix_attrib = x2_limit$
  elseif (ix_attrib == y_limit$) then
    lat%control(nc2+1:nc2+2) = contl(j)
    lat%control(nc2+1:nc2+2)%ix_lord = ix_overlay
    lat%control(nc2+1)%ix_attrib = y1_limit$
    lat%control(nc2+2)%ix_attrib = y2_limit$
  elseif (ix_attrib == aperture$) then
    lat%control(nc2+1:nc2+4) = contl(j)
    lat%control(nc2+1:nc2+4)%ix_lord = ix_overlay
    lat%control(nc2+1)%ix_attrib = x1_limit$
    lat%control(nc2+2)%ix_attrib = x2_limit$
    lat%control(nc2+3)%ix_attrib = y1_limit$
    lat%control(nc2+4)%ix_attrib = y2_limit$
    nc2 = nc2 + 4
  else
    ! If the slave attribute is a multipole component, make sure it exists.
    if (ix_attrib > n_attrib_maxx .and. .not. associated (slave%a_pole)) then
      call multipole_init(slave)
    endif
    free = attribute_free (slave, attribute_name(slave, ix_attrib), lat, err_print_flag, .true.)
    err = err .or. .not. free
    lat%control(nc2+1) = contl(j)
    lat%control(nc2+1)%ix_lord = ix_overlay
    nc2 = nc2 + 1
  endif
enddo

lord%n_slave = n_slave
lord%ix1_slave = nc0 + 1
lord%ix2_slave = nc0 + n_slave
lord%lord_status = overlay_lord$
lord%key = overlay$
lat%n_control_max = nc2

call str_upcase (at_name, attrib_name)
ix_attrib =  attribute_index (lord, at_name)
if (ix_attrib == 0) then
  print *, 'ERROR IN CREATE_OVERLAY: BAD ATTRIBUTE_NAME: ', attrib_name
  print *, '      TRYING TO CREATE OVERLAY: ', lord%name
  call err_exit
endif
lord%component_name = at_name
lord%ix_value = ix_attrib

! Loop over all slaves
! Free elements convert to overlay slaves.

do i = lord%ix1_slave, lord%ix2_slave

  ix_slave = lat%control(i)%ix_slave
  ix_branch = lat%control(i)%ix_branch

  if (ix_slave <= 0) then
    print *, 'ERROR IN CREATE_OVERLAY: INDEX OUT OF BOUNDS.', ix_slave
    call err_exit
  endif

  slave => lat%branch(ix_branch)%ele(ix_slave)

  if (slave%slave_status == free$ .or. slave%slave_status == group_slave$) &
                                                  slave%slave_status = overlay_slave$

  ! You cannot overlay super_slaves 

  if (slave%slave_status == super_slave$) then
    print *, 'ERROR IN CREATE_OVERLAY: ILLEGAL OVERLAY ON ', slave%name
    print *, '      BY: ', lord%name
    call err_exit
  endif

  ! update controller info for the slave ele

  slave%n_lord = slave%n_lord + 1
  call add_lattice_control_structs (lat, slave)
  lat%ic(slave%ic2_lord) = i

enddo

call control_bookkeeper (lat, ix_overlay)

end subroutine


