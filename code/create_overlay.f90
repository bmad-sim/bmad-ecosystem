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
!     %ix_branch     -- Index of branch element belongs to.
!     %ix_attrib     -- Index of attribute controlled
!     %coef          -- Coefficient
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!   err_print_flag -- Logical, optional: If present and False then suppress
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

use bmad_interface, except_dummy => create_overlay
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: slave, lord
type (control_struct)  contl(:)

integer i, j, nc0, ix_overlay, nc2, ix_con
integer ix_slave, n_slave, ix_attrib, ix_branch

character(*) attrib_name
character(40) at_name
character(16), parameter :: r_name = 'create_overlay'
logical err, free
logical, optional :: err_print_flag

! Error check

n_slave = size (contl)

do j = 1, n_slave
  ix_slave  = contl(j)%ix_slave
  ix_branch = contl(j)%ix_branch

  if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) then
    call out_io (s_fatal$, r_name,  'BRANCH INDEX OUT OF BOUNDS. \i0\ ', ix_branch)
    if (global_com%exit_on_error) call err_exit
  endif

  if (ix_slave <= 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    call out_io (s_fatal$, r_name,  'INDEX OUT OF BOUNDS. \i0\ ', ix_slave)
    if (global_com%exit_on_error) call err_exit
  endif
enddo

! Mark element as an overlay lord

lord => lat%ele(ix_overlay)
call check_controller_controls (contl, lord%name, err)
if (err) return

lord%lord_status = overlay_lord$
lord%key = overlay$
call set_ele_defaults(lord)

call str_upcase (at_name, attrib_name)
ix_attrib =  attribute_index (lord, at_name)
if (ix_attrib == 0) then
  call out_io (s_fatal$, r_name,  'BAD ATTRIBUTE_NAME: ' // attrib_name, &
                                  'TRYING TO CREATE OVERLAY: ' // lord%name)
  if (global_com%exit_on_error) call err_exit
endif
lord%component_name = at_name
lord%ix_value = ix_attrib

if (n_slave == 0) return ! If no slaves then nothing to do.

! Loop over all slaves.

nc0 = lat%n_control_max
nc2 = nc0

do j = 1, n_slave
  ix_attrib = contl(j)%ix_attrib
  ix_slave = contl(j)%ix_slave
  slave => lat%branch(contl(j)%ix_branch)%ele(ix_slave)

  if (nc2+4 > size(lat%control)) call reallocate_control (lat, nc2+100)

  ! If the slave attribute is a multipole component, make sure it exists.
  if (ix_attrib > num_ele_attrib$ .and. .not. associated (slave%a_pole)) then
    call multipole_init(slave)
  endif
  free = attribute_free (slave, attribute_name(slave, ix_attrib), err_print_flag, .true.)
  err = err .or. .not. free
  lat%control(nc2+1) = contl(j)
  lat%control(nc2+1)%ix_lord = ix_overlay
  nc2 = nc2 + 1

enddo

lord%n_slave = n_slave
lord%ix1_slave = nc0 + 1
lord%ix2_slave = nc0 + n_slave
lat%n_control_max = nc2

! Loop over all slaves
! Free elements convert to overlay slaves.

do i = 1, lord%n_slave

  slave => pointer_to_slave(lord, i, ix_con)

  if (slave%slave_status == free$) slave%slave_status = control_slave$

  ! You cannot overlay super_slaves 

  if (slave%slave_status == super_slave$) then
    call out_io (s_fatal$, r_name,  'ILLEGAL OVERLAY ON ' // slave%name, &
                                    ' BY: ' // lord%name)
    if (global_com%exit_on_error) call err_exit
  endif

  ! update controller info for the slave ele

  slave%n_lord = slave%n_lord + 1
  call add_lattice_control_structs (lat, slave)
  lat%ic(slave%ic2_lord) = ix_con

enddo

! Finish: Do control bookkeeping.

call control_bookkeeper (lat, lord)

end subroutine


