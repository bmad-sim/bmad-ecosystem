!+
! Subroutine create_overlay (ring, ix_overlay, attrib_name, contl)
!
! Subroutine to add the controller information to slave elements of
! an overlay_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring        -- Ring_struct: Ring to modify.
!   ix_overlay  -- Integer: Index of overlay element.
!   attrib_name -- Character(16): Name of attribute in the overlay that is
!                     to be varied.
!   contl(:)      -- Control_struct: control info. 1 element for each slave.
!     %ix_slave   -- Index of element to control
!     %ix_attrib  -- Index of attribute controlled
!     %coef       -- Coefficient
!
! Output:
!   ring    -- Ring_struct: Modified ring.
!
! Note: Use NEW_CONTROL to get an index for the overlay element
!
! Example:
!   call new_control (ring, ix_ovr)      ! get index of overlay in ring%ele_
!   ring%ele_(ix_ovr)%name = 'OVERLAY1'  ! overlay name
!   ring%ele_(ix_ovr)%value(k1$) = 0.1   ! starting value
!
!   contl(1)%ix_slave = 10   ! RING%ELE_(10) is, say, a quadrupole.
!   contl(1)%ix_attrib = k1$ ! The overlay controls the quadrupole strength.
!   contl(1)%coef = 0.1      ! A change in the overlay value of 1 produces
!                            !    a change of 0.1 in k1 of element 10.
!
!   contl(2)%ix_slave = 790  ! RING%ELE_(790) is, say, a sextupole.
!   contl(2)%ix_attrib = k2$ ! The overlay controls the sextupole strength.
!   contl(2)%coef = -0.1     ! make changes antisymmetric.
!
!   call create_overlay (ring, ix_ovr, 'K1', contl(1:2))  ! create the overlay
!-

#include "CESR_platform.inc"

subroutine create_overlay (ring, ix_overlay, attrib_name, contl)

  use bmad_struct
  use bmad_interface, except => create_overlay
  use bookkeeper_mod, only: control_bookkeeper

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: slave, lord
  type (control_struct)  contl(:)

  integer i, j, ix, ixc, ix_overlay, n_con2
  integer ix_slave, n_slave, ix_attrib, slave_type

  character(*) attrib_name
  character(16) at_name

! Mark element as an overlay lord

  lord => ring%ele_(ix_overlay)
  n_slave = size (contl)

  ix = ring%n_control_max
  n_con2 = ix + n_slave
  if (n_con2 > size(ring%control_)) &
                      ring%control_ => reallocate (ring%control_, n_con2+500)

  do j = 1, n_slave
    ring%control_(ix+j) = contl(j)
    ring%control_(ix+j)%ix_lord = ix_overlay
  enddo

  lord%n_slave = n_slave
  lord%ix1_slave = ix + 1
  lord%ix2_slave = ix + n_slave
  lord%control_type = overlay_lord$
  lord%key = overlay$
  ring%n_control_max = n_con2

  call str_upcase (at_name, attrib_name)
  ix_attrib =  attribute_index (lord, at_name)
  if (ix_attrib == 0) then
    print *, 'ERROR IN CREATE_OVERLAY: BAD ATTRIBUTE_NAME: ', attrib_name
    print *, '      TRYING TO CREATE OVERLAY: ', lord%name
    call err_exit
  endif
  lord%attribute_name = at_name
  lord%ix_value = ix_attrib

! Loop over all slaves
! Free elements convert to overlay slaves.

  do i = lord%ix1_slave, lord%ix2_slave

    ix_slave = ring%control_(i)%ix_slave
    if (ix_slave <= 0) then
      print *, 'ERROR IN CREATE_OVERLAY: INDEX OUT OF BOUNDS.', ix_slave
      call err_exit
    endif

    slave => ring%ele_(ix_slave)
    slave_type = slave%control_type

    if (slave_type == free$) slave%control_type = overlay_slave$

! You cannot overlay super_slaves 

    if (slave_type == super_slave$) then
      print *, 'ERROR IN CREATE_OVERLAY: ILLEGAL OVERLAY ON ', slave%name
      print *, '      BY: ', lord%name
      call err_exit
    endif

! update controller info for the slave ele

    slave%n_lord = slave%n_lord + 1
    call adjust_control_struct (ring, ix_slave)
    ixc = slave%ic2_lord
    ring%ic_(ixc) = i

  enddo

  call control_bookkeeper (ring, ix_overlay)

end subroutine


