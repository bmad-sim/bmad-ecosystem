!+
! Subroutine create_i_beam (ring, ix_i_beam, ix_slave_)
!
! Subroutine to add the controller information to slave elements of
! an i_beam_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring         -- Ring_struct: Ring to modify.
!   ix_i_beam    -- Integer: Index of i_beam element.
!   ix_slave_(:) -- Index of element to control
!
! Output:
!   ring    -- Ring_struct: Modified ring.
!
! Note: Use NEW_CONTROL to get an index for the i_beam element
!
! Example:
!   call new_control (ring, ix_ele)        ! get IX_ELE index
!   ring%ele_(ix_ele)%name = 'I_BEAM1'     ! i_beam name
!   n_control = 2                          ! control 2 elements
!
!   ix_slave_(1) = 10   ! RING%ELE_(10) is Q01W say.
!   ix_slave_(2) = 12   ! RING%ELE_(12) is Q02W say.
!
!   call create_i_beam (ring, ix_ele, ix_slave_(1:2))  ! create the i_beam
!-

#include "CESR_platform.inc"

subroutine create_i_beam (ring, ix_i_beam, ix_slave_)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer ::  slave, i_beam

  integer, intent(in) :: ix_i_beam, ix_slave_(:)
  integer i, j, ix, ix2, ixc, n_con2
  integer ix_slave, slave_type, idel, n_slave

  real(rdef) s_max, s_min

! Mark element as an i_beam lord

  i_beam => ring%ele_(ix_i_beam)

  n_slave = size (ix_slave_)
  ix = ring%n_control_max
  n_con2 = ix + n_slave

  if (n_con2 > size(ring%control_)) &
                      ring%control_ => reallocate (ring%control_, n_con2+500)

  do j = 1, n_slave
    ring%control_(ix+j)%ix_slave  = ix_slave_(j)
    ring%control_(ix+j)%ix_lord   = ix_i_beam
    ring%control_(ix+j)%coef      = 0
    ring%control_(ix+j)%ix_attrib = 0
  enddo

  i_beam%n_slave = n_slave
  i_beam%ix1_slave = ix + 1
  i_beam%ix2_slave = ix + n_slave
  i_beam%control_type = i_beam_lord$
  i_beam%key = i_beam$
  ring%n_control_max = n_con2

! Loop over all slaves
! Free elements convert to overlay slaves.

  s_max = -1e30  ! something large and negative
  s_min =  1e30  ! something large and positive

  do i = i_beam%ix1_slave, i_beam%ix2_slave

    ix_slave = ring%control_(i)%ix_slave
    if (ix_slave <= 0) then
      print *, 'ERROR IN CREATE_I_BEAM: INDEX OUT OF BOUNDS.', ix_slave
      call err_exit
    endif

    slave => ring%ele_(ix_slave)
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
    call adjust_control_struct (ring, ix_slave)
    ixc = slave%ic2_lord
    ring%ic_(ixc) = i

! compute min/max

    s_max = max(s_max, slave%s)
    s_min = min(s_min, slave%s-slave%value(l$))

  enddo

! center of i_beam

  i_beam%value(s_center$) = (s_max + s_min) / 2

end subroutine
