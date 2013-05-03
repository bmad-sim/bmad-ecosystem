!+
! Subroutine s_calc (lat)
!
! Subroutine to calculate the longitudinal distance S for the elements
! in a lattice.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat -- lat_struct:
!
! Output:
!   lat -- lat_struct:
!-

subroutine s_calc (lat)

use bmad_interface, except_dummy => s_calc

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, slave0
type (branch_struct), pointer :: branch

integer i, j, n, ic, icon, ix2
real(8) ss, s_end

! Just go through all the elements and add up the lengths.

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (.not. bmad_com%auto_bookkeeper .and. branch%param%bookkeeping_state%s_position /= stale$) cycle

  ! Branches that branch from another branch start from zero
  ele => branch%ele(0)
  if (branch%ix_from_branch > -1) ele%s = 0  
  if (ele%bookkeeping_state%s_position == stale$) ele%bookkeeping_state%s_position = ok$

  ss = ele%s
  do n = 0, branch%n_ele_track
    ele => branch%ele(n)
    if (ele%bookkeeping_state%s_position == stale$) ele%bookkeeping_state%s_position = ok$
    ss = ss + ele%value(l$)
    ele%s = ss
  enddo

  branch%param%total_length = ss - branch%ele(0)%s
  branch%param%bookkeeping_state%s_position = ok$
enddo

! Now fill in the s positions of the super_lords and zero everyone else.
! Exception: A null_ele lord element is the result of a superposition on a multipass section.
! We need to preserve the s value of this element.

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n)
  if (.not. bmad_com%auto_bookkeeper .and. lord%bookkeeping_state%s_position /= stale$) cycle
  if (lord%key == null_ele$) cycle
  if (lord%n_slave == 0) cycle  ! Can happen when manipulating a lattice.

  select case (lord%lord_status)
  case (super_lord$, overlay_lord$)
    slave => pointer_to_slave(lord, lord%n_slave)
    lord%s = slave%s - lord%value(lord_pad2$)
  case (girder_lord$)
    slave0 => pointer_to_slave(lord, 1)
    do
      if (slave0%n_slave == 0) exit
      slave0 => pointer_to_slave(slave0, 1)
    enddo
    slave0 => pointer_to_next_ele (slave0, -1)

    slave => pointer_to_slave(lord, lord%n_slave)
    do
      if (slave%n_slave == 0) exit
      slave => pointer_to_slave(slave, 1)
    enddo

    lord%s = slave%s
    lord%value(ds_path_length$) = slave%s - slave0%s
    if (lord%value(ds_path_length$) < 0) lord%value(ds_path_length$) = &
                              lord%value(ds_path_length$) + slave0%branch%param%total_length
  case default
    lord%s = 0
  end select

  lord%bookkeeping_state%s_position = ok$
enddo

lat%lord_state%s_position = ok$

end subroutine
