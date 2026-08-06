!+
! Subroutine s_calc (lat)
!
! Subroutine to calculate the longitudinal distance S for the elements
! in a lattice.
!
! Input:
!   lat -- lat_struct:
!
! Output:
!   lat -- lat_struct:
!-

subroutine s_calc (lat)

use bookkeeper_mod, except_dummy => s_calc

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, slave0
type (branch_struct), pointer :: branch, b2

integer i, j, n, ic, icon, ix2
real(rp) ss, tot_len, beta, time
logical s_shift, len_change

! Just go through all the elements and add up the lengths.

s_shift = .false.
len_change = .false.

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (branch%param%bookkeeping_state%s_position /= stale$) cycle
  s_shift = .true.

  ! Branches that branch from another branch start from s = 0
  ele => branch%ele(0)
  if (branch%ix_from_branch > -1) then
    ele%s_start = 0
    ele%s = 0  
  endif
  if (ele%bookkeeping_state%s_position == stale$) ele%bookkeeping_state%s_position = ok$

  ss = ele%s
  do n = 0, branch%n_ele_track
    ele => branch%ele(n)
    ele%s_start = ss
    ! Patch element have a length that is a dependent variable
    if (ele%key == patch$ .and. ele%bookkeeping_state%attributes /= ok$) call attribute_bookkeeper(ele)
    ss = ss + ele%value(l$)
    ele%s = ss
    if (ele%bookkeeping_state%s_position == stale$) ele%bookkeeping_state%s_position = ok$
  enddo

  tot_len = ss - branch%ele(0)%s
  len_change = (abs(branch%param%total_length - tot_len) > 1d-15 * tot_len)
  branch%param%total_length = tot_len
  branch%param%bookkeeping_state%s_position = ok$
enddo

! Now fill in the s positions of the lords.
! Note: The s-position of a overlay, group, or control lord will not make sense if the lord
! controls multiple disjoint elements.


if (.not. s_shift .and. lat%lord_state%s_position /= stale$) return
lat%lord_state%s_position = ok$

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n)
  lord%bookkeeping_state%s_position = ok$

  ! Important: Do not mangle null_eles since null_eles in the lord section are used by bmad_parser 
  ! to preserve information on placement of drifts or null_eles that have been superimposed upon.
  if (lord%key == null_ele$) cycle
  if (lord%n_slave == 0) cycle  ! Can happen when manipulating a lattice during parsing.

  select case (lord%lord_status)
  case (super_lord$, overlay_lord$, group_lord$)
    slave => pointer_to_slave(lord, 1)
    lord%s_start = slave%s_start + lord%value(lord_pad1$)
    slave => pointer_to_slave(lord, lord%n_slave)
    lord%s = slave%s - lord%value(lord_pad2$)
  case (girder_lord$)
    call find_element_ends (lord, slave0, slave)
    lord%s_start = slave0%s
    lord%s = slave%s
    lord%value(l$) = slave%s - slave0%s
    if (lord%value(l$) < 0) lord%value(l$) = lord%value(l$) + slave0%branch%param%total_length
  case default ! multipass_lord and control_lord elements do not have an s-position.
    lord%s_start = 0
    lord%s = 0
  end select
enddo

!

if (len_change) then
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do n = 0, branch%n_ele_max
      ele => branch%ele(n)
      if (ele%value(e_tot$) == 0) cycle
                          
      select case (ele%key)
      case (crab_cavity$, rfcavity$, rf_bend$)
        beta = ele%value(p0c$) / ele%value(e_tot$)
        b2 => pointer_to_branch(ele)
        time = branch%param%total_length / (c_light * beta)
        if (time == 0) cycle
        if (is_true(ele%value(harmon_master$))) then
          ele%value(rf_frequency$) = ele%value(harmon$) / time
        else
          ele%value(harmon$) = ele%value(rf_frequency$) * time
        endif
      end select

    enddo
  enddo
endif

end subroutine
