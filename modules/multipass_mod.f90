module multipass_mod

use bmad_struct
use bmad_interface

! Multipass_top_info_struct gives complete information about a single 
! multipass_lord and all its slaves ("from the top down"). 
! If the multipass_lord has super_lords as slaves, n_super will be the number of 
! super_slaves per super_lord.
! ix_slave(1:n_pass, 1:n_super_slave) is a matrix of slaves in the tracking lattice.
! If there are no super_lords then n_super_slave = 1.and 
!      ix_super_lord(1:n_pass) = ix_slave (1:n_pass, 1)

type multipass_top_info_struct
  integer ix_lord          ! Lord index in lat%ele(:)
  integer n_pass           ! Number of passes (= number of slaves)
  integer n_super_slave    ! Number of super_slaves per super_lord. 
  integer, allocatable :: ix_super_lord(:)  ! Super_lord list in lat%ele(:) if they exist.
  integer, allocatable :: ix_slave(:,:)     ! Slaves list in lat%ele(:) tracking part.
end type

! Multipass_bottom_info_struct gives information about a singe
! multipass_slave ("from the bottom up").

type multipass_bottom_info_struct
  logical multipass     ! True if involved in multipass. False otherwise
  integer ix_top        ! Pointer to top(:) array
  integer ix_pass       ! Pass number
  integer ix_super      ! Index to ix_slave(ix_pass, ix_super_slave) matrix
end type

! top(i), i = 1, ..., n = number of multipass_lords in the lattice.
! bottom(i), i = 1, ..., n = lat%n_ele_track.

type multipass_all_info_struct
  type (multipass_top_info_struct), allocatable :: top(:)      ! Array of lords
  type (multipass_bottom_info_struct), allocatable :: bottom(:)
end type

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine multipass_all_info (lat, info)
!
! Subroutine to put multipass to a multipass_all_info_struct structure.
!
! Modules needed:
!   use multipass_mod
!
! Input:
!   lat -- lat_struct: Lattice
!
! Output
!   info -- Multipass_all_info_struct: Multipass information.
!-

subroutine multipass_all_info (lat, info)

implicit none

type (lat_struct), target :: lat
type (multipass_all_info_struct) info
type (ele_struct), pointer :: m_lord, super_lord

integer j, k, ik, ie, nl, ixsl, ixss 
integer n_multi_lord, n_pass, ix_pass, ix_slave, n_super_slave

! First get the number of multipass_lords.

n_multi_lord = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie)
  if (m_lord%control_type /= multipass_lord$) cycle
  n_multi_lord = n_multi_lord + 1
enddo

if (allocated (info%top)) deallocate (info%top, info%bottom)
allocate (info%top(n_multi_lord), info%bottom(lat%n_ele_max))
info%bottom(:)%multipass = .false.
info%bottom(:)%ix_top = -1
info%bottom(:)%ix_pass = -1
info%bottom(:)%ix_super = -1

! Fill in rest of the information

nl = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie)
  if (m_lord%control_type /= multipass_lord$) cycle
  nl = nl + 1

  info%top(nl)%ix_lord = ie
  n_pass = m_lord%n_slave
  info%top(nl)%n_pass = n_pass
  info%bottom(ie)%multipass = .true.
  info%bottom(ie)%ix_top = nl

  allocate (info%top(nl)%ix_super_lord(n_pass))

  ix_slave = lat%control(m_lord%ix1_slave)%ix_slave

  if (lat%ele(ix_slave)%control_type == super_lord$) then
    n_super_slave = lat%ele(ix_slave)%n_slave
    info%top(nl)%n_super_slave = n_super_slave
    allocate (info%top(nl)%ix_slave(n_pass, n_super_slave))
    do j = m_lord%ix1_slave, m_lord%ix2_slave
      ix_pass = j + 1 - m_lord%ix1_slave
      ixsl = lat%control(j)%ix_slave
      super_lord => lat%ele(ixsl)
      info%top(nl)%ix_super_lord(ix_pass) = ixsl
      info%bottom(ixsl)%multipass = .true.
      info%bottom(ixsl)%ix_top = nl
      info%bottom(ixsl)%ix_pass = ix_pass
      do k = super_lord%ix1_slave, super_lord%ix2_slave
        ik = k + 1 - super_lord%ix1_slave
        ixss = lat%control(k)%ix_slave
        info%top(nl)%ix_slave(ix_pass, ik) = ixss
        info%bottom(ixss)%multipass = .true.
        info%bottom(ixss)%ix_top = nl
        info%bottom(ixss)%ix_pass = ix_pass
        info%bottom(ixss)%ix_super = ik
      enddo
    enddo

  else
    info%top(nl)%n_super_slave = 1
    allocate (info%top(nl)%ix_slave(n_pass, 1))
    do j = m_lord%ix1_slave, m_lord%ix2_slave
      ix_pass = j + 1 - m_lord%ix1_slave
      ixss = lat%control(j)%ix_slave
      info%top(nl)%ix_super_lord(ix_pass) = ixss
      info%top(nl)%ix_slave(ix_pass, 1) = ixss
      info%bottom(ixss)%multipass = .true.
      info%bottom(ixss)%ix_top = nl
      info%bottom(ixss)%ix_pass = ix_pass
      info%bottom(ixss)%ix_super = 1
    enddo
  endif

enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function multipass_lord_index (ix_ele, lat, ix_pass, ix_super_lord) result (ix_multi_lord)
!
! Routine to find the index of the multipass lord of a lattice element.
! ix_multi_lord will be positive for elements:
!   multipass_slaves
!   super_lords that are slaves of a multipass_lord
!   super_slaves whose super_lord is a slave of a multipass_lord
!
! Modules needed:
!   use bmad
!
! Input:
!   ix_ele   -- Integer: Index of the lattice element.
!   lat      -- Lat_struct: Lattice containing the element.
!
! Output:
!   ix_pass       -- Integer, optional: Multipass turn number.
!   ix_super_lord -- Integer, optional: Index of the super_lord of the element.
!                      Set to -1 if the element is not a super_slave.
!   ix_multi_lord -- Integer: Index of the multipass_lord if there is one.
!                      Set to -1 if there is no multipass_lord.
!-

function multipass_lord_index (ix_ele, lat, ix_pass, ix_super_lord) result (ix_multi_lord)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer ix_ele, ix_multi_lord
integer, optional :: ix_super_lord, ix_pass
integer ix_sup, ix_mult, ic

!

ele => lat%ele(ix_ele)

ix_multi_lord = -1
if (present(ix_super_lord)) ix_super_lord = -1

if (ele%control_type == multipass_slave$) then
  ic = lat%ic(ele%ic1_lord)
  ix_multi_lord = lat%control(ic)%ix_lord
  if (present(ix_pass)) ix_pass = ic + 1 - lat%ele(ix_multi_lord)%ix1_slave
  return
endif

if (ele%control_type == super_slave$) then
  ic = lat%ic(ele%ic1_lord)
  ix_sup = lat%control(ic)%ix_lord
  if (lat%ele(ix_sup)%n_lord == 0) return
  ic = lat%ic(lat%ele(ix_sup)%ic1_lord)
  ix_mult = lat%control(ic)%ix_lord
  if (lat%ele(ix_mult)%control_type /= multipass_lord$) return
  ix_multi_lord = ix_mult
  if (present(ix_super_lord)) ix_super_lord = ix_sup
  if (present(ix_pass)) ix_pass = ic + 1 - lat%ele(ix_multi_lord)%ix1_slave
  return
endif

if (ele%control_type == super_lord$) then
  if (ele%n_lord == 0) return
  ic = lat%ic(ele%ic1_lord)
  ix_mult = lat%control(ic)%ix_lord
  if (lat%ele(ix_mult)%control_type /= multipass_lord$) return
  ix_multi_lord = ix_mult
  if (present(ix_pass)) ix_pass = ic + 1 - lat%ele(ix_multi_lord)%ix1_slave
endif

end function

end module
