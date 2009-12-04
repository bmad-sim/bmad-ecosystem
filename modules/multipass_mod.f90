module multipass_mod

use bmad_struct
use bmad_interface

! Multipass_top_info_struct gives complete information about a single 
! multipass_lord and all its slaves ("from the top down"). 
! If the multipass_lord has super_lords as slaves, n_super will be the number of 
! super_slaves per super_lord.
! slave(1:n_pass, 1:n_super_slave) is a matrix of slaves in the tracking lattice.
! If there are no super_lords then n_super_slave = 1 and 
!      super_lord(1:n_pass) = slave (1:n_pass, 1)

type multipass_top_info_struct
  type (ele_struct), pointer :: lord          ! Lord element
  integer n_pass           ! Number of passes (= number of slaves)
  integer n_super_slave    ! Number of super_slaves per super_lord. 
  type (ele_pointer_struct), allocatable :: super_lord(:)  ! Super_lord list if they exist.
  type (ele_pointer_struct), allocatable :: slave(:,:)     ! Slaves list in tracking part.
end type

! Multipass_bottom_info_struct gives information about a singe
! multipass_slave ("from the bottom up").

type multipass_bottom_info_struct
  logical multipass     ! True if involved in multipass. False otherwise
  integer ix_pass       ! Pass number
  integer, allocatable :: ix_top(:)   ! Pointers to top(:) array
  integer, allocatable :: ix_super(:) ! Indexes to slave(ix_pass, super_slave%ix_ele) matrix
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
! Subroutine to put multipass information into a multipass_all_info_struct structure.
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
type (ele_struct), pointer :: m_lord, super_lord, slave1, slave

integer i, j, k, n, ik, ie, nl
integer n_multi_lord, n_pass, ix_pass, n_super_slave

! First get the number of multipass_lords.

n_multi_lord = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie)
  if (m_lord%lord_status /= multipass_lord$) cycle
  n_multi_lord = n_multi_lord + 1
enddo

if (allocated (info%top)) deallocate (info%top, info%bottom)
allocate (info%top(n_multi_lord), info%bottom(lat%n_ele_max))
info%bottom(:)%multipass = .false.
info%bottom(:)%ix_pass = -1
do i = 1, lat%n_ele_max
  allocate (info%bottom(i)%ix_top(0))
  allocate (info%bottom(i)%ix_super(0))
enddo

! Fill in rest of the information

nl = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie) 
  if (m_lord%lord_status /= multipass_lord$) cycle
  nl = nl + 1

  info%top(nl)%lord => m_lord
  n_pass = m_lord%n_slave
  info%top(nl)%n_pass = n_pass
  info%bottom(ie)%multipass = .true.
  n = size(info%bottom(ie)%ix_top)
  call re_allocate(info%bottom(ie)%ix_top, n+1)
  info%bottom(ie)%ix_top(n+1) = nl

  allocate (info%top(nl)%super_lord(n_pass))

  slave1 => pointer_to_slave(lat, m_lord, 1)

  if (slave1%lord_status == super_lord$) then
    n_super_slave = slave1%n_slave
    info%top(nl)%n_super_slave = n_super_slave
    allocate (info%top(nl)%slave(n_pass, n_super_slave))
    do j = 1, m_lord%n_slave
      super_lord => pointer_to_slave(lat, m_lord, j)
      info%top(nl)%super_lord(j)%ele => super_lord
      info%bottom(super_lord%ix_ele)%multipass = .true.
      n = size(info%bottom(super_lord%ix_ele)%ix_top)
      call re_allocate(info%bottom(super_lord%ix_ele)%ix_top, n+1)
      info%bottom(super_lord%ix_ele)%ix_top(n+1) = nl
      info%bottom(super_lord%ix_ele)%ix_pass = j
      do k = 1, super_lord%n_slave
        slave => pointer_to_slave(lat, super_lord, k)
        info%top(nl)%slave(j, k)%ele => slave
        info%bottom(slave%ix_ele)%multipass = .true.
        info%bottom(slave%ix_ele)%ix_pass = j
        n = size(info%bottom(slave%ix_ele)%ix_top)
        call re_allocate(info%bottom(slave%ix_ele)%ix_top, n+1)
        call re_allocate(info%bottom(slave%ix_ele)%ix_super, n+1)
        info%bottom(slave%ix_ele)%ix_top(n+1) = nl
        info%bottom(slave%ix_ele)%ix_super(n+1) = k
      enddo
    enddo

  else
    info%top(nl)%n_super_slave = 1
    allocate (info%top(nl)%slave(n_pass, 1))
    do j = 1, m_lord%n_slave
      slave => pointer_to_slave(lat, m_lord, j)
      info%top(nl)%super_lord(j)%ele => slave
      info%top(nl)%slave(j, 1)%ele => slave
      info%bottom(slave%ix_ele)%multipass = .true.
      info%bottom(slave%ix_ele)%ix_pass = j
      n = size(info%bottom(slave%ix_ele)%ix_top)
      call re_allocate(info%bottom(slave%ix_ele)%ix_top, n+1)
      call re_allocate(info%bottom(slave%ix_ele)%ix_super, n+1)
      info%bottom(slave%ix_ele)%ix_top(n+1) = nl
      info%bottom(slave%ix_ele)%ix_super(n+1) = 1
    enddo
  endif

enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function pointer_to_multipass_lord (ele, lat, ix_pass, super_lord) result (multi_lord)
!
! Routine to find the multipass lord of a lattice element.
! A multi_lord will be found for:
!   multipass_slaves
!   super_lords that are slaves of a multipass_lord
!   super_slaves whose super_lord is a slave of a multipass_lord
!
! Modules needed:
!   use multipass_mod
!
! Input:
!   ele   -- Ele_struct: Lattice element.
!   lat   -- Lat_struct: Lattice containing the element.
!
! Output:
!   ix_pass    -- Integer, optional: Multipass turn number.
!                      Set to -1 if element is not a multipass slave
!   super_lord -- Ele_struct, pointer, optional: super_lord of the element.
!                      Set to NULL if ele is not a super_slave.
!   multi_lord -- Ele_struct, pointer: multipass_lord if there is one.
!                      Set to NULL if there is no multipass_lord.
!-

function pointer_to_multipass_lord (ele, lat, ix_pass, super_lord) result (multi_lord)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ele
type (ele_struct), pointer :: multi_lord, sup_lord
type (ele_struct), pointer, optional :: super_lord

integer, optional :: ix_pass
integer ix_con

!

if (present(super_lord)) nullify (super_lord)
if (present(ix_pass)) ix_pass = -1

if (ele%slave_status == multipass_slave$) then
  multi_lord => pointer_to_lord (lat, ele, 1, ix_con)
  if (present(ix_pass)) ix_pass = ix_con + 1 - multi_lord%ix1_slave
  return
endif

if (ele%slave_status == super_slave$) then
  sup_lord => pointer_to_lord (lat, ele, 1)
  if (present(super_lord)) super_lord = sup_lord

  if (sup_lord%n_lord /= multipass_slave$) return
  multi_lord => pointer_to_lord (lat, sup_lord, 1, ix_con)
  if (present(ix_pass)) ix_pass = ix_con + 1 - multi_lord%ix1_slave
  return
endif

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine multipass_chain (ele, lat, ix_pass, n_links, chain_ele)
!
! Routine to return the chain of elements that represent the same physical element
! when there is multipass.
!
! Modules needed:
!   use multipass_mod
!
! Input:
!   ele  -- Ele_pointer_struct: Element in a multipass chain.
!   lat  -- Lat_struct: Lattice structure.
!
! Output
!   ix_pass      -- Integer: Multipass pass number of the input element. 
!                     Set to -1 if input element is not in a multipass section.
!   n_links      -- Integer, optional: Number of times the physical element is passed through.
!   chain_ele(:) -- Ele_pointer_struct, optional, allocatable: pointers to the
!                    elements of the chain. Note: chain_ele(ix_pass)%ele => ele
!-

subroutine multipass_chain (ele, lat, ix_pass, n_links, chain_ele)

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: m_lord, s_lord, slave
type (ele_pointer_struct), allocatable, optional :: chain_ele(:)

integer i, j, k, ix_pass, ic, ix_lord, ix_off, ix_c, n_links

! Init

ix_pass = -1
n_links = 0

! element is a multipass_slave case

if (ele%slave_status == multipass_slave$) then
  m_lord => pointer_to_lord(lat, ele, 1)
  if (present(chain_ele)) call re_allocate_eles (chain_ele, m_lord%n_slave, .false.)
  n_links = m_lord%n_slave
  do j = 1, m_lord%n_slave
    if (present(chain_ele)) chain_ele(j)%ele => pointer_to_slave(lat, m_lord, j)
    if (chain_ele(j)%ele%ix_ele  == ele%ix_ele) ix_pass = j
  enddo
endif

! element is a super_slave case

if (ele%slave_status == super_slave$) then
  s_lord => pointer_to_lord(lat, ele, 1)
  if (s_lord%slave_status /= multipass_slave$) return

  ! Find offset in super_lord

  do j = 1, s_lord%n_slave
    slave => pointer_to_slave(lat, ele, j)
    if (slave%ix_ele == ele%ix_ele) ix_off = j
  enddo

  ! Construct chain

  m_lord => pointer_to_lord (lat, ele, 1)
  if (present(chain_ele)) call re_allocate_eles (chain_ele, m_lord%n_slave, .false.)
  n_links = m_lord%n_slave
  do j = 1, m_lord%n_slave
    s_lord => pointer_to_slave(lat, m_lord, j)
    slave => pointer_to_slave(lat, s_lord, ix_off)
    if (present(chain_ele)) chain_ele(j)%ele => slave
    if (ix_c == ele%ix_ele) ix_pass = j
  enddo

endif

end subroutine multipass_chain

end module
