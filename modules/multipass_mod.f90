module multipass_mod

use bmad_struct
use bmad_interface

! Multipass_lord_info_struct gives complete information about a single 
! multipass_lord and all its slaves ("from the top down"). 
! If the multipass_lord has super_lords as slaves, n_super will be the number of 
! super_slaves per super_lord.
! slave(1:n_pass, 1:n_super_slave) is a matrix of slaves in the tracking lattice.
! If there are no super_lords then n_super_slave = 1 and 
!      super_lord(1:n_pass) = slave (1:n_pass, 1)

type multipass_lord_info_struct
  type (ele_struct), pointer :: lord          ! Lord element
  integer n_pass           ! Number of passes (= number of slaves)
  integer n_super_slave    ! Number of super_slaves per super_lord. 
  type (ele_pointer_struct), allocatable :: super_lord(:)  ! Super_lord list if they exist.
  type (ele_pointer_struct), allocatable :: slave(:,:)     ! Slaves list in tracking part.
end type

! Multipass_ele_info_struct gives information about a singe element in the lattice
! ("from the bottom up").

type multipass_ele_info_struct
  logical multipass     ! True if involved in multipass. False otherwise
  integer ix_pass       ! Pass number
  integer, allocatable :: ix_lord(:)   ! Pointers to lord(:) array
  integer, allocatable :: ix_super(:) ! Indexes to slave(ix_pass, super_slave%ix_ele) matrix
end type

type multipass_branch_info_struct
  type (multipass_ele_info_struct), allocatable :: ele(:)
end type

! %lord(i), i = 1, ..., n = number of multipass_lords in the lattice.
! %branch(i), i = 0, ..., n = ubound(lat%branch)

type multipass_all_info_struct
  type (multipass_lord_info_struct), allocatable :: lord(:)      ! Array of lords
  type (multipass_branch_info_struct), allocatable :: branch(:)
end type

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine multipass_all_info (lat, info)
!
! Routine to put multipass information into a multipass_all_info_struct structure.
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
type (multipass_all_info_struct), target :: info
type (ele_struct), pointer :: m_lord, super_lord, slave1, slave
type (multipass_ele_info_struct), pointer :: s_info
type (branch_struct), pointer :: branch

integer i, j, k, n, ib, ik, ie, nl
integer n_multi_lord, n_pass, ix_pass, n_super_slave

! First get the number of multipass_lords.

call deallocate_multipass_all_info_struct (info)

n_multi_lord = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie)
  if (m_lord%lord_status /= multipass_lord$) cycle
  n_multi_lord = n_multi_lord + 1
enddo

allocate (info%lord(n_multi_lord), info%branch(0:ubound(lat%branch, 1)))
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  allocate (info%branch(ib)%ele(branch%n_ele_max))
  info%branch(ib)%ele(:)%multipass = .false.
  info%branch(ib)%ele(:)%ix_pass = -1
  do i = 1, branch%n_ele_max
    s_info => info%branch(ib)%ele(i)
    allocate (s_info%ix_lord(0))
    allocate (s_info%ix_super(0))
  enddo
enddo

! Fill in rest of the information

nl = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie) 
  if (m_lord%lord_status /= multipass_lord$) cycle
  nl = nl + 1

  info%lord(nl)%lord => m_lord
  n_pass = m_lord%n_slave
  info%lord(nl)%n_pass = n_pass
  s_info => info%branch(0)%ele(ie)
  s_info%multipass = .true.
  n = size(s_info%ix_lord)
  call re_allocate(s_info%ix_lord, n+1)
  s_info%ix_lord(n+1) = nl

  allocate (info%lord(nl)%super_lord(n_pass))

  slave1 => pointer_to_slave(m_lord, 1)

  if (slave1%lord_status == super_lord$) then
    n_super_slave = slave1%n_slave
    info%lord(nl)%n_super_slave = n_super_slave
    allocate (info%lord(nl)%slave(n_pass, n_super_slave))
    do j = 1, m_lord%n_slave
      super_lord => pointer_to_slave(m_lord, j)
      info%lord(nl)%super_lord(j)%ele => super_lord
      s_info => info%branch(0)%ele(super_lord%ix_ele)
      s_info%multipass = .true.
      n = size(s_info%ix_lord)
      call re_allocate(s_info%ix_lord, n+1)
      s_info%ix_lord(n+1) = nl
      s_info%ix_pass = j
      do k = 1, super_lord%n_slave
        slave => pointer_to_slave(super_lord, k)
        info%lord(nl)%slave(j, k)%ele => slave
        s_info => info%branch(slave%ix_branch)%ele(slave%ix_ele)
        s_info%multipass = .true.
        s_info%ix_pass = j
        n = size(s_info%ix_lord)
        call re_allocate(s_info%ix_lord, n+1)
        call re_allocate(s_info%ix_super, n+1)
        s_info%ix_lord(n+1) = nl
        s_info%ix_super(n+1) = k
      enddo
    enddo

  else
    info%lord(nl)%n_super_slave = 1
    allocate (info%lord(nl)%slave(n_pass, 1))
    do j = 1, m_lord%n_slave
      slave => pointer_to_slave(m_lord, j)
      info%lord(nl)%super_lord(j)%ele => slave
      info%lord(nl)%slave(j, 1)%ele => slave
      s_info => info%branch(slave%ix_branch)%ele(slave%ix_ele)
      s_info%multipass = .true.
      s_info%ix_pass = j
      n = size(s_info%ix_lord)
      call re_allocate(s_info%ix_lord, n+1)
      call re_allocate(s_info%ix_super, n+1)
      s_info%ix_lord(n+1) = nl
      s_info%ix_super(n+1) = 1
    enddo
  endif

enddo

end subroutine multipass_all_info

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine deallocate_multipass_all_info_struct (info)
!
! Routine to deallocate the allocatable arrays in an multipass_all_info_struct.
!
! Modules needed:
!   use multipass_mod
!
! Input:
!   info -- Multipass_all_info_struct: Variable to deallocate
!
! Output
!   info -- Multipass_all_info_struct: Multipass information.
!-

subroutine deallocate_multipass_all_info_struct (info)

implicit none

type (multipass_all_info_struct) info

!

if (allocated(info%lord)) deallocate(info%lord)
if (allocated(info%branch)) deallocate(info%branch)

end subroutine deallocate_multipass_all_info_struct 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function pointer_to_multipass_lord (ele, ix_pass, super_lord) result (multi_lord)
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
!
! Output:
!   ix_pass    -- Integer, optional: Multipass turn number.
!                      Set to -1 if element is not a multipass slave
!   super_lord -- Ele_struct, pointer, optional: super_lord of the element.
!                      Set to NULL if ele is not a super_slave.
!   multi_lord -- Ele_struct, pointer: multipass_lord if there is one.
!                      Set to NULL if there is no multipass_lord.
!-

function pointer_to_multipass_lord (ele, ix_pass, super_lord) result (multi_lord)

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: multi_lord, sup_lord
type (ele_struct), pointer, optional :: super_lord

integer, optional :: ix_pass

!

nullify (multi_lord)
if (present(super_lord)) nullify (super_lord)
if (present(ix_pass)) ix_pass = -1

if (ele%slave_status == multipass_slave$) then
  multi_lord => pointer_to_lord(ele, 1, ix_slave = ix_pass)
  return
endif

if (ele%slave_status == super_slave$) then
  sup_lord => pointer_to_lord(ele, 1)
  if (present(super_lord)) super_lord => sup_lord

  if (sup_lord%slave_status /= multipass_slave$) return
  multi_lord => pointer_to_lord(sup_lord, 1, ix_slave = ix_pass)
endif

end function pointer_to_multipass_lord

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine multipass_chain (ele, ix_pass, n_links, chain_ele)
!
! Routine to return the chain of elements that represent the same physical element
! when there is multipass.
!
! Modules needed:
!   use multipass_mod
!
! Input:
!   ele  -- Ele_pointer_struct: Element in a multipass chain.
!
! Output
!   ix_pass              -- Integer: Multipass pass number of the input element. 
!                             Set to -1 if input element is not in a multipass section.
!   n_links              -- Integer: Number of times the physical element is passed through.
!   chain_ele(1:n_links) -- Ele_pointer_struct, optional, allocatable: pointers to the
!                             elements of the chain. Note: chain_ele(ix_pass)%ele => ele
!-

subroutine multipass_chain (ele, ix_pass, n_links, chain_ele)

implicit none

type (ele_struct) :: ele
type (ele_struct), pointer :: m_lord, s_lord, slave
type (ele_pointer_struct), allocatable, optional :: chain_ele(:)

integer :: n_links
integer i, j, ix_pass, ix_off

! Init

ix_pass = -1
n_links = 0

! element is a multipass_slave case

if (ele%slave_status == multipass_slave$) then
  m_lord => pointer_to_lord(ele, 1)
  if (present(chain_ele)) call re_allocate_eles (chain_ele, m_lord%n_slave)
  n_links = m_lord%n_slave
  do j = 1, m_lord%n_slave
    slave => pointer_to_slave(m_lord, j)
    if (present(chain_ele)) chain_ele(j)%ele => slave
    if (slave%ix_ele == ele%ix_ele .and. slave%ix_branch == ele%ix_branch) ix_pass = j
  enddo
endif

! element is a super_slave case

if (ele%slave_status == super_slave$) then
  s_lord => pointer_to_lord(ele, 1)
  if (s_lord%slave_status /= multipass_slave$) return

  ! Find offset in super_lord

  do j = 1, s_lord%n_slave
    slave => pointer_to_slave(s_lord, j)
    if (slave%ix_ele /= ele%ix_ele .or. slave%ix_branch /= ele%ix_branch) cycle
    ix_off = j
    exit
  enddo

  ! Construct chain

  m_lord => pointer_to_lord(s_lord, 1)
  if (present(chain_ele)) call re_allocate_eles (chain_ele, m_lord%n_slave)
  n_links = m_lord%n_slave
  do j = 1, m_lord%n_slave
    s_lord => pointer_to_slave(m_lord, j)
    slave => pointer_to_slave(s_lord, ix_off)
    if (present(chain_ele)) chain_ele(j)%ele => slave
    if (slave%ix_ele == ele%ix_ele .and. slave%ix_branch == ele%ix_branch) ix_pass = j
  enddo

endif

end subroutine multipass_chain

end module
