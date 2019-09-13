!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine multipass_chain (ele, ix_pass, n_links, chain_ele, use_super_lord)
!
! Routine to return the chain of elements that represent the same physical element
! when there is multipass.
!
! Input:
!   ele                  -- Ele_pointer_struct: Element in a multipass chain.
!   use_super_lord       -- logical, optional: If present and True and if ele is a super_slave, 
!                             construct the chain_ele(:) array using the corresponding super_lords.
!
! Output
!   ix_pass              -- Integer: Multipass pass number of the input element. 
!                             Set to -1 if input element is not in a multipass section.
!   n_links              -- Integer: Number of times the physical element is passed through.
!   chain_ele(1:n_links) -- Ele_pointer_struct, optional, allocatable: pointers to the
!                             elements of the chain. Note: chain_ele(ix_pass)%ele => ele
!-

subroutine multipass_chain (ele, ix_pass, n_links, chain_ele, use_super_lord)

use bmad_routine_interface, except_dummy => multipass_chain

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: m_lord, s_lord, slave, ele2
type (ele_pointer_struct), allocatable, optional :: chain_ele(:)

integer :: n_links
integer i, j, ix_pass, ix_off
logical, optional :: use_super_lord

! Init

ix_pass = -1
n_links = 0

if (logic_option(.false., use_super_lord) .and. ele%slave_status == super_slave$) then
  ele2 => pointer_to_lord(ele, 1)
else
  ele2 => ele
endif

! element is a multipass_slave case

if (ele2%slave_status == multipass_slave$) then
  m_lord => pointer_to_lord(ele2, 1)
  if (present(chain_ele)) call re_allocate_eles (chain_ele, m_lord%n_slave)
  n_links = m_lord%n_slave
  do j = 1, m_lord%n_slave
    slave => pointer_to_slave(m_lord, j)
    if (present(chain_ele)) chain_ele(j)%ele => slave
    if (slave%ix_ele == ele2%ix_ele .and. slave%ix_branch == ele2%ix_branch) ix_pass = j
  enddo
endif

! element is a super_slave case

if (ele2%slave_status == super_slave$) then
  s_lord => pointer_to_lord(ele2, 1)
  if (s_lord%slave_status /= multipass_slave$) return

  ! Find offset in super_lord

  do j = 1, s_lord%n_slave
    slave => pointer_to_slave(s_lord, j)
    if (slave%ix_ele /= ele2%ix_ele .or. slave%ix_branch /= ele2%ix_branch) cycle
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
    if (slave%ix_ele == ele2%ix_ele .and. slave%ix_branch == ele2%ix_branch) ix_pass = j
  enddo
endif

end subroutine multipass_chain

