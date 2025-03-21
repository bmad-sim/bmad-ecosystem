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

use equal_mod, except_dummy => multipass_all_info

implicit none

type (lat_struct), target :: lat
type (multipass_all_info_struct), target :: info
type (ele_struct), pointer :: m_lord, super_lord, slave1, slave
type (multipass_ele_info_struct), pointer :: s_info
type (branch_struct), pointer :: branch

integer i, j, k, n, ib, ik, ie, nl
integer n_multi_lord, n_pass, ix_pass, n_super_slave

! First get the number of multipass_lords.

if (allocated(info%lord)) deallocate(info%lord)
if (allocated(info%branch)) deallocate(info%branch)

n_multi_lord = 0
do ie = lat%n_ele_track+1, lat%n_ele_max
  m_lord => lat%ele(ie)
  if (m_lord%lord_status /= multipass_lord$) cycle
  n_multi_lord = n_multi_lord + 1
enddo

allocate (info%lord(n_multi_lord), info%branch(0:ubound(lat%branch, 1)))
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  allocate (info%branch(ib)%ele(0:branch%n_ele_max))
  info%branch(ib)%ele(:)%multipass = .false.
  info%branch(ib)%ele(:)%ix_pass = -1
  do i = 0, branch%n_ele_max
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

