!+
! Subroutine ele_order_calc (lat, order)
!
! Routine to collect the information needed to construct unique lattice element name.
!
! The output of this routine is used in ele_unique_name.
!
! Input:
!   lat       -- lat_struct: Lattice to analyze.
!
! Output:
!   order     -- lat_ele_order_struct: Structure holding the element order information.
!-

subroutine ele_order_calc (lat, order)

use bmad_routine_interface, dummy => ele_order_calc

implicit none

type this_order1_struct
  integer ix_ele   ! element index in branch
  integer ix_nt    ! Nametable index
end type

type this_order_struct
  type (this_order1_struct), allocatable :: ele(:)
end type

type (this_order_struct), allocatable :: ord_br(:)
type (this_order1_struct), allocatable :: temp_ord(:)

type (lat_struct), target :: lat
type (lat_ele_order_struct) order
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele2
type (nametable_struct), pointer :: nt

integer nn, i0, i1, ib, ie, in, ix, ib_ub, ie_ub, it1, sum_count
integer, allocatable :: n_count(:), indx(:)

character(40) name0

!

ib_ub = ubound(lat%branch, 1)
allocate (n_count(0:ib_ub), ord_br(0:ib_ub))

if (allocated(order%branch)) then
  if (ubound(order%branch,1) /= ib_ub) deallocate(order%branch)
endif
if (.not. allocated(order%branch)) allocate (order%branch(0:ib_ub))

do ib = 0, ib_ub
  allocate(ord_br(ib)%ele(10))

  ie_ub = lat%branch(ib)%n_ele_max
  if (allocated(order%branch(ib)%ele)) then
    if (ubound(order%branch(ib)%ele,1) /= ie_ub) deallocate(order%branch(ib)%ele)
  endif
  if (.not. allocated(order%branch(ib)%ele)) allocate (order%branch(ib)%ele(0:ie_ub))
end do

!

nt => lat%nametable

i0 = 0

do

  n_count = 0
  name0 = nt%name(nt%index(i0))
  i1 = i0 - 1

  do
    i1 = i1 + 1
    if (i1 > nt%n_max) exit
    it1 = nt%index(i1)
    if (nt%name(it1) /= name0) exit
    ele => pointer_to_ele(lat, it1)
    ele2 => ele
    if (ele%lord_status == super_lord$) ele2 => pointer_to_slave(ele, 1)
    ib = ele2%ix_branch
    if (n_count(ib) + 1 > size(ord_br(ib)%ele)) then
      call move_alloc(ord_br(ib)%ele, temp_ord)
      allocate (ord_br(ib)%ele(2*n_count(ib)))
      ord_br(ib)%ele(1:n_count(ib)) = temp_ord
    endif
    n_count(ib) = n_count(ib) + 1
    nn = n_count(ib)
    ord_br(ib)%ele(nn)%ix_ele = ele2%ix_ele
    ord_br(ib)%ele(nn)%ix_nt = it1
  enddo

  sum_count = sum(n_count)
  do ib = 0, ib_ub
    nn = n_count(ib)
    if (nn == 0) cycle

    if (sum_count == 1 .or. nn == 1) then
      in = ord_br(ib)%ele(1)%ix_nt
      ele => pointer_to_ele(lat, in)
      order%branch(ele%ix_branch)%ele(ele%ix_ele)%ix_branch = ib
      if (sum_count == 1) then
        order%branch(ele%ix_branch)%ele(ele%ix_ele)%ix_order = -1  ! Unique in lattice
      else
        order%branch(ele%ix_branch)%ele(ele%ix_ele)%ix_order = 0  ! Unique in branch
      endif
      cycle
    endif

    call re_allocate(indx, nn, .false.)
    call indexer(ord_br(ib)%ele(1:nn)%ix_ele, indx(1:nn))
    do ix = 1, nn
      in = ord_br(ib)%ele(indx(ix))%ix_nt
      ele => pointer_to_ele(lat, in)
      order%branch(ele%ix_branch)%ele(ele%ix_ele)%ix_branch = ib
      order%branch(ele%ix_branch)%ele(ele%ix_ele)%ix_order = ix
    enddo
  enddo

  if (i1 > nt%n_max) exit
  i0 = i1
enddo

end subroutine ele_order_calc
