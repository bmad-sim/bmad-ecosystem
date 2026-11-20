!+
! Subroutine tao_shape_init (shape, err, print_err)
!
! Routine to split shape%ele_id in the form str = "xxx::yyy" into an element class
! and an element name. Example: 
!   shape%ele_id = "quad::q*".
! gives
!   ix_class = quadrupole$
!   ele_name = "Q*"
!
! If "::" is not found then ix_class is set to 0 (all classes).
! If str is of the form: "*::yyy" then ix_class is set to 0 (all classes).
! Class abbreviations and lower case names accepted 
! ele_name will be converted to upper case
!
! Input:
!   shape     -- tao_ele_shape_struct: Shape
!   print_err -- Logical, optional: If True then print an error message if 
!                 there is a problem. Default is True.
!
! Output:
!   shape     -- tao_ele_shape_struct: Shape
!   err       -- Set true if there is a problem translating the element class.
!-

subroutine tao_shape_init (shape, err, print_err)

use tao_interface, dummy => tao_shape_init

implicit none

type (tao_ele_shape_struct) shape
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele, lord
type (lat_struct), pointer :: lat

type in_list_struct
  logical, allocatable :: ele(:)
end type
type (in_list_struct), allocatable :: in_list(:)

integer i, ix, ios, ib, ie, il, n

character(*), parameter :: r_name = 'tao_shape_init'
character(20) class, b_str
character(60) ele_id

logical, optional :: print_err
logical err

!

err = .false.
shape%ix_key = -1
ele_id = shape%ele_id

ix = index(ele_id, '::')
if (ix /= 0) then
  if (ele_id(1:ix-1) == 'data') return
  if (ele_id(1:ix-1) == 'var') return
  if (ele_id(1:ix-1) == 'lat') return
  if (ele_id(1:ix-1) == 'alias') return
  if (ele_id(1:ix-1) == 'type') return
  if (ele_id(1:ix-1) == 'building_wall') return
endif

!

ix = index(ele_id, '>>')
b_str = ''
if (ix /= 0) then
  b_str = ele_id(1:ix-1)
  ele_id = ele_id(ix+2:)
endif

!

shape%ix_key = 0
class = ''
ix = index(ele_id, '::')
if (ix /= 0) then
  class = ele_id(:ix-1)
  ele_id = ele_id(ix+2:)
  if (class /= '*') then
    shape%ix_key = key_name_to_key_index (class, .true.)
    if (shape%ix_key < 1) then
      if (logic_option (.true., print_err)) call out_io (s_error$, r_name, 'BAD ELEMENT CLASS NAME: ' // class)
      err = .true.
    endif
  endif
endif

shape%name_ele = upcase(ele_id)

! A position dependent name to match to is something like "ele1:ele2" which is a range construct
! or "name##N" signifying the N^th element. A position dependent ele_id name means that simple element name
! matching will not work. In this case store the results of lat_ele_locator.
! The reason why lat_ele_locator is not used in all cases is that with something like ele_id = "Q*" will match
! to many elements so in large lattices this might be significant time consideration.

if (allocated(shape%uni)) deallocate(shape%uni)
if (index(ele_id, ':') == 0 .and. index(ele_id, '##') == 0 .and. b_str == '') return

allocate (shape%uni(ubound(s%u,1)))

! A range construct like "q1:q3" gives the tracking elements in the range so a super_lord in the range 
! is represented by its super_slaves. But the lat_layout and floor_plan drawing routines skip super_slaves 
! and draw the super_lord instead. So replace any super_slave in the list with its super_lord(s).
! Also %loc is set for all elements since tao_pointer_to_ele_shape matches by element index.

do i = 1, ubound(s%u, 1)
  lat => s%u(i)%model%lat
  call lat_ele_locator(shape%ele_id, lat, eles, n)

  if (allocated(in_list)) deallocate(in_list)
  allocate (in_list(0:ubound(lat%branch, 1)))
  do ib = 0, ubound(lat%branch, 1)
    allocate (in_list(ib)%ele(0:lat%branch(ib)%n_ele_max))
    in_list(ib)%ele = .false.
  enddo

  allocate (shape%uni(i)%eles(max(1, n)))
  shape%uni(i)%n_loc = 0

  do ie = 1, n
    ele => eles(ie)%ele
    if (ele%slave_status == super_slave$) then
      do il = 1, ele%n_lord
        lord => pointer_to_lord(ele, il)
        if (lord%lord_status /= super_lord$) cycle
        call add_to_list(shape%uni(i), lord)
      enddo
    else
      call add_to_list(shape%uni(i), ele)
    endif
  enddo
enddo

!------------------------------------------------------------------------
contains

subroutine add_to_list(su, ele)

type (tao_ele_pointer_struct) su
type (ele_struct), target :: ele
integer nn

!

if (in_list(ele%ix_branch)%ele(ele%ix_ele)) return  ! Already in list.
in_list(ele%ix_branch)%ele(ele%ix_ele) = .true.

nn = su%n_loc + 1
if (nn > size(su%eles)) call re_allocate_eles(su%eles, 2*nn, .true.)
su%eles(nn)%ele => ele
su%eles(nn)%loc = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)
su%n_loc = nn

end subroutine add_to_list

end subroutine tao_shape_init
