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

integer i, ix, ios

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

if (index(ele_id, ':') == 0 .and. index(ele_id, '##') == 0 .and. b_str == '') then
  if (allocated(shape%uni)) deallocate(shape%uni)
  return
endif

allocate (shape%uni(ubound(s%u,1)))
do i = 1, ubound(s%u, 1)
  call lat_ele_locator(shape%ele_id, s%u(i)%model%lat, shape%uni(i)%eles, shape%uni(i)%n_loc)
enddo

end subroutine tao_shape_init
