!+
! Subroutine tao_ele_shape_info (ix_uni, ele, ele_shapes, e_shape, label_name, y1, y2, ix_shape_min)
!
! Routine to return info on an ele shape for a particular element.
!
! Input:
!   ix_uni        -- integer: Universe index.
!   ele           -- ele_struct: Lattice element.
!   ele_shapes(:) -- tao_ele_shape_struct: Array of shapes to search.
!   ix_shape_min  -- integer, optional: Index of minimum ele_shape(:) index to start search from. Default is 1.
!
! Output:
!   e_shape       -- tao_ele_shape_struct, pointer :: e_shape
!   label_name    -- character(*): Label name.
!   y1, y2        -- real(rp): shape transverse sizes.
!   ix_shape_min  -- integer, optional: Ele_shape(:) index to start next search if multiple shapes are associated with ele.
!-

subroutine tao_ele_shape_info (ix_uni, ele, ele_shapes, e_shape, label_name, y1, y2, ix_shape_min)

use tao_interface, dummy => tao_ele_shape_info

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: lord
type (tao_ele_shape_struct) ele_shapes(:)
type (tao_ele_shape_struct), pointer :: e_shape

real(rp) y1, y2, y, dat_var_value
integer, optional :: ix_shape_min
integer ix_uni, ix_shape, ix

character(*) label_name
character(40) dat_var_name, shape, prefix
character(*), parameter :: r_name = 'tao_ele_shape_info'

!

ix_shape = integer_option(1, ix_shape_min)
e_shape => tao_pointer_to_ele_shape (ix_uni, ele, ele_shapes, dat_var_name, dat_var_value, ix_shape)
if (present(ix_shape_min)) ix_shape_min = ix_shape

if (.not. associated(e_shape)) return

ix = index(e_shape%shape, ':')
if (ix == 0) then
  prefix = ''
  shape = e_shape%shape
else
  prefix = e_shape%shape(1:ix-1)
  shape = e_shape%shape(ix+1:)
endif

! offsets

y = e_shape%size

select case (prefix)
case ('var', 'asym_var')
  select case (ele%key)
  case (sbend$)
    y2 = y * ele%value(g$)
  case (quadrupole$)
    y2 = y * ele%value(k1$)
  case (sextupole$)
    y2 = y * ele%value(k2$)
  case (octupole$)
    y2 = y * ele%value(k3$)
  case (solenoid$)
    y2 = y * ele%value(ks$)
  end select
  y1 = y2
  if (prefix == 'asym_var') y2 = 0
case ('vvar')
  y1 = dat_var_value
  y2 = dat_var_value
case ('asym_vvar')
  y1 = dat_var_value
  y2 = 0
case ('', 'pattern')
  y1 = y
  y2 = y
case default
  call out_io (s_error$, r_name, 'Unknown shape prefix: ' // e_shape%shape)
  y1 = y
  y2 = y
endselect

! label_name

if (e_shape%label == 'name') then
  label_name = dat_var_name
  if (label_name == '') label_name = ele%name
  if (ele%slave_status == multipass_slave$ .and. label_name == ele%name) then
    lord => pointer_to_lord(ele, 1)
    label_name = lord%name
  endif

elseif (e_shape%label == 's') then
  write (label_name, '(f16.2)') ele%s - ele%value(l$) / 2
  call string_trim (label_name, label_name, ix)

elseif (e_shape%label == 'none') then
  label_name = ''

else
  call out_io (s_error$, r_name, 'BAD ELEMENT LABEL TYPE: ' // e_shape%label)
  label_name = '???'
endif 

end subroutine
