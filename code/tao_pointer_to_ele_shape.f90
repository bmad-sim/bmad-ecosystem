!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape_min) result (e_shape)
!
! Routine to return the shape associated with a lattice element
!
! Input:
!   ix_uni        -- integer: Universe index.
!   ele           -- ele_struct: Lattice element.
!   ele_shape(:)  -- tao_ele_shape_struct: Array of shapes to search.
!   ix_shape_min  -- integer, optional: Index of minimum ele_shape(:) index to start search from. Default is 1.
!
! Output:
!   e_shape       -- tao_ele_shape_struct, pointer: Associated shape. 
!                       Nullified if there is no associated shape.
!   dat_var_name  -- character(*), optional: Name of datum or variable associated with e_shape. 
!                       Will be set to "" if there is no associated datum or variable.
!   dat_var_value -- real(rp), optional: Value of datum or variable associated with e_shape.
!                       Will be set to zero if there is no associated datum or variable.
!   ix_shape_min  -- integer, optional: Ele_shape(:) index to start next search if multiple shapes are associated with ele.
!-

function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape_min) result (e_shape)

use tao_interface, dummy => tao_pointer_to_ele_shape

implicit none

type (ele_struct) ele

type (tao_ele_shape_struct), target :: ele_shape(:)
type (tao_ele_shape_struct), pointer :: e_shape, es
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)
type (tao_data_struct), pointer :: datum
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_var_struct), pointer :: var
type (tao_real_pointer_struct), allocatable :: re_array(:)

real(rp), optional :: dat_var_value

integer, optional :: ix_shape_min
integer ix_uni
integer j, j2, k, n_ele_track

character(*), optional :: dat_var_name
character(*), parameter :: r_name = 'tao_pointer_to_ele_shape'

logical err

!

nullify(e_shape)
if (present(dat_var_name)) dat_var_name = ''
if (present(dat_var_value)) dat_var_value = 0

if (ele%lord_status == group_lord$) return
if (ele%lord_status == overlay_lord$) return
if (ele%slave_status == super_slave$) return

do k = integer_option(1, ix_shape_min), size(ele_shape)
  es => ele_shape(k)
  if (present(ix_shape_min)) ix_shape_min = k + 1
  if (.not. es%draw) cycle

  ! Data

  if (index(es%ele_id, 'dat::') /= 0) then
    call out_io (s_error$, r_name, 'SHAPE USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // es%ele_id)
    call err_exit
  endif

  if (es%ele_id(1:6) == 'data::') then
    call tao_find_data (err, es%ele_id, d_array = d_array, log_array = logic_array, re_array = re_array)
    if (err) cycle
    do j = 1, size(d_array)
      datum => d_array(j)%d
      if (datum%d1%d2%ix_universe /= ix_uni) cycle
      if (size(logic_array) /= 0) then
        if (.not. logic_array(j)%l) cycle
      endif
      if (ele%ix_branch /= datum%ix_branch .or. ele%ix_ele /=  datum%ix_ele) cycle
      e_shape => es
      if (present(dat_var_name)) dat_var_name = tao_datum_name(datum)
      if (present(dat_var_value)) then
        if (size(re_array) > 0) then
          dat_var_value = re_array(j)%r
        else
          dat_var_value = datum%model_value
        endif
      endif
      return
    enddo
    cycle
  endif

  ! Variables

  if (es%ele_id(1:5) == 'var::') then
    call tao_find_var (err, es%ele_id, v_array = v_array, log_array = logic_array, re_array = re_array)
    if (err) cycle

    do j = 1, size(v_array)
      var => v_array(j)%v
      if (size(logic_array) /= 0) then
        if (.not. logic_array(j)%l) cycle
      endif
      do j2 = 1, size(var%slave)
        if (var%slave(j2)%ix_uni /= ix_uni) cycle
        if (var%ele_name == 'PARTICLE_START') then
          if (ele%ix_branch /= 0 .or. ele%ix_ele /= 0) cycle
        else
          if (ele%ix_branch /= var%slave(j2)%ix_branch .or. ele%ix_ele /=  var%slave(j2)%ix_ele) cycle
        endif
        e_shape => es
        if (present(dat_var_name)) dat_var_name = tao_var1_name(var)
        if (present(dat_var_value)) then
          if (size(re_array) > 0) then
            dat_var_value = re_array(j)%r
          else
            dat_var_value = var%slave(j2)%model_value
          endif
        endif
        return
      enddo
    enddo
    cycle
  endif

  ! All else

  if (es%ele_id == '') cycle
  if (es%ix_ele_key == -1) cycle
  if (es%ix_ele_key /= 0 .and. es%ix_ele_key /= ele%key) cycle
  if (.not. match_wild(ele%name, es%name_ele)) cycle

  e_shape => es
  return
enddo

end function tao_pointer_to_ele_shape 
