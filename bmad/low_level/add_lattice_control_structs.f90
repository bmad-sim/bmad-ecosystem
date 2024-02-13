!+
! Subroutine add_lattice_control_structs (ele, n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field, add_at_end)
! 
! Subroutine to adjust the control structure of a lat so that extra control elements can be added.
! Note: Control struct arrays cannot be reduced by this routine.
!
! Input:
!   ele               -- ele_struct: Lord or slave element that needs extra control elements.
!   n_add_slave       -- integer, optional: Number of field slaves to add to lord. Default is zero.
!   n_add_lord        -- integer, optional: Number of field lords to add to slave. Default is zero.
!   n_add_slave_field -- integer, optional: Number of field slaves to add to lord. Default is zero.
!   n_add_lord_field  -- integer, optional: Number of field lords to add to slave. Default is zero.
!   add_at_end        -- logical, optional: Used when n_add_slave or n_add_slave_field is non-zero.
!                          If True then new space is added at the end of the array.
!                          If False then new space is added at the front of the array.
!                          Default is True.
!
! Output:
!   lat -- lat_struct: Lat with control structure fixed.
!-

subroutine add_lattice_control_structs (ele, n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field, add_at_end)

use bmad_interface, except_dummy => add_lattice_control_structs

implicit none

type (lat_struct), pointer :: lat
type (ele_struct) ele
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele2

integer, optional ::  n_add_slave, n_add_lord, n_add_slave_field, n_add_lord_field
integer nadd_slave, nadd_lord, n_con, i1, i2, i2_field, n_con2, n_ic, n_ic2, ib, ie, nadd_field
integer i1_field

logical, optional :: add_at_end

character(40), parameter :: r_name = 'add_lattice_control_structs'

! Fix slave bookkeeping in lord

lat => ele%branch%lat
nadd_slave = integer_option(0, n_add_slave)
nadd_field = integer_option(0, n_add_slave_field)

if (nadd_slave < 0 .or. nadd_field < 0) then
  call out_io (s_fatal$, r_name, 'Number of slaves to add is negative! ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Add field slaves

if (nadd_field > 0) then

  n_con = lat%n_control_max
  n_con2 = lat%n_control_max + nadd_field
  if (n_con2 > size(lat%control)) call reallocate_control(lat, nint(1.2*n_con2)+10)
  lat%n_control_max = n_con2

  ! If no existing slaves of this lord then just put the new lat%control elements
  ! at the end of the array.

  if (ele%ix1_slave < 1) then
    ele%ix1_slave = n_con + 1
    lat%control(n_con+1:n_con+nadd_field)%lord = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)

  ! Else we need to make room in lat%control for the new slaves by moving
  ! a slice of lat%control.

  else
    i2 = ele%ix1_slave + ele%n_slave - 1
    i1_field = i2 + 1
    i2_field = i2 + ele%n_slave_field

    if (logic_option(.true., add_at_end)) then
      lat%control(i2_field+nadd_field+1:n_con+nadd_field) = lat%control(i2_field+1:n_con)
      lat%control(i2_field+1:i2_field+nadd_field)%lord  = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)
      lat%control(i2_field+1:i2_field+nadd_field)%slave = lat_ele_loc_struct()
      lat%control(i2_field+1:i2_field+nadd_field)%ix_attrib = 0
      lat%control(i2_field+1:i2_field+nadd_field)%attribute = ''
      where (lat%ic > i2_field) lat%ic = lat%ic + nadd_field

    else
      lat%control(i1_field+nadd_field:n_con+nadd_field) = lat%control(i1_field:n_con)
      lat%control(i1_field:i1_field+nadd_field-1)%lord  = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)
      lat%control(i1_field:i1_field+nadd_field-1)%slave = lat_ele_loc_struct()
      lat%control(i1_field:i1_field+nadd_field-1)%ix_attrib = 0
      lat%control(i1_field:i1_field+nadd_field-1)%attribute = ''
      where (lat%ic >= i1_field) lat%ic = lat%ic + nadd_field
    endif

    do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      do ie = 1, branch%n_ele_max
        ele2 => branch%ele(ie)
        if (ele2%ix1_slave <= ele%ix1_slave) cycle
        ele2%ix1_slave = ele2%ix1_slave + nadd_field
      enddo
    enddo
  endif

  ele%n_slave_field = ele%n_slave_field + nadd_field
endif

! Add regular (non-field) slaves to lord

if (nadd_slave > 0) then

  n_con = lat%n_control_max
  n_con2 = lat%n_control_max + nadd_slave
  if (n_con2 > size(lat%control)) call reallocate_control(lat, nint(1.2*n_con2)+10)
  lat%n_control_max = n_con2

  ! If no existing slaves of this lord then just put the new lat%control elements
  ! at the end of the array.

  if (ele%ix1_slave < 1) then  
    ele%ix1_slave = n_con + 1
    lat%control(n_con+1:n_con+nadd_slave)%lord = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)

  ! Else we need to make room in lat%control for the new slaves by moving
  ! a slice of lat%control.

  else
    i2 = ele%ix1_slave + ele%n_slave - 1

    if (logic_option(.true., add_at_end)) then
      lat%control(i2+1+nadd_slave:n_con+nadd_slave) = lat%control(i2+1:n_con)
      lat%control(i2+1:i2+nadd_slave)%lord = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)
      lat%control(i2+1:i2+nadd_slave)%slave = lat_ele_loc_struct()
      lat%control(i2+1:i2+nadd_slave)%ix_attrib = 0
      lat%control(i2+1:i2+nadd_slave)%attribute = ''
      where (lat%ic > i2) lat%ic = lat%ic + nadd_slave

    else
      i1 = ele%ix1_slave
      lat%control(i1+nadd_slave:n_con+nadd_slave) = lat%control(i1:n_con)
      lat%control(i1:i1+nadd_slave-1)%lord = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)
      lat%control(i1:i1+nadd_slave-1)%slave = lat_ele_loc_struct()
      lat%control(i1:i1+nadd_slave-1)%ix_attrib = 0
      lat%control(i1:i1+nadd_slave-1)%attribute = ''
      where (lat%ic >= i1) lat%ic = lat%ic + nadd_slave
    endif

    do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      do ie = 1, branch%n_ele_max
        ele2 => branch%ele(ie)
        if (ele2%ix1_slave <= ele%ix1_slave) cycle
        ele2%ix1_slave = ele2%ix1_slave + nadd_slave
      enddo
    enddo
  endif

  ele%n_slave = ele%n_slave + nadd_slave      
endif

!--------------------------------------------------------------------------------                                    
! Fix lord bookkeeping in slave

nadd_lord   = integer_option(0, n_add_lord)    
nadd_field = integer_option(0, n_add_lord_field)

if (nadd_lord < 0 .or. nadd_field < 0) then
  call out_io (s_fatal$, r_name, 'Number of lords to add is negative! ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Add field lords

if (nadd_field > 0) then
  n_ic = lat%n_ic_max
  n_ic2 = n_ic + nadd_field
  lat%n_ic_max = n_ic2

  if (n_ic2 > size(lat%ic)) call reallocate_control(lat, nint(1.2*(n_ic2)) + 10)

  if (ele%ic1_lord < 1) then
    ele%ic1_lord = n_ic + 1
  else
    i2 = ele%ic1_lord + ele%n_lord - 1
    i2_field = i2 + ele%n_lord_field
    lat%ic(i2_field+1+nadd_field:n_ic2) = lat%ic(i2_field+1:n_ic)
    lat%ic(i2_field+1:i2_field+nadd_field) = 0
    do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      do ie = 1, branch%n_ele_max
        ele2 => branch%ele(ie)
        if (ele2%ic1_lord <= ele%ic1_lord) cycle
        ele2%ic1_lord = ele2%ic1_lord + nadd_field
      enddo
    enddo
  endif

  ele%n_lord_field = ele%n_lord_field + nadd_field
endif

! Add regular (non-field) lords to slave

if (nadd_lord > 0) then
  n_ic = lat%n_ic_max
  n_ic2 = n_ic + nadd_lord
  lat%n_ic_max = n_ic2

  if (n_ic2 > size(lat%ic)) call reallocate_control(lat, nint(1.2*(n_ic2)) + 10)

  if (ele%ic1_lord < 1) then
    ele%ic1_lord = n_ic + 1
  else
    i2 = ele%ic1_lord + ele%n_lord - 1
    lat%ic(i2+1+nadd_lord:n_ic2) = lat%ic(i2+1:n_ic)
    lat%ic(i2+1:i2+nadd_lord) = 0
    do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      do ie = 1, branch%n_ele_max
        ele2 => branch%ele(ie)
        if (ele2%ic1_lord <= ele%ic1_lord) cycle
         ele2%ic1_lord = ele2%ic1_lord + nadd_lord
      enddo
    enddo
  endif

  ele%n_lord = ele%n_lord + nadd_lord
endif

end subroutine
