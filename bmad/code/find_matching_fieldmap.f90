!+
! Subroutine find_matching_fieldmap (file_name, ele, fm_type, match_ele, ix_field, ignore_slaves)
!
! Routine to find an equivalent fieldmap to the one given.
! Only lattice elements before ele are searched.
!
! Input:
!   file_name     -- character(*): File name associated with field to match to.
!   ele           -- ele_struct: Element holding the field to be matched.
!   fm_type       -- integer: Type of fieldmap: cartesian_map$, cylindircal_map$, or gen_grad_map$, grid_field$
!   ignore_slaves -- logical, optional: If True, ignore any multipass slaves. Default is False.
!
! Output:
!   match_ele   -- ele_struct, pointer: Pointer to element with matched field. Nullified if no match found.
!   ix_field    -- integer: index of field. For example: matching field => match_ele%cartesian_map(ix_field)
!                   Set to -1 if no match found.
!-

subroutine find_matching_fieldmap (file_name, ele, fm_type, match_ele, ix_field, ignore_slaves)

use bmad_struct

implicit none

type (ele_struct), target ::  ele
type (ele_struct), pointer :: match_ele
type (lat_struct), pointer :: lat

integer fm_type, ix_field
integer ib, ie

logical, optional :: ignore_slaves

character(*) file_name

!

lat => ele%branch%lat

do ib = 0, ele%ix_branch
  do ie = 1, lat%branch(ib)%n_ele_max

    if (ib == ele%ix_branch .and. ie == ele%ix_ele) exit
    match_ele => lat%branch(ib)%ele(ie)
    if (logic_option(.false., ignore_slaves) .and. match_ele%slave_status == multipass_slave$) cycle

    select case (fm_type)
    case (cartesian_map$)
      if (.not. associated(match_ele%cartesian_map)) cycle
      do ix_field = 1, size(match_ele%cartesian_map)
        if (match_ele%cartesian_map(ix_field)%ptr%file == file_name) return
      enddo
  
    case (cylindrical_map$)
      if (.not. associated(match_ele%cylindrical_map)) cycle
      do ix_field = 1, size(match_ele%cylindrical_map)
        if (match_ele%cylindrical_map(ix_field)%ptr%file == file_name) return
      enddo

    case (gen_grad_map$)
      if (.not. associated(match_ele%gen_grad_map)) cycle
      do ix_field = 1, size(match_ele%gen_grad_map)
        if (match_ele%gen_grad_map(ix_field)%file == file_name) return
      enddo

    case (grid_field$)
      if (.not. associated(match_ele%grid_field)) cycle
      do ix_field = 1, size(match_ele%grid_field)
        if (match_ele%grid_field(ix_field)%ptr%file == file_name) return
      enddo
    end select  

  enddo
enddo

nullify(match_ele)
ix_field = -1

end subroutine find_matching_fieldmap
