!+
! Subroutine transfer_fieldmap (ele_in, ele_out, who)
!
! Subroutine to transfer the field info from one element to another.
! In the end will have:
!     ele_out%cartesian_map    => ele_in%cartesian_map
!     ele_out%cylindrical_map  => ele_in%cylindrical_map
!     ele_out%grid_field       => ele_in%grid_field
!
! Input:
!   ele_in -- Ele_struct: Input element.
!   who    -- integer: Possibilities are: all$, cartesian_map$, cylindrical_map$, or grid_field$
!
! Output:
!   ele_out -- Ele_struct: Output element.
!-

subroutine transfer_fieldmap (ele_in, ele_out, who)

use bmad_routine_interface, dummy => transfer_fieldmap

implicit none

type (ele_struct) :: ele_in, ele_out

integer who
integer i, nm

! Cartesian_map

if (who == all$ .or. who == cartesian_map$) then
  if (associated(ele_in%cartesian_map) .and. associated(ele_out%cartesian_map)) then
    if (size(ele_in%cartesian_map) /= size(ele_out%cartesian_map)) then
      call unlink_fieldmap (cartesian_map = ele_out%cartesian_map)
      nm = size(ele_in%cartesian_map)
      allocate (ele_out%cartesian_map(nm))
      do i = 1, nm
        ele_out%cartesian_map(i) = ele_in%cartesian_map(i)
        ele_out%cartesian_map(i)%ptr%n_link = ele_out%cartesian_map(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%cartesian_map)
        if (associated(ele_out%cartesian_map(i)%ptr, ele_in%cartesian_map(i)%ptr)) then
          ele_out%cartesian_map(i) = ele_in%cartesian_map(i) ! Make sure same info
        else
          ele_out%cartesian_map(i)%ptr%n_link = ele_out%cartesian_map(i)%ptr%n_link - 1
          if (ele_out%cartesian_map(i)%ptr%n_link == 0) deallocate (ele_out%cartesian_map(i)%ptr)
          ele_out%cartesian_map(i) = ele_in%cartesian_map(i)
          ele_out%cartesian_map(i)%ptr%n_link = ele_out%cartesian_map(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%cartesian_map) .and. .not. associated(ele_out%cartesian_map)) then
    nm = size(ele_in%cartesian_map)
    allocate (ele_out%cartesian_map(nm))
    do i = 1, nm
      ele_out%cartesian_map(i) = ele_in%cartesian_map(i)
      ele_out%cartesian_map(i)%ptr%n_link = ele_in%cartesian_map(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%cartesian_map) .and. associated(ele_out%cartesian_map)) then
    call unlink_fieldmap (cartesian_map = ele_out%cartesian_map)
  endif
endif

! Cylindrical_map

if (who == all$ .or. who == cylindrical_map$) then
  if (associated(ele_in%cylindrical_map) .and. associated(ele_out%cylindrical_map)) then
    if (size(ele_in%cylindrical_map) /= size(ele_out%cylindrical_map)) then
      call unlink_fieldmap (cylindrical_map = ele_out%cylindrical_map)
      nm = size(ele_in%cylindrical_map)
      allocate (ele_out%cylindrical_map(nm))
      do i = 1, nm
        ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
        ele_out%cylindrical_map(i)%ptr%n_link = ele_out%cylindrical_map(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%cylindrical_map)
        if (associated(ele_out%cylindrical_map(i)%ptr, ele_in%cylindrical_map(i)%ptr)) then
          ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
        else
          ele_out%cylindrical_map(i)%ptr%n_link = ele_out%cylindrical_map(i)%ptr%n_link - 1
          if (ele_out%cylindrical_map(i)%ptr%n_link == 0) deallocate (ele_out%cylindrical_map(i)%ptr)
          ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
          ele_out%cylindrical_map(i)%ptr%n_link = ele_out%cylindrical_map(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%cylindrical_map) .and. .not. associated(ele_out%cylindrical_map)) then
    nm = size(ele_in%cylindrical_map)
    allocate (ele_out%cylindrical_map(nm))
    do i = 1, nm
      ele_out%cylindrical_map(i) = ele_in%cylindrical_map(i)
      ele_out%cylindrical_map(i)%ptr%n_link = ele_in%cylindrical_map(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%cylindrical_map) .and. associated(ele_out%cylindrical_map)) then
    call unlink_fieldmap (cylindrical_map = ele_out%cylindrical_map)
  endif
endif

! Gen_grad

if (who == all$ .or. who == gen_grad_map$) then
  if (associated(ele_in%gen_grad_map) .and. associated(ele_out%gen_grad_map)) then
    if (size(ele_in%gen_grad_map) /= size(ele_out%gen_grad_map)) then
      call unlink_fieldmap (gen_grad_map = ele_out%gen_grad_map)
      nm = size(ele_in%gen_grad_map)
      allocate (ele_out%gen_grad_map(nm))
      do i = 1, nm
        ele_out%gen_grad_map(i) = ele_in%gen_grad_map(i)
      enddo

    else
      do i = 1, size(ele_in%gen_grad_map)
        ele_out%gen_grad_map(i) = ele_in%gen_grad_map(i)
      enddo
    endif

  elseif (associated(ele_in%gen_grad_map) .and. .not. associated(ele_out%gen_grad_map)) then
    nm = size(ele_in%gen_grad_map)
    allocate (ele_out%gen_grad_map(nm))
    do i = 1, nm
      ele_out%gen_grad_map(i) = ele_in%gen_grad_map(i)
    enddo

  elseif (.not. associated(ele_in%gen_grad_map) .and. associated(ele_out%gen_grad_map)) then
    call unlink_fieldmap (gen_grad_map = ele_out%gen_grad_map)
  endif
endif

! Grid_field

if (who == all$ .or. who == grid_field$) then
  if (associated(ele_in%grid_field) .and. associated(ele_out%grid_field)) then
    if (size(ele_in%grid_field) /= size(ele_out%grid_field)) then
      call unlink_fieldmap (grid_field = ele_out%grid_field)
      nm = size(ele_in%grid_field)
      allocate (ele_out%grid_field(nm))
      do i = 1, nm
        ele_out%grid_field(i) = ele_in%grid_field(i)
        ele_out%grid_field(i)%ptr%n_link = ele_out%grid_field(i)%ptr%n_link + 1
      enddo

    else
      do i = 1, size(ele_in%grid_field)
        if (associated(ele_out%grid_field(i)%ptr, ele_in%grid_field(i)%ptr)) then
          ele_out%grid_field(i) = ele_in%grid_field(i)
        else
          ele_out%grid_field(i)%ptr%n_link = ele_out%grid_field(i)%ptr%n_link - 1
          if (ele_out%grid_field(i)%ptr%n_link == 0) deallocate (ele_out%grid_field(i)%ptr)
          ele_out%grid_field(i) = ele_in%grid_field(i)
          ele_out%grid_field(i)%ptr%n_link = ele_out%grid_field(i)%ptr%n_link + 1
        endif
      enddo
    endif

  elseif (associated(ele_in%grid_field) .and. .not. associated(ele_out%grid_field)) then
    nm = size(ele_in%grid_field)
    allocate (ele_out%grid_field(nm))
    do i = 1, nm
      ele_out%grid_field(i) = ele_in%grid_field(i)
      ele_out%grid_field(i)%ptr%n_link = ele_in%grid_field(i)%ptr%n_link + 1
    enddo

  elseif (.not. associated(ele_in%grid_field) .and. associated(ele_out%grid_field)) then
    call unlink_fieldmap (grid_field = ele_out%grid_field)
  endif
endif

end subroutine transfer_fieldmap

