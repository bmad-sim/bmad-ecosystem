!+
! subroutine tao_use_data (action, data_type, locations)
!
! Veto, restore or use specified datums. range syntax is just like
!    indexing in fortran: 1:34, 46, 58:78
!
! Input:
!   action	         -- character(*): veto, use or restore
!   data_type        -- character(*): the selected data name
!   locations        -- character(*): the index location expression
!
! Output:
!-

subroutine tao_use_data (action, data_type, locations)

use tao_mod

implicit none

type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr

character(*) :: action
character(*) :: data_type
character(*) :: locations

logical, allocatable :: action_logic(:) !which elements do we take action on?
logical, automatic :: picked(size(s%u))
logical err

integer which, i, iu, n1, n2, ix1, ix2
integer err_num

character(12) :: r_name = "tao_use_data"
character(16) d_name
character(200) line

! decipher action

call match_word (action, name$%use_veto_restore, which)

! loop over the universes to do.

call tao_pick_universe (data_type, d_name, picked, err)
if (err) return

do iu = 1, size(s%u)

  ! find data name and name

  if (.not. picked(iu)) cycle
  call tao_find_data (err, s%u(iu), d_name, d2_ptr, d1_ptr)
  if (err) return

  ! find locations

  n1 = lbound(d2_ptr%d1(1)%d, 1)
  n2 = ubound(d2_ptr%d1(1)%d, 1)
  allocate(action_logic(n1:n2))
  call location_decode (locations, action_logic, n1, err_num) 
  if (err_num == -1) return

  ! set d%good_user based on action and action_logic

  if (associated(d1_ptr)) then
    call use (d1_ptr)
  else
    do i = 1, size(d2_ptr%d1)
      call use (d2_ptr%d1(i))
      if (err) return
    enddo
  endif

  ! Optimizer bookkeeping and Print out changes.

  call tao_set_data_useit_opt()
  call tao_data_show_use (d2_ptr)

  deallocate(action_logic)

enddo

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

subroutine use (d1)

type (tao_d1_data_struct) d1
integer i

!

select case (which)
case (use$)
  d1%d(:)%good_user = action_logic
case (veto$)
  where (action_logic) d1%d(:)%good_user = .not. action_logic
case (restore$)
  where (action_logic) d1%d(:)%good_user = action_logic
case default
  call out_io (s_error$, r_name, &
                    "Internal error picking name$%use_veto_restore")
  err = .true.
end select

end subroutine

end subroutine tao_use_data
