!+
! subroutine tao_use_data (do_all_universes, action, data_name, locations)
!
! Veto, restore or use specified datums. range syntax is just like
!    indexing in fortran: 1:34, 46, 58:78
!
! Input:
!   do_all_universes -- Logical: Apply to all universes?
!                         if not just use s%u(s%global%u_view)
!   action	         -- character(*): veto, use or restore
!   data_name        -- charatcer(*): the selected data name
!   data_name        -- character(*): the selected data name
!   locations        -- character(*): the index location expression
!
! Output:
!-

subroutine tao_use_data (do_all_universes, action, data_name, locations)

use tao_mod

implicit none

character(*)                :: action
character(*)                :: data_name
character(*)                :: locations

type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr

logical do_all_universes
logical, allocatable :: action_logic(:) !which elements do we take action on?

integer which, i, iu, n1, n2
logical err

integer err_num

character(12) :: r_name = "tao_use_data"
character(200) line

! decipher action

call match_word (action, name$%use_veto_restore, which)

! loop over the universes to do.

if (do_all_universes) then
  do iu = 1, size(s%u)
    call use_data (s%u(i))
  enddo
else
  call use_data (s%u(s%global%u_view))
endif

!----------------------------------------------------------------
contains

subroutine use_data (u)

type (tao_universe_struct) u

! find data name and name

call tao_find_data (err, u, data_name, d2_ptr, d1_ptr)
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

end subroutine

!----------------------------------------------------------------
! contains

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
