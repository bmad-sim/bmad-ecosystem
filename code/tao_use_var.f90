!+
! subroutine tao_use_var (action, var_name, locations)
!
! Veto, restore or use specified datums. range syntax: 1:34 46 58:78
!
! Input:
!   var_name  -- charatcer(*): the selected variable name or all
!   locations -- character(*): the index location expression
!
! Output:
!-

subroutine tao_use_var (action, var_name, locations)

use tao_mod

implicit none

character(*)                     :: action
character(*)                     :: var_name
character(*)                     :: locations

type (tao_v1_var_struct), pointer :: v1_ptr

logical, allocatable :: action_logic(:) !which elements do we take action on?

integer which, i, n1, n2, j
integer err_num

character(12) :: r_name = "tao_use_var"
character(200) line
character(3) select_all

logical do_all_universes
logical err, all_selected

call match_word (action, name$%use_veto_restore, which)

! Are we selecting all variables?
call str_upcase(select_all, var_name(1:3))
all_selected = .false.
if (index(select_all, 'ALL') .ne. 0) all_selected = .true.

if (all_selected) then
  !loop over all variables
  if (locations .eq. ' ') locations = "all"
  do j = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(j)
    call use_var ()
    if (err_num .eq. -1) return
  enddo
else
  ! find variable
  call tao_find_var (err, var_name, v1_ptr)
  if (err) return
  call use_var ()
  if (err_num .eq. -1) return
endif

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

!make sure v1_ptr is set properly before calling this!
subroutine use_var ()

! find locations

err_num = 0
n1 = lbound(v1_ptr%v, 1)
n2 = ubound(v1_ptr%v, 1)
allocate(action_logic(n1:n2))
line = locations
call location_decode (line, action_logic, n1, err_num) 
if (err_num .eq. -1) return

! set d%good_user based on action and action_logic

select case (which)
case (use$)
  v1_ptr%v(:)%good_user = action_logic
case (veto$)
  where (action_logic) v1_ptr%v(:)%good_user = .not. action_logic
case (restore$)
  where (action_logic) v1_ptr%v(:)%good_user = action_logic
case default
  call out_io (s_error$, r_name, &
                    "Internal error picking name$%use_veto_restore")
  err = .true.
end select


! optimizer bookkeeping and print out changes

call tao_set_var_useit_opt()
call tao_var_show_use (v1_ptr)

deallocate(action_logic)

end subroutine use_var  

end subroutine tao_use_var
