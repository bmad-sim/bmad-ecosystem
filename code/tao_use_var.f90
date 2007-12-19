!+
! subroutine tao_use_var (action, var_name)
!
! Veto, restore or use specified datums. 
!
! Input:
!   action    -- Character(*): 'use', 'veto', or 'restore'
!   var_name  -- charatcer(*): the selected variable name or all
!-

subroutine tao_use_var (action, var_name)

use tao_mod

implicit none

character(*) :: action
character(*) :: var_name

type (tao_var_array_struct), allocatable, save :: var(:)

logical, allocatable :: action_logic(:) !which elements do we take action on?

integer which, i

character(12) :: r_name = "tao_use_var"

logical err

! Use, veto or restore?

call match_word (action, name$%use_veto_restore, which)

! If "use" then must veto everything first

if (which == use$) then
  call tao_find_var (err, var_name, v_array = var, all_elements = .true.)
  if (err) return
  do i = 1, size(var)
    var(i)%v%good_user = .false.
  enddo
endif

! now do the set.

call tao_find_var (err, var_name, v_array = var)
if (err) return

do i = 1, size(var)
  select case (which)
  case (use$, restore$)
    var(i)%v%good_user = .true.
  case (veto$)
    var(i)%v%good_user = .false.
  case default
    call out_io (s_error$, r_name, &
                    "Internal error picking name$%use_veto_restore")
    err = .true.
  end select
enddo

! optimizer bookkeeping and print out changes

call tao_set_var_useit_opt()

call tao_var_show_use (var(1)%v%v1)
do i = 2, size(var)
  if (.not. associated (var(i)%v%v1, var(i-1)%v%v1)) &
                          call tao_var_show_use (var(i)%v%v1)
enddo

end subroutine tao_use_var
