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
type (tao_v1_var_array_struct), allocatable, save :: v1(:)

logical, allocatable :: action_logic(:) !which elements do we take action on?

integer which, i

character(12) :: r_name = "tao_use_var"

logical err

! Use, veto or restore?

call match_word (action, name$%use_veto_restore, which)

call tao_find_var (err, var_name, v1_array = v1, v_array = var)
if (err) return

! If "use" then must veto everything first


if (which == use$) then
  do i = 1, size(v1)
    v1(i)%v1%v%good_user = .false.
  enddo
endif

! now do the set.

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

do i = 1, size(s%v1_var)
  if (s%v1_var(i)%name == ' ') cycle
  call tao_var_show_use (s%v1_var(i))
enddo

end subroutine tao_use_var
