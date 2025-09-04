!+
! Subroutine tao_fixer (switch, word1, word2)
!
! Fixer commands which is one of:
!   <ele_name>
!   -active {<who>}
!   -store {<who>}
!
! Input:
!   switch    -- character(*): Action to take. One on : '-activate', '-save'.
!   word1     -- character(*): First word of command.
!   word2     -- character(*): Secton word of command.
!-

subroutine tao_fixer (switch, word1, word2)

use tao_interface, dummy => tao_fixer
use tao_command_mod, only: tao_next_switch
use fixer_mod

implicit none

type (ele_struct), pointer :: fixer
type (tao_universe_struct), pointer :: u

character(*) switch, word1, word2
character(60) ele_str
character(40) action
character(*), parameter :: r_name = 'tao_fixer'
logical is_ok, err

!

ele_str = word1
u => tao_pointer_to_universe(ele_str)
if (.not. associated(u)) return
fixer => pointer_to_ele (u%model%lat, word1)
if (.not. associated(fixer)) return
if (fixer%key /= fixer$ .and. fixer%key /= beginning_ele$) then
  call out_io(s_error$, r_name, 'Element if not a fixer nor a beginning element: ' // ele_str)
  return
endif

call tao_next_switch (switch, [character(20):: '-activate', '-on', '-save'], .false., action, err);  if (err) return

select case (action)
case ('-activate', '-on')
  if (word2 /= '') then
    call out_io(s_error$, r_name, 'Extra stuff on line: ' // word2)
    return
  endif

  call set_active_fixer(fixer, .true.)
  u%calc%lattice = .true.

case ('-save')
  is_ok = transfer_fixer_params(fixer, .true., word2)

case default
  call out_io(s_error$, r_name, 'Switch not recognized: ' // action)
end select

end subroutine
