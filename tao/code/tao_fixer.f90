!+
! Subroutine tao_fixer (switch, word1, word2)
!
! Input:
!   switch    -- character(*): Action to take. One on : 'activate', 'save', 'write'.
!   word1     -- character(*): First word of command.
!   word2     -- character(*): Secton word of command.
!-

subroutine tao_fixer (switch, word1, word2)

use tao_interface, dummy => tao_fixer
use tao_command_mod, only: tao_next_switch
use fixer_mod

implicit none

type (ele_struct), pointer :: fixer
type (ele_struct) :: dflt_fixer
type (tao_universe_struct), pointer :: u
type (ele_attribute_struct) attrib
type (branch_struct), pointer :: branch

character(*) switch, word1, word2
character(40) action, name, str_val
character(60) ele_str
character(200) file_name, line1, line2
character(*), parameter :: r_name = 'tao_fixer'

real(rp) val
integer j, iu
logical is_ok, err, is_default

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
branch => fixer%branch

call tao_next_switch (switch, [character(20):: 'activate', 'on', 'save', 'write'], &
                                        .false., action, err);  if (err) return

select case (action)
case ('activate', 'on')
  if (word2 /= '') then
    call out_io(s_error$, r_name, 'Extra stuff on line: ' // word2)
    return
  endif

  call set_active_fixer(fixer, .true.)
  u%calc%lattice = .true.

case ('save')
  is_ok = transfer_fixer_params(fixer, branch%particle_start, .true., word2)

case ('write')
  file_name = word2
  if (file_name == '') file_name = trim(fixer%name) // '.bmad'
  call init_ele(dflt_fixer, fixer$)

  iu = lunget()
  open (iu, file = file_name)
  write (line1, '(2a)') trim(fixer%name), ': fixer'

  do j = 1, num_ele_attrib$
    attrib = attribute_info(fixer, j)
    if (attrib%state /= is_free$) cycle
    val = fixer%value(j)
    if (val == dflt_fixer%value(j)) cycle
    if (attribute_type(attrib%name) == is_switch$) then
      name = switch_attrib_value_name (attrib%name, val, fixer, is_default)
      if (is_default) cycle
    endif

    write (iu, '(2a)') trim(line1), ','

    select case (attribute_type(attrib%name))
    case (is_logical$)
      write (line1, '(3a, l1)') '      ', trim(attrib%name), ' = ', (val /= 0)
    case (is_integer$)
      write (line1, '(3a, i0)') '      ', trim(attrib%name), ' = ', int(val)
    case (is_real$)
      write(str_val, '(es25.17e3)') val
      write (line1, '(4a)') '      ', trim(attrib%name), ' = ', str_val
    case (is_switch$)
      write (line1, '(4a)') '      ', trim(attrib%name), ' = ', name
    end select
  enddo

  write (iu, '(2a)') trim(line1)
  close(iu)
  call out_io(s_info$, r_name, 'Created lattice file: ' // file_name)

case default
  call out_io(s_error$, r_name, 'Switch not recognized: ' // action)
end select

end subroutine
