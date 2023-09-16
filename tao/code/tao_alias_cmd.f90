!+
! Subroutine tao_alias_cmd (alias, string)
!
! Routine define an alias command.
! 
! Input:
!   alias  -- Name of the tao command file.
!   string(9) -- Command file arguments.
!
! Output:
!   
!-

subroutine tao_alias_cmd (alias, string)

use tao_interface, dummy => tao_alias_cmd

implicit none

character(*) alias
character(*) string
character(20) :: r_name = 'tao_alias_cmd'

integer i

! print aliases?

if (alias == ' ') then
  do i = 1, s%com%n_alias
    call out_io (s_blank$, r_name, s%com%alias(i)%name // s%com%alias(i)%expanded_str)
  enddo
  return
endif

! Check to see if alias already defined. If so just overwrite it.

do i = 1, s%com%n_alias
  if (s%com%alias(i)%name == alias) then
    s%com%alias(i)%expanded_str = string
    return
  endif
enddo

if (s%com%n_alias >= size(s%com%alias)) then
  call out_io (s_error$, r_name, 'NUMBER OF ALIASES TOO LARGE!')
  return
endif

s%com%n_alias = s%com%n_alias + 1
s%com%alias(s%com%n_alias)%name = alias
s%com%alias(s%com%n_alias)%expanded_str = string

end subroutine 
