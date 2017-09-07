!+
! Subroutine tao_read_cmd (who, file_name)
!
! Routine to read in stuff.
!
! Input:
!   who       -- Character(*): Must be 'lattice'.
!   file_name -- Character(*): Input file name.
!-

subroutine tao_read_cmd (who, file_name)

use tao_mod, dummy => tao_read_cmd
use madx_ptc_module, only: m_u, m_t, read_universe_pointed

implicit none

type (tao_universe_struct), pointer :: u

character(*) who, file_name
character(20) action
character(20) :: names(2) = ['lattice', 'ptc    ']
character(*), parameter :: r_name = 'tao_read_cmd'

integer i, j, ix, iu, nd, ii
logical err

!

call string_trim (who, action, ix)
call match_word (action, names, ix)
if (ix == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED "WHAT": ' // action)
  return
elseif (ix < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS "WHAT": ' // action)
  return
endif
action = names(ix)

select case (action)

!-------------------------
! lattice

case ('lattice')

  u => tao_pointer_to_universe(-1)
  call bmad_parser2 (file_name, u%model%lat)
  call tao_var_repoint()
  u%calc%lattice = .true.

  do i = 0, ubound(u%model%lat%branch, 1)
    if (u%model%lat%branch(i)%n_ele_track /= u%design%lat%branch(i)%n_ele_track .or. &
        u%model%lat%branch(i)%n_ele_max /= u%design%lat%branch(i)%n_ele_max) then
      call out_io (s_fatal$, r_name, &
              'IT IS FORBIDDEN TO USE THE "read Lattice" COMMAND TO MODIFY THE NUMBER OF LATTICE ELEMENTS.', &
              'WILL STOP HERE.')
      stop
    endif
  enddo

case ('ptc')
  call read_universe_pointed (m_u, m_t, file_name)

end select

end subroutine
