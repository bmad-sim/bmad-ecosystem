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

use tao_interface, dummy => tao_read_cmd
use madx_ptc_module, only: m_u, m_t, read_universe_pointed

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_struct), pointer :: var

character(*) who, file_name
character(20) action
character(20) :: names(2) = ['lattice', 'ptc    ']
character(*), parameter :: r_name = 'tao_read_cmd'

integer i, j, iv, is, ix, iu, nd, ii
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

  ! If the lattice is perturbed in terms of the number of elements, tao_var_repoint will fail. 
  ! The reason is that tao_var_repoint will not be able to find the variable slave elements if they have moved.
  ! Also Tao routines are allowed to assume that elements in the design and model lattices are commensurate.
  ! Therefore, it is forbidden for any Tao command (for example the "read lattice" command) to modify element positions.

  do i = 0, ubound(u%model%lat%branch, 1)
    if (u%model%lat%branch(i)%n_ele_track /= u%design%lat%branch(i)%n_ele_track .or. &
        u%model%lat%branch(i)%n_ele_max /= u%design%lat%branch(i)%n_ele_max) then
      call out_io (s_fatal$, r_name, &
              'IT IS FORBIDDEN TO USE THE "read Lattice" COMMAND TO MODIFY THE NUMBER OF LATTICE ELEMENTS.', &
              'WILL STOP HERE.')
      stop
    endif
  enddo

  ! Check for consistancy of variable slaves.

  do iv = 1, s%n_var_used
    var => s%var(iv)
    do is = 2, size(var%slave)
      if (var%slave(is)%model_value == var%slave(1)%model_value) cycle
      call out_io (s_error$, r_name, 'TAO VARIABLE: ' // tao_var1_name(var), &
              'WHICH HAS MULTIPLE SLAVE PARAMETERS NOW DOES NOT HAVE ALL SLAVE PARAMETER VALUES THE SAME.', &
              'THIS WILL CAUSE STRANGE BEHAVIOR. RECOMMENDATION: USE THE "set ele -update" COMMAND INSTEAD.')
      exit
    enddo
  enddo

case ('ptc')
  call read_universe_pointed (m_u, m_t, file_name)

end select

end subroutine
