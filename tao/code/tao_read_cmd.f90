!+
! Subroutine tao_read_cmd (who, unis, file_name, silent)
!
! Routine to read in stuff.
!
! Input:
!   who       -- character(*): Must be 'lattice'.
!   unis      -- character(*): Universes to apply to
!   file_name -- character(*): Input file name.
!   silent    -- logical: Silent 
!-

subroutine tao_read_cmd (who, unis, file_name, silent)

use tao_interface, dummy => tao_read_cmd
use madx_ptc_module, only: m_u, m_t, read_universe_pointed

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_struct), pointer :: var

character(*) who, unis, file_name
character(20) action
character(20) :: names(2) = ['lattice', 'ptc    ']
character(*), parameter :: r_name = 'tao_read_cmd'

integer i, j, iv, is, ix, iu, nd, ii
logical silent, err
logical, allocatable :: u_pick(:)

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

  call tao_pick_universe(trim(unis) // '@', action, u_pick, err);  if (err) return

  do iu = 1, ubound(s%u, 1)
    if (.not. u_pick(iu)) cycle
    u => s%u(iu)

    call bmad_parser2 (file_name, u%model%lat, err_flag = err)
    u%calc%lattice = .true.
    if (err) call tao_abort_command_file()

    ! If the lattice is perturbed in terms of the number of elements, tao_var_repoint will fail. 
    ! The reason is that tao_var_repoint will not be able to find the variable slave elements if they have moved.
    ! Also Tao routines are allowed to assume that elements in the design and model lattices are commensurate.
    ! Therefore, it is forbidden for any Tao command (for example the "read lattice" command) to modify element positions.

    if (.not. silent) then
      do i = 0, ubound(u%model%lat%branch, 1)
        if (u%model%lat%branch(i)%n_ele_track /= u%design%lat%branch(i)%n_ele_track .or. &
            u%model%lat%branch(i)%n_ele_max /= u%design%lat%branch(i)%n_ele_max) then
          call out_io (s_fatal$, r_name, &
                  'IT IS FORBIDDEN TO USE THE "read Lattice" COMMAND TO MODIFY THE NUMBER OF LATTICE ELEMENTS.', &
                  'WILL STOP HERE.')
          stop
        endif
      enddo
    endif

  enddo

  ! Check for consistancy of variable slaves.

  call tao_var_repoint()

  do iv = 1, s%n_var_used
    var => s%var(iv)
    do is = 2, size(var%slave)
      if (var%slave(is)%model_value == var%slave(1)%model_value) cycle
      call out_io (s_error$, r_name, 'TAO VARIABLE: ' // tao_var1_name(var), &
              'WHICH HAS MULTIPLE SLAVE PARAMETERS NOW DOES NOT HAVE ALL SLAVE PARAMETER VALUES THE SAME.', &
              'THIS CAN CAUSE BEHAVIOR MOST STRANGE.', &
              'NOTE: "set var <var_name>|model = <value>" CAN BE USED TO RECTIFY THE SITUATION.')
      exit
    enddo
  enddo

! This is experimental!
case ('ptc')
  call read_universe_pointed (m_u, m_t, file_name)

end select

end subroutine
