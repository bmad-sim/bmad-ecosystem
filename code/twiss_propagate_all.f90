!+
! Subroutine twiss_propagate_all (lat, ix_branch)
!
! Subroutine to propagate the twiss parameters from the start to the end.
!
! Note: lat%ele(:)%a Twiss parameters are associated with the "A" mode and
! the lat%ele(:)%b Twiss parameters are associated with the "B" mode.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat%ele(0) -- lat_struct: Twiss parameters at the start
!   ix_branch  -- Integer, optional: Branch index. Default is 0 (main lattice).
!   bmad_status  -- Common block status structure:
!       %type_out      -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   lat         -- lat_struct: Lat
!   bmad_status  -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

subroutine twiss_propagate_all (lat, ix_branch)

use bmad_struct
use bmad_interface, except_dummy => twiss_propagate_all

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer n, n_use
integer, optional :: ix_branch

! Propagate twiss

branch => lat%branch(integer_option(0, ix_branch))
n_use = branch%n_ele_track

bmad_status%ok = .true.

do n = 1, n_use
  call twiss_propagate1 (branch%ele(n-1), branch%ele(n))
  if (.not. bmad_status%ok) return
enddo

! Make sure final mode is same as initial mode

if (branch%param%lattice_type == circular_lattice$) then
  if (branch%ele(0)%mode_flip .neqv. branch%ele(n_use)%mode_flip) then
    call do_mode_flip (branch%ele(n_use))
    if (bmad_status%type_out .and. .not. bmad_status%ok) then
      print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
      print *, '      EQUAL TO THE INITIAL'
    endif
  endif
endif

end subroutine
