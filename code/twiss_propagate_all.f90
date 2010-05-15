!+
! Subroutine twiss_propagate_all (lat, ix_branch)
!
! Subroutine to propagate the twiss, coupling, and dispersion parameters from 
! the start to the end of a lattice branch.
!
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat        -- Lat_struct: lattice.
!     %branch(ix_branch%rele(0) -- Branch beginning element with the starting parameters.
!   ix_branch  -- Integer, optional: Branch index. Default is 0 (main lattice).
!   bmad_status  -- Common block status structure:
!       %type_out      -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   lat          -- lat_struct: Lattice with parameters computed for the branch.
!   bmad_status  -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

subroutine twiss_propagate_all (lat, ix_branch)

use bmad_struct
use bmad_interface, except_dummy => twiss_propagate_all

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: lord, slave

integer n, n_track
integer, optional :: ix_branch

! Propagate twiss

branch => lat%branch(integer_option(0, ix_branch))
n_track = branch%n_ele_track

bmad_status%ok = .true.

do n = 1, n_track
  call twiss_propagate1 (branch%ele(n-1), branch%ele(n))
  if (.not. bmad_status%ok) return
enddo

! Make sure final mode is same as initial mode

if (branch%param%lattice_type == circular_lattice$) then
  if (branch%ele(0)%mode_flip .neqv. branch%ele(n_track)%mode_flip) then
    call do_mode_flip (branch%ele(n_track))
    if (bmad_status%type_out .and. .not. bmad_status%ok) then
      print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
      print *, '      EQUAL TO THE INITIAL'
    endif
  endif
endif

! Super_lord elements get the twiss parameters at the exit end

do n = lat%n_ele_track + 1, lat%n_ele_max
  lord => lat%ele(n)
  if (lord%lord_status /= super_lord$) cycle
  slave => pointer_to_slave (lat, lord, lord%n_slave)
  if (slave%ix_branch /= branch%ix_branch) cycle
  call transfer_twiss (slave, lord)
enddo

end subroutine
