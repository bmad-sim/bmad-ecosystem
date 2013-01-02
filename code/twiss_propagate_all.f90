!+
! Subroutine twiss_propagate_all (lat, ix_branch, err_flag)
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
!   err_flag     -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine twiss_propagate_all (lat, ix_branch, err_flag)

use bmad_struct
use bmad_interface, except_dummy => twiss_propagate_all

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: lord, slave, ele

integer n, n_track
integer, optional :: ix_branch

logical, optional :: err_flag
logical err

! Make sure gamma for ele(0) is correct.

branch => lat%branch(integer_option(0, ix_branch))
ele => branch%ele(0)

if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta

! Propagate twiss

n_track = branch%n_ele_track

if (present(err_flag)) err_flag = .true.

do n = 1, n_track
  call twiss_propagate1 (branch%ele(n-1), branch%ele(n), err)
  if (err) return
enddo

! Make sure final mode is same as initial mode

if (branch%param%lattice_type == circular_lattice$) then
  if (branch%ele(0)%mode_flip .neqv. branch%ele(n_track)%mode_flip) then
    call do_mode_flip (branch%ele(n_track), err)
    if (err .and. global_com%type_out) then
      print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
      print *, '      EQUAL TO THE INITIAL'
    endif
  endif
endif

! Super_lord elements get the twiss parameters at the exit end

do n = lat%n_ele_track + 1, lat%n_ele_max
  lord => lat%ele(n)
  select case (lord%lord_status) 
  case (super_lord$, overlay_lord$, group_lord$)
    if (lord%n_slave == 0) cycle
    slave => pointer_to_slave(lord, lord%n_slave)
    if (slave%ix_branch /= branch%ix_branch) cycle
    call transfer_twiss (slave, lord)
  end select
enddo

if (present(err_flag)) err_flag = .false.

end subroutine
