!+
! Subroutine twiss_propagate_all (lat)
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
!   bmad_status  -- Common block status structure:
!       %type_out      -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   lat         -- lat_struct: Lat
!   bmad_status  -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

#include "CESR_platform.inc"

subroutine twiss_propagate_all (lat)

use bmad_struct
use bmad_interface, except_dummy => twiss_propagate_all

implicit none

type (lat_struct) :: lat

integer n, n_use

! Propagate twiss

n_use = lat%n_ele_track

bmad_status%ok = .true.

do n = 1, n_use
  call twiss_propagate1 (lat%ele(n-1), lat%ele(n))
  if (.not. bmad_status%ok) return
enddo

! Make sure final mode is same as initial mode

if (lat%param%lattice_type == circular_lattice$) then
  if (lat%ele(0)%mode_flip .neqv. lat%ele(n_use)%mode_flip) then
    call do_mode_flip (lat%ele(n_use))
    if (bmad_status%type_out .and. .not. bmad_status%ok) then
      print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
      print *, '      EQUAL TO THE INITIAL'
    endif
  endif
endif

end subroutine
