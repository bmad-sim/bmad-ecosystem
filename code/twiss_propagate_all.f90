!+
! Subroutine twiss_propagate_all (ring)
!
! Subroutine to propagate the twiss parameters from the start to the end.
!
! Note: ring%ele_(:)%x Twiss parameters are associated with the "A" mode and
! the ring%ele_(:)%y Twiss parameters are associated with the "B" mode.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring%ele_(0) -- Ring_struct: Twiss parameters at the start
!   bmad_status  -- Common block status structure:
!       %type_out      -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   ring         -- Ring_struct: Ring
!   bmad_status  -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

#include "CESR_platform.inc"

subroutine twiss_propagate_all (ring)

  use bmad_struct
  use bmad_interface, except => twiss_propagate_all

  implicit none

  type (ring_struct)  ring

  integer n, n_use

! Propagate twiss

  n_use = ring%n_ele_use

  bmad_status%ok = .true.

  do n = 1, n_use
    call twiss_propagate1 (ring%ele_(n-1), ring%ele_(n))
    if (.not. bmad_status%ok) return
  enddo

! Make sure final mode is same as initial mode

  if (ring%param%lattice_type == circular_lattice$) then
    if (ring%ele_(0)%mode_flip .neqv. ring%ele_(n_use)%mode_flip) then
      call do_mode_flip (ring%ele_(n_use), ring%ele_(n_use))
      if (bmad_status%type_out .and. .not. bmad_status%ok) then
        print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
        print *, '      EQUAL TO THE INITIAL'
      endif
    endif
  endif

end subroutine
