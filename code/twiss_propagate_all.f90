!+
! Subroutine twiss_propagate_all (ring)
!
! Subroutine to propagate the twiss parameters from the start to the end.
!
! Note: ring%ele_()%x Twiss parameters are associated with the "A" mode and
! the ring%ele_()%y Twiss parameters are associated with the "B" mode.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring.ELE_(0) -- Ring_struct: Twiss parameters at the start
!   bmad_status -- Common block status structure:
!       %type_out -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   ring         -- Ring_struct: Ring
!   bmad_status -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

#include "CESR_platform.inc"

subroutine twiss_propagate_all (ring)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring

  integer n, n_use
  logical temp_type_out

! init phase

  ring%ele_(0)%x%phi = 0
  ring%ele_(0)%y%phi = 0

! propagate twiss

  n_use = ring%n_ele_use

  temp_type_out = bmad_status%type_out
  bmad_status%ok = .true.

  do n = 1, n_use
    call twiss_propagate1 (ring%ele_(n-1), ring%ele_(n))
    if (.not. bmad_status%ok) return
  enddo

  bmad_status%type_out = temp_type_out

! make sure final mode is same as initial mode

  if (n_use == ring%n_ele_ring) then
    if (ring%ele_(0)%mode_flip .neqv. ring%ele_(n_use)%mode_flip) then
      call do_mode_flip (ring%ele_(n_use), ring%ele_(n_use))
      if (bmad_status%type_out .and. .not. bmad_status%ok) then
        print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
        print *, '      EQUAL TO THE INITIAL'
      endif
    endif
  endif

  return
  end
