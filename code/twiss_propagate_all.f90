!+
! Subroutine TWISS_PROPAGATE_ALL (RING)
!
! Subroutine to propagate the twiss parameters from the start to the end.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING.ELE_(0) -- Ring_struct: Twiss parameters at the start
!     STATUS       -- Common block status structure
!       .TYPE_OUT     -- Logical: If .true. then will type a message if the
!                        modes are flipped.
!
! Output:
!     RING     -- Ring_struct: Ring
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine twiss_propagate_all (ring)

  use bmad_struct
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
  if (ring%param%symmetry == mobius_symmetry$) bmad_status%type_out = .false.

  do n = 1, n_use
    call twiss_propagate1 (ring%ele_(n-1), ring%ele_(n))
  enddo

  bmad_status%type_out = temp_type_out

! make sure final mode is same as initial mode

  if (n_use == ring%n_ele_ring) then
    if (ring%ele_(0)%mode_flip .xor. ring%ele_(n_use)%mode_flip) then
      call do_mode_flip (ring%ele_(n_use), ring%ele_(n_use))
      if (bmad_status%type_out .and. .not. bmad_status%ok) then
        type *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
        type *, '      EQUAL TO THE INITIAL'
      endif
    endif
  endif

  return
  end
