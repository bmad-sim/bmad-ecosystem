!+
! Subroutine twiss_at_s (ring, s, ele)
! 
! Subroutine to return the twiss parameters at a specific longitudinal 
! position.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring -- Ring_struct: Ring holding the lattice.
!   s    -- Real: Longitudinal position.
!
! Output:
!   ele -- Ele_struct: Element structure holding the twiss_parameters.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine twiss_at_s (ring, s, ele)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) :: ring
  type (ele_struct) :: ele

  real s

  integer i

! error_check

  if (s < 0 .or. s > ring%param%total_length) then
    type *, 'ERROR IN TWISS_AT_S: S POSITION OUT OF BOUNDS.', s
    call err_exit
  endif

!

  do i = 1, ring%n_ele_ring
    if (abs(ring%ele_(i)%s - s) < 1e-5) then
      ele = ring%ele_(i)
      return
    elseif (ring%ele_(i)%s > s) then
      call twiss_and_track_partial (ring%ele_(i-1), ring%ele_(i), &
                                         ring%param, s-ring%ele_(i-1)%s, ele)
      return
    endif
  enddo

end subroutine
