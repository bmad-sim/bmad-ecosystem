!+
! Subroutine twiss_at_s (ring, s, ele)
! 
!------------------------------------------------------------------
! THIS ROUTINE IS OBSOLETE. PLEASE USE TWISS_AND_TRACK_AT_S INSTEAD.
!------------------------------------------------------------------
!
! Subroutine to return the twiss parameters at a specific longitudinal 
! position.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Ring holding the lattice.
!   s    -- Real(rdef): Longitudinal position. If s is negative the
!            the position is taken to be ring%param%total_length - s.
!
! Output:
!   ele -- Ele_struct: Element structure holding the twiss_parameters.
!-

#include "CESR_platform.inc"

subroutine twiss_at_s (ring, s, ele)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) :: ring
  type (ele_struct) :: ele

  real(rp) s, s_use

  integer i

! For negative s use ring%param%total_length - s
! Actually take into account that the ring may start from some non-zero s.

  s_use = s
  if (s < ring%ele_(0)%s) s_use = ring%param%total_length - s

! error_check

  i = ring%n_ele_ring
  if (s_use < ring%ele_(0)%s .or. s > ring%ele_(i)%s) then
    print *, 'ERROR IN TWISS_AT_S: S POSITION OUT OF BOUNDS.', s
    call err_exit
  endif

! Propagate to position

  do i = 1, ring%n_ele_ring
    if (abs(ring%ele_(i)%s - s_use) < 1e-5) then
      ele = ring%ele_(i)
      return
    elseif (ring%ele_(i)%s > s_use) then
      call twiss_and_track_partial (ring%ele_(i-1), ring%ele_(i), &
                                     ring%param, s_use-ring%ele_(i-1)%s, ele)
      return
    endif
  enddo

end subroutine
