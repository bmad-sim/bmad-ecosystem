!+
! Subroutine twiss_and_track_at_s (ring, s, ele, orb_, here)
! 
! Subroutine to return the twiss parameters and particle orbit at a 
! given longitudinal position. See also twiss_and_track_partial
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Ring holding the lattice.
!   s    -- Real(rp): Longitudinal position. If s is negative the
!            the position is taken to be ring%param%total_length - s.
!   orb_(0:n_ele_max) -- Coord_struct: Orbit through the ring.
!                           (usually the closed orbit)
!
! Output:
!   ele  -- Ele_struct: Element structure holding the Twiss parameters.
!             if orb_ is not given then the Twiss parameters are calculated
!             with respect to the zero orbit.
!   here -- Coord_struct, optional: Particle position at the position s.
!-

#include "CESR_platform.inc"

subroutine twiss_and_track_at_s (ring, s, ele, orb_, here)

  use bmad_struct
  use bmad_interface, except => twiss_and_track_at_s

  implicit none

  type (ring_struct) :: ring
  type (ele_struct) :: ele
  type (coord_struct), optional :: orb_(0:)
  type (coord_struct), optional :: here

  real(rp) s, s_use

  integer i

! For negative s use ring%param%total_length - s
! Actually take into account that the ring may start from some non-zero s.

  s_use = s
  if (s < ring%ele_(0)%s) s_use = ring%param%total_length - s

! error_check

  i = ring%n_ele_use
  if (s_use < ring%ele_(0)%s .or. s > ring%ele_(i)%s) then
    print *, 'ERROR IN TWISS_AND_TRACK_AT_S: S POSITION OUT OF BOUNDS.', s
    call err_exit
  endif

! Propagate to position

  do i = 1, ring%n_ele_use
    if (abs(ring%ele_(i)%s - s_use) < 1e-5) then
      ele = ring%ele_(i)
      here = orb_(i)
      return
    elseif (ring%ele_(i)%s > s_use) then
      if (present(orb_)) then
        call twiss_and_track_partial (ring%ele_(i-1), ring%ele_(i), &
                   ring%param, s_use-ring%ele_(i-1)%s, ele, orb_(i-1), here)
      else
        call twiss_and_track_partial (ring%ele_(i-1), ring%ele_(i), &
                   ring%param, s_use-ring%ele_(i-1)%s, ele)
      endif
      return
    endif
  enddo

end subroutine
