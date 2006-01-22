!+
! Subroutine twiss_and_track_at_s (ring, s, ele, orb, orb_at_s)
! 
! Subroutine to return the twiss parameters and particle orbit at a 
! given longitudinal position. See also twiss_and_track_partial.
!
! Note: When calculating the Twiss parameters, this routine assumes 
! that the lattice elements already contain the Twiss parameters calculated
! for the ends of the elements. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Ring holding the lattice.
!   s    -- Real(rp): Longitudinal position. If s is negative the
!            the position is taken to be ring%param%total_length - s.
!   orb(0:n_ele_max) -- Coord_struct, optional: Orbit through the ring.
!                           (usually the closed orbit)
!
! Output:
!   ele  -- Ele_struct, optional: Element structure holding the Twiss parameters.
!             if orb is not given then the Twiss parameters are calculated
!             with respect to the zero orbit.
!   orb_at_s -- Coord_struct, optional: Particle position at the position s.
!             If orb_at_s is present then this routine assumes that orb is
!             present.
!-

#include "CESR_platform.inc"

subroutine twiss_and_track_at_s (ring, s, ele, orb, orb_at_s)

  use bmad_struct
  use bmad_interface, except => twiss_and_track_at_s

  implicit none

  type (ring_struct) :: ring
  type (ele_struct), optional :: ele
  type (coord_struct), optional :: orb(0:)
  type (coord_struct), optional :: orb_at_s

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

    ! Test if we have the correct element. The factor of 1e-5 is for roundoff.

    if (s_use - ring%ele_(i)%s < 1e-5) cycle

    ! If close enough to edge of element just use element info.

    if (s_use - ring%ele_(i)%s < -1e-5) then
      if (present(ele)) ele = ring%ele_(i)
      if (present(orb_at_s)) orb_at_s = orb(i)
      return
    endif

    ! Normal case where we need to partially track through

    if (present(orb)) then
      call twiss_and_track_partial (ring%ele_(i-1), ring%ele_(i), &
                   ring%param, s_use-ring%ele_(i-1)%s, ele, orb(i-1), orb_at_s)
    else
      call twiss_and_track_partial (ring%ele_(i-1), ring%ele_(i), &
                   ring%param, s_use-ring%ele_(i-1)%s, ele)
    endif

    return

  enddo

end subroutine
