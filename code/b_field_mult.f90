!+
! Subroutine B_FIELD_MULT (RING, COORD, FIRST, LAST, S_POS, B_VECTOR)
!
!   Subroutine to calculate the magnetic field vector due to multiple circular
! current loops.
! -- Created by Daniel Fromowitz, September 1999.
!
! Input:
!     RING  -- Ring_struct: Ring
!     COORD -- Coord_struct: TRUE cartesian coordinates of particle, i.e. z is
!                            relative to the (linac) origin; it is not a
!                            displacement!
!     FIRST -- Integer: Index of first element in section
!     LAST  -- Integer: Index of last element in section
!     S_POS -- Real: Array of longitudinal positions of coil components
!
! Output:
!     B_VECTOR(3) -- Real: (Cartesian) Magnetic field vector x, y, and z
!                          components (in units of mu_0 / 2)
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine b_field_mult (ring, coord, first, last, s_pos, b_vector)
  use bmad_struct
  implicit none

  type (ring_struct)  ring
  type (coord_struct)  coord
  integer first, last
  real s_pos(n_comp_maxx), b_vector(3)

  real b_loop(3)
  integer i, x$, y$, z$
  parameter (x$ = 1, y$ = 2, z$ = 3)

  do i = 1, 3
    b_vector(i) = 0.0
  enddo

  do i = first, last
    if (ring%ele_(i)%key == loop$) then
      call b_field_loop (coord, ring%ele_(i), s_pos(i - first + 1), b_loop)
      b_vector(x$) = b_vector(x$) + b_loop(x$)
      b_vector(y$) = b_vector(y$) + b_loop(y$)
      b_vector(z$) = b_vector(z$) + b_loop(z$)
    elseif (ring%ele_(i)%key /= drift$) then
      type *, 'ERROR IN B_FIELD_MULT.F:'
      type *, 'COIL CONTAINER DOES NOT KNOW HOW TO USE A ',  &
        key_name(ring%ele_(i)%key),'.'
      type *, 'EXITING.'
      call err_exit
    endif
  enddo

  return
  end
