#include "CESR_platform.inc"

module dynamic_aperture_mod

  use bmad_struct
  use bmad_interface

  type aperture_struct
    type (coord_struct)  closed_orbit
    real(rp) x, y
    integer plane
    integer ix_ring
    integer i_turn
  end type

  type track_input_struct
    integer n_turn
    real(rp) x_init, y_init
    real(rp) accuracy
  end type


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

contains

!+
! Subroutine dynamic_aperture (ring, orb0, theta_xy, track_input, aperture)
!
! Subroutine to determine the dynamic aperture of a lattice by tracking.
! The subroutine works by determining where on a radial line y = const * x
! the aperture is. Here x and y are deviations from the closed orbit.
!
! Modules Needed:
!   use dynamic_aperture_mod
!
! Input:
!   ring        -- Ring_struct: Ring containing the lattice.
!   orb0        -- Coord_struct: Closed orbit at the start.
!     %vec(6)      -- Energy offset.
!   theta_xy    -- Real(rp): Angle of radial line (in radians) in x-y space.
!                    Angle is "normalized" by %x_init, %y_init.
!   track_input -- Track_input_struct: Structure holding the input data:
!     %n_turn     -- Number of turns tracked.
!     %x_init     -- Initial x coordinate to start with for theta_xy = 0.
!     %y_init     -- Initial y coordinate to start with for theta_xy = pi/2.
!     %accuracy   -- Accuracy needed of aperture results.
!
! Output:
!     aperture  -- Aperture_struct:
!       %closed_orbit -- Closed orbit coordinates
!       %x            -- X at aperture limit
!       %y            -- Y at aperture limit
!       %plane        -- Plane in which lost (X_PLANE$ or Y_PLANE$)
!       %ix_ring      -- Index where lost
!       %i_turn       -- Turn where lost
!
! Note: The radial lines are spaced equally in angle using coordinates
!       normalized by %X_INIT and %Y_INIT
!-

subroutine dynamic_aperture (ring, orb0, theta_xy, track_input, aperture)

  implicit none

  type (ring_struct)  ring
  type (param_struct)  param_save
  type (coord_struct)  orb0
  type (coord_struct), save, allocatable :: orbit_(:)
  type (aperture_struct)  aperture
  type (track_input_struct)  track_input

  integer it, i, turn_lost, ixr, ie_max

  real(rp) eps_rel(4), eps_abs(4)
  real(rp) e_init, theta_xy                                   
  real(rp) x0, x1, x2, y0, y1, y2

  logical aperture_bracketed

! init setup

  if (track_input%x_init == 0 .or. track_input%y_init == 0) then
    print *, 'ERROR IN DYNAMIC_APERTURE: TRACK_INPUT.X_INIT OR',  &
                                             ' TRACK_INPUT.Y_INIT = 0'
    call err_exit
  endif

  param_save = ring%param
  ring%param%aperture_limit_on = .true.

  call reallocate_coord (orbit_, ring%n_ele_max)

! Find starting point

  x0 = 0
  y0 = 0
  x1 = cos(theta_xy) * track_input%x_init
  y1 = sin(theta_xy) * track_input%y_init

  aperture_bracketed = .false.

! use a binary search to find where the aparture is along the line

  test_loop: do

    orbit_(0) = orb0
    orbit_(0)%vec(1) = orbit_(0)%vec(1) + x1
    orbit_(0)%vec(3) = orbit_(0)%vec(3) + y1

! track n_turns

    do it = 1, track_input%n_turn
      call track_all (ring, orbit_)
      if (ring%param%lost) exit
      orbit_(0) = orbit_(ring%n_ele_use)
    enddo

! change search interval end

    if (ring%param%lost) then
      x2 = x1
      y2 = y1
      turn_lost = it
      aperture_bracketed = .true.
    else
      x0 = x1
      y0 = y1
    endif

! calculate new starting point

    if (aperture_bracketed) then
      x1 = (x0 + x2) / 2
      y1 = (y0 + y2) / 2
    else
      x1 = 2 * x1
      y1 = 2 * y1
    endif

! check to see if there is an aperture

    if (x1 > 1000*track_input%x_init .or.  &
                              y1 > 1000*track_input%y_init) then
      print *, 'ERROR IN DYNAMIC_APERTURE: CANNOT FIND APERTURE LIMIT'
      call err_exit
    endif

! see if we are accurate enough

    if (aperture_bracketed) then
      if (sqrt((x2-x0)**2 + (y2-y0)**2) <=  &
                                   2*track_input%accuracy) exit test_loop
    endif

  enddo test_loop

! fill in the info

  aperture%x = x1
  aperture%y = y1
  aperture%i_turn = turn_lost
  aperture%closed_orbit = orb0

  ixr = ring%param%ix_lost
  if (ring%ele_(ixr)%value(x_limit$) /= 0 .and.  &
                abs(orbit_(ixr)%vec(1)) > ring%ele_(ixr)%value(x_limit$)) then
    aperture%plane = x_plane$
  else
    aperture%plane = y_plane$
  endif

  if (ring%ele_(ixr)%key == drift$) ixr = ixr + 1
  aperture%ix_ring = ixr

!

  ring%param = param_save

end subroutine

end module
