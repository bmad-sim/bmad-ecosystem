!+
! Subroutine DYNAMIC_APERTURE (RING, TRACK_INPUT, APERTURE)
!
! Subroutine to determine the dynamic aperture of a lattice by tracking.
! The subroutine works by determining where on a radial line y = const * x
! the aperture is.  Here x and y are deviations from the closed orbit.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   RING    -- Ring_struct: Ring containing the lattice.
!   TRACK_INPUT   -- Track_input_struct: Structure holding the input data:
!     %N_TURN         -- Number of turns tracked.
!     %X_INIT         -- Suggested initial x coordinate to start with.
!     %Y_INIT         -- Suggested initial y coordinate to start with.
!     %N_XY_PTS       -- Number of radial lines in a quarter plane.
!     %E_MAX          -- Maximum energy deviation of initial coordinates.
!     %N_ENERGY_PTS   -- Number of off energy runs to do.
!     %ACCURACY       -- Accuracy needed of aperture results.
!
! Output:
!     APERTURE_(*) -- Aperture_struct: Array of Structures for the results.
!                     Each element of the array is the results for one
!                     particular energy.
!     %dE_E       -- dE/E
!     %CLOSED_ORBIT -- Closed orbit coordinates
!     %X_(I)      -- X at aperture limit
!     %Y_(I)      -- Y at aperture limit
!     %PLANE_(I)  -- Plane in which lost (X_PLANE$ or Y_PLANE$)
!     %IX_RING(I) -- Index where lost
!     %I_TURN(I)  -- Turn where lost
!
! Note: The radial lines are spaced equally in angle using coordinates
!       normalized by %X_INIT and %Y_INIT
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:51  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine dynamic_aperture (ring, track_input, aperture_)

  use bmad_struct                                      
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (param_struct)  ring_param_state
  type (coord_struct)  closed_(0:n_ele_maxx), orbit_(0:n_ele_maxx), co
  type (aperture_struct)  aperture_(*)
  type (track_input_struct)  track_input
  type (mat627_struct) mats627(n_ele_maxx)

  integer i_e, i_xy, it, i, turn_lost, ixr, i_e_max

  real eps_rel(4), eps_abs(4)
  real e_init, theta                                   
  real x0, x1, x2, y0, y1, y2

  logical aperture_bracketed, track_on

!

  if (track_input%x_init == 0 .or. track_input%y_init == 0) then
    type *, 'ERROR IN DYNAMIC_APERTURE: TRACK_INPUT.X_INIT OR',  &
                                             ' TRACK_INPUT.Y_INIT = 0'
    call err_exit
  endif

  ring_param_state = ring%param

  call twiss_at_start (ring)
  call closed_orbit_at_start (ring, co, 4, .true.)
  eps_rel(:) = 0.000001
  eps_abs(:) = 0.000001
  call closed_orbit_from_tracking (ring, closed_, 4, eps_rel, eps_abs, co)
  ring%param%aperture_limit_on = .true.

  call ring_make_mat627 (ring, -1, +1, mats627)

  i_e_max = max(1, track_input%n_energy_pts)
  do i_e = 1, i_e_max

    e_init = (i_e - 1) * track_input%e_max / max (1, i_e_max - 1)
    aperture_(i_e)%dE_E = e_init
    aperture_(i_e)%closed_orbit = closed_(0)

    do i_xy = 1, track_input%n_xy_pts

      theta = (i_xy - 1) * pi / max(1, track_input%n_xy_pts - 1)

      x0 = 0
      y0 = 0
      x1 = cos(theta) * track_input%x_init
      y1 = sin(theta) * track_input%y_init

      aperture_bracketed = .false.
      track_on = .true.

      do while (track_on)

        do i = 1, 6
          orbit_(0)%vec(i) = closed_(0)%vec(i)
        enddo
        orbit_(0)%x%pos = orbit_(0)%x%pos + x1
        orbit_(0)%y%pos = orbit_(0)%y%pos + y1
        orbit_(0)%z%vel = orbit_(0)%z%vel + e_init

! track n_turns

        do it = 1, track_input%n_turn
          call track_long (ring, orbit_, 0, +1, mats627)
          if (ring%param%lost) exit
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
          type *, 'ERROR IN DYNAMIC_APERTURE: CANNOT FIND APERTURE LIMIT'
          call err_exit
        endif

! see if we are accurate enough

        if (aperture_bracketed) then
          if (sqrt((x2-x0)**2 + (y2-y0)**2) <=  &
                                        2*track_input%accuracy) then

            aperture_(i_e)%x_(i_xy) = x1
            aperture_(i_e)%y_(i_xy) = y1
            aperture_(i_e)%i_turn(i_xy) = turn_lost
            ixr = ring%param%ix_lost

            if (ring%ele_(ixr)%value(x_limit$) /= 0 .and.  &
                orbit_(ixr)%x%pos > ring%ele_(ixr)%value(x_limit$)) then
              aperture_(i_e)%plane_(i_xy) = x_plane$
            else
              aperture_(i_e)%plane_(i_xy) = y_plane$
            endif

            if (ring%ele_(ixr)%key == drift$) ixr = ixr + 1
            aperture_(i_e)%ix_ring(i_xy) = ixr

            track_on = .false.

          endif
        endif

!

      enddo ! while (track_on)

    enddo ! i_xy

  enddo ! i_e

!

  ring%param = ring_param_state
  return
  end
