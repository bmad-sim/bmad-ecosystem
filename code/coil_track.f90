!+
! Subroutine COIL_TRACK (START, ELE_INDEX, RING, END)
!
!   Subroutine to track a particle through current loops via iteration of
! difference equations.  (Some extremely low-energy particles---below about
! 1.8 MeV---may track incorrectly.)
! -- Created by Daniel Fromowitz, January 2000.
!
! Input:
!     START     -- Coord_struct: Coordinates before iteration
!     ELE_INDEX -- Integer: Index of coil
!     RING      -- Ring_struct: Ring
!
! Output:
!     END       -- Coord_struct: Coordinates after iteration
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine coil_track (start, ele_index, ring, end)

  use bmad_struct

  implicit none

  type (ring_struct)  ring
  type (coord_struct)  start, end
  integer ele_first, ele_last
  integer ele_index


  real dt, dt_partial, dt_temp, g, g_partial, g_temp, time_disp, time_total
  real z_end, z_start, b_vector(3), mu_0_over_2, original_z_prime
  real x_dot_temp, y_dot_temp, s_pos(n_comp_maxx), s_cumul
  real x_lim, y_lim
  parameter (mu_0_over_2 = twopi * 1.0e-7)
  integer i, ic, x$, y$, z$
  parameter (x$ = 1, y$ = 2, z$ = 3)

! check

  if (ring%param%lattice_type /= linac_lattice$) then
    type *, 'ERROR IN COIL_TRACK:'
    type *, 'COIL ELEMENT REQUIRES LATTICE_TYPE := LINAC_LATTICE.'
    type *, 'EXITING.'
    call err_exit
  endif

!

  end = start     ! transfer start to end

  ic = ring%ele_(ele_index)%ic1_lord
  ele_first = ring%control_(ic)%ix_lord
  ic = ring%ele_(ele_index)%ic2_lord
  ele_last  = ring%control_(ic)%ix_lord
  if (ele_last - ele_first + 1  > n_comp_maxx) then
    type *, 'ERROR IN COIL_TRACK:'
    type '(a,a,i4)', ' THERE ARE TOO MANY COIL COMPONENTS.  RESET',  &
      ' ''n_comp_maxx'' TO AT LEAST ',ele_last - ele_first + 1
    type *,'IN BMAD_STRUCT.'
    type *, 'EXITING.'
    call err_exit
  endif

  original_z_prime = end%z%vel
  z_start = ring%ele_(ele_index - 1)%s

!     Create array of actual longitudinal positions of coil components
  s_cumul = z_start
  do i = ele_first, ele_last
    s_cumul = s_cumul + ring%ele_(i)%value(l$)
    s_pos(i - ele_first + 1) = s_cumul
  enddo

  call change_basis (end, ring%param%energy, z_start, .true., time_disp)

  z_end = ring%ele_(ele_index)%s
  time_total = 0.0
  dt = 5.0e-13
  g = ring%param%particle * mu_0_over_2 * c_light**2 * dt /  &
      (ring%param%energy * 1.0e9 * (start%z%vel + 1))
  dt_partial = dt * c_light / 20.
  g_partial = ring%param%particle * mu_0_over_2 * c_light**2 /  &
              (ring%param%energy * 1.0e9 * (start%z%vel + 1))

  do while (end%z%pos < z_end)

    call b_field_mult (ring, end, ele_first, ele_last, s_pos, b_vector)
    x_dot_temp = g * (end%y%vel * b_vector(z$) - end%z%vel * b_vector(y$))  &
                 + end%x%vel
    y_dot_temp = g * (end%z%vel * b_vector(x$) - end%x%vel * b_vector(z$))  &
                 + end%y%vel
    end%z%vel  = g * (end%x%vel * b_vector(y$) - end%y%vel * b_vector(x$))  &
                 + end%z%vel
    end%x%vel = x_dot_temp
    end%y%vel = y_dot_temp
    end%x%pos = end%x%pos + x_dot_temp * dt
    end%y%pos = end%y%pos + y_dot_temp * dt
    time_total = time_total + dt
    if (end%z%vel >= c_light / 20.) then
      end%z%pos = end%z%pos + end%z%vel * dt
    elseif (end%z%vel > 0.0) then
!         Use a larger time step to speed the particle tracking
      end%z%pos = end%z%pos + end%z%vel * dt
      dt_temp = dt_partial / end%z%vel
      g_temp = g_partial * dt_temp
      call b_field_mult (ring, end, ele_first, ele_last, s_pos, b_vector)
      x_dot_temp = (end%y%vel * b_vector(z$) - end%z%vel * b_vector(y$))  &
                   * g_temp + end%x%vel
      y_dot_temp = (end%z%vel * b_vector(x$) - end%x%vel * b_vector(z$))  &
                   * g_temp + end%y%vel
      end%z%vel  = (end%x%vel * b_vector(y$) - end%y%vel * b_vector(x$))  &
                   * g_temp + end%z%vel
      end%x%vel = x_dot_temp
      end%y%vel = y_dot_temp
      end%x%pos = end%x%pos + x_dot_temp * dt_temp
      end%y%pos = end%y%pos + y_dot_temp * dt_temp
      time_total = time_total + dt_temp
      end%z%pos = end%z%pos + end%z%vel * dt_temp
    else
!         Disregard particles reflected from a magnetic mirror
      end%z%pos = z_end
      ring%param%lost = .true.
    endif

  enddo

! Back up so particle is exactly at the end of the section

  dt = (z_end - end%z%pos) / end%z%vel
  time_total = time_total + dt
  g = ring%param%particle * mu_0_over_2 * c_light**2 * dt /  &
      (ring%param%energy * 1.0e9 * (start%z%vel + 1))
  end%x%pos = end%x%pos + end%x%vel * dt
  end%y%pos = end%y%pos + end%y%vel * dt
  end%z%pos = z_end
  call b_field_mult (ring, end, ele_first, ele_last, s_pos, b_vector)
  x_dot_temp = g * (end%y%vel * b_vector(z$) - end%z%vel * b_vector(y$))  &
               + end%x%vel
  y_dot_temp = g * (end%z%vel * b_vector(x$) - end%x%vel * b_vector(z$))  &
               + end%y%vel
  end%z%vel  = g * (end%x%vel * b_vector(y$) - end%y%vel * b_vector(x$))  &
               + end%z%vel
  time_disp = time_disp + time_total - (z_end - z_start) / c_light

  call change_basis (end, ring%param%energy, z_end, .false., time_disp)
! Force particle energy to remain the same despite rounding errors
  end%z%vel = original_z_prime

! Check for particles outside aperture

  if (ring%param%aperture_limit_on) then
    x_lim = ring%ele_(ele_index)%value(x_limit$)
    y_lim = ring%ele_(ele_index)%value(y_limit$)
    if (x_lim > 0 .and. abs(end%x%pos) > x_lim)  &
      ring%param%lost = .true.
    if (y_lim > 0 .and. abs(end%y%pos) > y_lim)  &
      ring%param%lost = .true.
  endif

  return
  end
