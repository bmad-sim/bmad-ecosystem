!+
! Subroutine track_runge_kutta (start, end, s_start, s_end, eps, 
!                                del_s_step, del_s_min, func_type, param)
!
! Subroutine to do tracking using Runge-Kutta integration. 
! The core Runge-Kutta routine used here is odeint from Numerical Recipes.
! See the Numerical Recipes in F90 book for more information.
!
! Note: To use track_runge_kutta you need to link in a subroutine called:
!   field_rk (position, field)
! Arguments are:
!   Input:  position(3)  -- Real: (x, y, s) position.
!   Output: field(3)     -- Real: See func_type for more info.
! 
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!   use ode_path    ! This only needed if you want to see intermediate values
!
! Input:
!   start      -- Coord_struct: Starting coords.
!
! Output:
!   end        -- Coord_struct: Ending coords.
!   s_start    -- Real: Starting point.
!   s_end      -- Real: Ending point.
!   eps        -- Real: Error tollerance. A good value would be, say, 1e-5.
!   del_s_step -- Real: Initial step size.
!   del_s_min  -- Real: Minimum step size.
!   func_type  -- Character*(*): Descripter for field_rk function:
!                   = 'B_FIELD'     -> field_rk returns the B field (Tesla).
!                   = 'KICK_FIELD'  -> field_rk returns the kicks.
!   param      -- Param_struct: [Optional] Beam parameters.
!                 [Needed if func_type = 'B_FIELD']
!     %enegy     -- Energy in GeV
!     %particle  -- Particle type [positron$, or electron$]
!
!
! Ode_path module Input:
!   save_steps -- Logical: Set True if you want the intermediate results saved.
!
! Ode_path module Output:
!   kout    -- Integer: Number of data points.
!   xp(:)   -- Real: Array of s locations.
!   yp(:,:) -- Real: y(:,i) holds the coords at s = xp(i).
!- 

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:55  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine track_runge_kutta (start, end, s_start, s_end, eps, &
                                  del_s_step, del_s_min, func_type, param)


  use bmad_struct
  use bmad_interface
  use nr, only: odeint

  implicit none

  type (coord_struct) start, end
  type (param_struct), optional :: param

  real s_start, s_end, eps, del_s_step, del_s_min

  character*(*) func_type

  external derivs

! init

  bmad_common%func_type = func_type

  if (func_type == 'B_FIELD') then
    bmad_common%factor = 0.2997 * param%particle / &
                                           param%energy (1 + start%z%vel)
  elseif (func_type == 'KICK_FIELD') then
    bmad_common%factor = 1.0
  else
    print *, 'ERROR IN TRACK_RUNGE_KUTTA: UNKNOWN "FUNC_TYPE": ', func_type
  endif

!

  end = start
  end%x%vel = end%x%vel / (1 + end%z%vel)
  end%y%vel = end%y%vel / (1 + end%z%vel)

  call odeint (end%vec, x_start, s_end, eps, del_s_step, del_s_min, &
                                                            derivs, rkqs)

  end%x%vel = end%x%vel * (1 + end%z%vel)
  end%y%vel = end%y%vel * (1 + end%z%vel)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine derivs (x, y, dydx)

  use bmad_struct
  use nrtype

  implicit none
                                   
  real(sp), intent(in) :: x
  real(sp), intent(in) :: y(:)
  real(sp), intent(out) :: dydx(;)

  type (coord_struct) here

  real field(3)

!

  
  call field_rk ((/ y(1), y(3), x /), field)

  dydx(1) = y(2)    ! dx/ds = 
  dydx(2) = 
  dydx(3) = y(4)
  dydx(4) = 
  dydx(5) = 1 - sqrt(1 + y(2)**2 + y(4)**2)
  dydx(6) = 0       ! dE/ds = 0

end subroutine


