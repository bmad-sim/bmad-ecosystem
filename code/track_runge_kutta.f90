!+
! Subroutine track_runge_kutta (start, end, s_start, s_end, rel_eps, abs_eps,
!                                      del_s_step, del_s_min, func_type, param)
!
! Subroutine to do tracking using Runge-Kutta integration. 
! The core Runge-Kutta routine used here is odeint2 which is
! a modified version of odeint from Numerical Recipes.
! See the "Numerical Recipes in F90" book and the odeint.f90 file
! for more information.
!
! Note: To use track_runge_kutta you need to link in a subroutine called:
!       field_rk (position, field)
! Arguments are:
!   Input:  position  -- Coord_struct: Position of particle. 
!                          Here position%z%pos = s coordinate.
!   Output: field(3)  -- Real(rdef): See func_type for more info.
! 
! See the html bmad programming notes for more details on how to use this
! subroutine.
!
! Modules needed:
!   use bmad
!   use ode_path    ! This only needed if you want to see intermediate values
!
! Input:
!   start      -- Coord_struct: Starting coords.
!   s_start    -- Real(rdef): Starting point.
!   s_end      -- Real(rdef): Ending point.
!   rel_eps    -- Real(rdef): % Error tollerance. A good value would be, say, 1e-5.
!   abs_eps    -- Real(rdef): absolute error. A good value would be, say, 1e-8.
!   del_s_step -- Real(rdef): Initial guess for a step size.
!   del_s_min  -- Real(rdef): Minimum step size (can be zero).
!   func_type  -- Character*(*): Descripter for field_rk function:
!                   = 'B_FIELD'     -> field_rk returns the B field (Tesla).
!                   = 'KICK_FIELD'  -> field_rk returns the kicks.
!   param      -- Param_struct: [Optional] Beam parameters.
!                 [Needed if func_type = 'B_FIELD']
!     %enegy     -- Energy in GeV
!     %particle  -- Particle type [positron$, or electron$]
!
! Output:
!   end        -- Coord_struct: Ending coords.
!     %z%pos     = s_len - path_length where s_len = s_end-s_start. 
!                  Thus for a wiggler, where the "zero" orbit path_length 
!                  is not s_len, there needs to be a correction term:
!                    end%z%pos = end%z%pos + zero_orbit_path_length - s_len
!
!
! Ode_path module Input:
!   save_steps -- Logical: Set True if you want the intermediate results saved.
!
! Ode_path module Output:
!   kount   -- Integer: Number of data points.
!   xp(:)   -- Real(rdef): Array of s locations.
!   yp(:,:) -- Real(rdef): y(:,i) holds the coords at s = xp(i).
!- 

!$Id$
!$Log$
!Revision 1.6  2002/02/23 20:32:26  dcs
!Double/Single Real toggle added
!
!Revision 1.5  2002/02/01 16:03:07  dcs
!*** empty log message ***
!
!Revision 1.4  2002/02/01 15:52:34  dcs
!bmad_common%factor set for 'KICK_FIELD'
!
!Revision 1.3  2002/01/29 16:44:44  dcs
!Fix comments
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine track_runge_kutta (start, end, s_start, s_end, rel_eps, abs_eps, &
                                       del_s_step, del_s_min, func_type, param)


  use bmad
  use nr

  implicit none

  type (coord_struct) start, end
  type (param_struct), optional :: param

  real(rdef) s_start, s_end, rel_eps, abs_eps, del_s_step, del_s_min, fac

  character*(*) func_type

  external derivs

! init

  bmad_common%func_type = func_type
  bmad_common%factor = 0.2997 * param%particle / &
                                       (param%energy * (1 + start%z%vel))

  if (func_type /= 'B_FIELD' .and. func_type/= 'KICK_FIELD') then
    print *, 'ERROR IN TRACK_RUNGE_KUTTA: UNKNOWN "FUNC_TYPE": ', func_type
    call err_exit
  endif

!

  end = start

! convert from cononical momentum P_x, P_y to x', y'

  fac = sqrt((1 + end%z%vel)**2 - end%x%vel**2 - end%y%vel**2)
  end%x%vel = end%x%vel / fac
  end%y%vel = end%y%vel / fac

! track

  call odeint2 (end%vec, s_start, s_end, rel_eps, abs_eps, &
                                  del_s_step, del_s_min, derivs, rkqs)

! convert back to cononical momentum

  fac = sqrt((1 + end%z%vel)**2 - end%x%vel**2 - end%y%vel**2)
  end%x%vel = end%x%vel * fac
  end%y%vel = end%y%vel * fac

end subroutine
