!+
! Subroutine track1_runge_kutta (start, ele, param, end)
!
! Subroutine to do tracking using Runge-Kutta integration. 
! The core Runge-Kutta routine used here is odeint_bmad which is
! a modified version of odeint from Numerical Recipes.
! See the "Numerical Recipes in F90" book.
!
! Modules needed:
!   use bmad
!
! Input:
!   start      -- Coord_struct: Starting coords.
!   ele        -- Ele_struct
!     %value(rel_tol$)    -- Real(rp): % Error tollerance. 
!                             A good value would be, say, 1e-5.
!     %value(abs_tol$)    -- Real(rp): absolute error. 
!                             A good value would be, say, 1e-8.
!   param      -- Param_struct: Beam parameters.
!     %enegy     -- Energy in GeV
!     %particle  -- Particle type [positron$, or electron$]
!
! Output:
!   end        -- Coord_struct: Ending coords.
!     %vec(5)     = s_len - path_length where s_len = s_end-s_start. 
!                  Thus for a wiggler, where the "zero" orbit path_length 
!                  is not s_len, there needs to be a correction term:
!                    end%vec(5) = end%vec(5) + zero_orbit_path_length - s_len
!
!   odeint_com -- common structure holding the path from the last tracking.
!     %n_pts     -- Integer: Number of data points.
!     %s(:)      -- Real(rp): Array of s locations.
!     %orb(:)    -- Coord_struct: Array of coordinates
!       %vec(1)     -- X position, etc.
!- 

#include "CESR_platform.inc"

subroutine track1_runge_kutta (start, ele, param, end)

  use runge_kutta_mod
  use bmad_interface

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (param_struct), target, intent(inout) :: param
  type (ele_struct), target, intent(inout) :: ele

  real(rp) s_start, s_end, rel_tol, abs_tol, del_s_step, del_s_min, fac

! init

  if (ele%value(rel_tol$) == 0) then
    rel_tol = 1e-6
  else
    rel_tol = ele%value(rel_tol$)
  endif
  
  if (ele%value(abs_tol$) == 0) then
    abs_tol = 1e-7
  else
    abs_tol = ele%value(abs_tol$)
  endif
  
! Track.
! Remember: offset_particle converts from cononical momentum P_x, P_y to x', y'

  del_s_step = 1e-3
  del_s_min = 1e-8
  
  call offset_particle (ele, param, end, set$)
  call odeint_bmad (start, ele, param, end, 0.0_rp, ele%value(l$), &
                                  rel_tol, abs_tol, del_s_step, del_s_min)
  call offset_particle (ele, param, end, unset$)

end subroutine
