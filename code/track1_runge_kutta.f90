!+
! Subroutine track1_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)
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
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct
!   param      -- lat_param_struct: Beam parameters.
!     %enegy     -- Energy in GeV
!     %particle  -- Particle type [positron$, or electron$]
!
!   bmad_com -- Bmad common block.
!     %rel_tol_adaptive_tracking -- Relative tolerance. Default is 1e-6.
!     %abs_tol_adaptive_tracking -- Absolute tolerance. Default is 1e-7.
!
! Output:
!   end_orb    -- Coord_struct: Ending coords.
!   err_flag   -- Logical: Set True if there is an error. False otherwise.
!   track      -- Track_struct, optional: Structure holding the track information.
!- 

subroutine track1_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)

use runge_kutta_mod, except_dummy => track1_runge_kutta

implicit none

type (coord_struct) :: start_orb, end_orb, start2_orb
type (lat_param_struct), target, intent(inout) :: param
type (ele_struct), target, intent(inout) :: ele
type (track_struct), optional :: track

real(rp) rel_tol, abs_tol, del_s_step, del_s_min, dref_time, beta0

logical err_flag

! Init 

del_s_step = 1e-3
del_s_min = 1e-8

start2_orb = start_orb

! Convert to element coords

call offset_particle (ele, param, start2_orb, set$, set_canonical = .false., set_hvkicks = .false., set_multipoles = .false.)

! Track.

call odeint_bmad (start2_orb, ele, param, end_orb, 0.0_rp, ele%value(l$), bmad_com%rel_tol_adaptive_tracking, &
                  bmad_com%abs_tol_adaptive_tracking, del_s_step, del_s_min, .true., err_flag, track)
if (err_flag) return

! The z value computed in odeint_bmad is off for elements where the particle changes energy is not 
! constant (see odeint_bmad for more details). In this case make the needed correction.
! dref_time is reference time for transversing the element under the assumption, used by odeint_bmad, that 
! the reference velocity is constant and equal to the velocity at the final enegy.

beta0 = ele%value(p0c$) / ele%value(e_tot$)
dref_time = ele%value(l$) / (beta0 * c_light)
end_orb%vec(5) = end_orb%vec(5) + (ele%value(delta_ref_time$) - dref_time) * end_orb%beta * c_light

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

! convert to lab coords.

call offset_particle (ele, param, end_orb, unset$, set_canonical = .false., set_hvkicks = .false., set_multipoles = .false.)

end subroutine
