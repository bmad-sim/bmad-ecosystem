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

real(rp) rel_tol, abs_tol, dref_time, beta0, s0, s1, ds_ref

logical err_flag

! Convert to element coords.
! For a patch, convert to the downstream coords so that the downstream face 
! can be simply described as being at s = s1. Additionally, with a patch, s 
! is the distance before the downstream face so the s0 starting position is 
! negative and the s1 stopping position is 0.

start2_orb = start_orb
beta0 = ele%value(p0c$) / ele%value(e_tot$)

if (ele%key == patch$) then
  call track_a_patch (ele, start2_orb, .false., s0, ds_ref)
  s0 = s0 * start2_orb%direction * ele%orientation
  s1 = 0
  start2_orb%vec(5) = start2_orb%vec(5) + (ds_ref + s0) * start2_orb%beta / beta0 
else
  call offset_particle (ele, param, set$, start2_orb, set_hvkicks = .false., set_multipoles = .false.)
  s0 = 0; s1 = ele%value(l$)
endif

! Track.

call odeint_bmad (start2_orb, ele, param, end_orb, s0, s1, .true., err_flag, track)
if (err_flag) return

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

! convert to lab coords.

if (ele%key /= patch$) then
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false.)
endif

! The z value computed in odeint_bmad is off for elements where the particle changes energy is not 
! constant (see odeint_bmad for more details). In this case make the needed correction.
! dref_time is reference time for transversing the element under the assumption, used by 
! odeint_bmad, that the reference velocity is constant and equal to the velocity at the final enegy.

if (ele%key /= patch$) then
  dref_time = ele%value(l$) / (beta0 * c_light)
  end_orb%vec(5) = end_orb%vec(5) + (ele%value(delta_ref_time$) - dref_time) * end_orb%beta * c_light
endif

end subroutine
