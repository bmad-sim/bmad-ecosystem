!+
! Subroutine track1_runge_kutta (start_orb, ele, param, end_orb, err_flag, track, mat6, make_matrix)
!
! Subroutine to do tracking using Runge-Kutta integration. 
! The core Runge-Kutta routine used here is odeint_bmad which is a modified version of odeint from Numerical Recipes.
! See the "Numerical Recipes in F90" book.
!
! Input:
!   start_orb   -- Coord_struct: Starting coords.
!   ele         -- Ele_struct
!   param       -- lat_param_struct: Lattice parameters.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
!   bmad_com -- Bmad common block (not an argument).
!     %rel_tol_adaptive_tracking -- Relative tolerance. Default is 1d-8.
!     %abs_tol_adaptive_tracking -- Absolute tolerance. Default is 1d-10.
!     %max_num_runge_kutta_step  -- Maximum number of steps before particle is considered lost.
!
! Output:
!   end_orb     -- coord_struct: Ending coords.
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!   track       -- track_struct, optional: Structure holding the track information.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix propagated through the element.
!- 

subroutine track1_runge_kutta (start_orb, ele, param, end_orb, err_flag, track, mat6, make_matrix)

use runge_kutta_mod, except_dummy => track1_runge_kutta

implicit none

type (coord_struct) :: start_orb, end_orb
type (lat_param_struct), target, intent(inout) :: param
type (ele_struct), target, intent(inout) :: ele
type (track_struct), optional :: track

real(rp), optional :: mat6(6,6)
real(rp) rel_tol, abs_tol, beta_ref, s0_body, s1_body, ds_ref, dref_time

logical err_flag, set_spin
logical, optional :: make_matrix

character(*), parameter :: r_name = 'track1_runge_kutta'

! Runge Kutta is not able to handle a zero length element with a non-zero multipole.
! Nor can RK handle particles at rest (EG e_gun). 

if (ele%key /= patch$ .and. ele%value(l$) == 0) then
  call track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)
  return
endif

end_orb = start_orb

if (start_orb%vec(6) == -1) then
  call out_io (s_error$, r_name, 'Runge Kutta is not able to handle particles with no energy in: ' // ele%name)
  end_orb%state = lost$
  return
endif

! Convert to element coords.
! For a patch, convert to the downstream coords so that the downstream face 
! can be simply described as being at s = s1_body. Additionally, with a patch, s 
! is the distance before the downstream face so the s0_body starting position is 
! negative and the s1_body stopping position is 0.

set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)

if (ele%key == patch$) then
  call track_a_patch (ele, end_orb, .false., s0_body, ds_ref)
  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  if (ele%orientation*end_orb%direction == 1) then
    end_orb%vec(5) = end_orb%vec(5) + (ds_ref + s0_body) * end_orb%beta / beta_ref 
    s0_body = ele%value(l$) + s0_body
    s1_body = ele%value(l$)
  else
    end_orb%vec(5) = end_orb%vec(5) + (ds_ref - s0_body) * end_orb%beta / beta_ref 
    s0_body = s0_body
    s1_body = 0
  endif
else
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., &
                                      set_spin = set_spin, mat6 = mat6, make_matrix = make_matrix)
  if (ele%orientation*end_orb%direction == 1) then
    s0_body = 0; s1_body = ele%value(l$)
  else
    s0_body = ele%value(l$); s1_body = 0
  endif
endif

! Track.
! Note that if ele is a slave, ele%field_calc = refer_to_lords$ and no error message is printed. 

if ((ele%key == lcavity$ .or. ele%key == rfcavity$) .and. &
                ele%field_calc == bmad_standard$ .and. ele%value(l$) < ele%value(l_hard_edge$)) then
  call out_io (s_error$, r_name, 'RUNGE-KUTTA TRACKING THROUGH RF CAVITY: ' // ele%name, &
                          'WILL NOT BE ACCURATE SINCE THE LENGTH IS LESS THAN THE HARD EDGE MODEL LENGTH.')
endif

call odeint_bmad (end_orb, ele, param, s0_body, s1_body, err_flag, track, mat6, make_matrix)
if (err_flag) return

! convert to lab coords.

if (ele%key /= patch$) then
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., &
                                        set_spin = set_spin, mat6 = mat6, make_matrix = make_matrix)
endif

end subroutine
