!+
! Subroutine track1_runge_kutta (orbit, ele, param, err_flag, track, mat6, make_matrix)
!
! Subroutine to do tracking using Runge-Kutta integration. 
! The core Runge-Kutta routine used here is odeint_bmad which is a modified version of odeint from Numerical Recipes.
! See the "Numerical Recipes in F90" book.
!
! Input:
!   orbit       -- Coord_struct: Starting coords.
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
!   orbit       -- coord_struct: Ending coords.
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!   track       -- track_struct, optional: Structure holding the track information.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix propagated through the element.
!- 

subroutine track1_runge_kutta (orbit, ele, param, err_flag, track, mat6, make_matrix)

use runge_kutta_mod, except_dummy => track1_runge_kutta

implicit none

type (coord_struct) :: orbit
type (lat_param_struct), target :: param
type (ele_struct), target :: ele
type (track_struct), optional :: track

real(rp), optional :: mat6(6,6)
real(rp) rel_tol, abs_tol, beta_ref, s0_body, s1_body, ds_ref, dref_time, length

logical err_flag, set_spin
logical, optional :: make_matrix

character(*), parameter :: r_name = 'track1_runge_kutta'

! Runge Kutta is not able to handle a zero length element with a non-zero multipole.
! Nor can RK handle particles at rest (EG e_gun). 

if (ele%key /= patch$ .and. ele%value(l$) == 0) then
  call track_a_zero_length_element (orbit, ele, param, err_flag, track)
  return
endif

if (orbit%vec(6) == -1) then
  call out_io (s_error$, r_name, 'Runge Kutta is not able to handle particles with no energy in: ' // ele%name)
  orbit%state = lost$
  return
endif

! Convert to element coords.
! For a patch, convert to the downstream coords so that the downstream face 
! can be simply described as being at s = s1_body. Additionally, with a patch, s 
! is the distance before the downstream face so the s0_body starting position is 
! negative and the s1_body stopping position is 0.

set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)

if (ele%key == patch$) then
  call track_a_patch (ele, orbit, .false., s0_body, ds_ref)
  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  if (ele%orientation*orbit%direction*orbit%time_dir == 1) then
    orbit%vec(5) = orbit%vec(5) + orbit%time_dir * (ds_ref + s0_body) * orbit%beta / beta_ref 
    length = patch_length(ele, exit_end$)
    s0_body = length + s0_body
    s1_body = length
  else
    orbit%vec(5) = orbit%vec(5) + orbit%time_dir * (ds_ref - s0_body) * orbit%beta / beta_ref 
    s0_body = s0_body
    s1_body = 0
  endif

else
  if (ele%bookkeeping_state%has_misalign) call offset_particle (ele, set$, orbit, set_hvkicks = .false., &
                                      set_spin = set_spin, mat6 = mat6, make_matrix = make_matrix)
  if (ele%orientation*orbit%direction*orbit%time_dir == 1) then
    s0_body = 0; s1_body = ele%value(l$)
  else
    s0_body = ele%value(l$); s1_body = 0
  endif
endif

! Track.
! Note that if ele is a slave, ele%field_calc = refer_to_lords$ and no error message is printed. 

if ((ele%key == lcavity$ .or. ele%key == rfcavity$) .and. ele%field_calc == bmad_standard$ .and. &
                ele%value(l$) < ele%value(l_active$) - 0.1*bmad_com%significant_length) then
  call out_io (s_error$, r_name, 'RUNGE-KUTTA TRACKING THROUGH A CAVITY: ' // ele%name, &
                          'WILL NOT BE ACCURATE SINCE THE LENGTH IS LESS THAN THE HARD EDGE MODEL LENGTH.')
endif

call odeint_bmad (orbit, ele, param, s0_body, s1_body, err_flag, track, mat6, make_matrix)
if (err_flag) return

! convert to lab coords.

if (ele%key /= patch$) then
  if (ele%bookkeeping_state%has_misalign) call offset_particle (ele, unset$, orbit, set_hvkicks = .false., &
                                        set_spin = set_spin, mat6 = mat6, make_matrix = make_matrix)
endif

end subroutine
