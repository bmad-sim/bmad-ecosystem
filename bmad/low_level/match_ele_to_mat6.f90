!+
! Subroutine match_ele_to_mat6 (ele, start_orb, mat6, vec0, err_flag, twiss_ele, include_delta_time, set_trombone)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of a match element. 
!
! Note: ele%taylor%term will be deallocated if the xfer map has changed.
!
! Input:
!   ele                 -- ele_struct: Match element.
!   start_orb           -- coord_struct: Starting orbit.
!   twiss_ele           -- ele_struct, optional: If present and ele%value(match_end$) is True, 
!                            use the Twiss parameters in this element as the upstream end Twiss.
!   include_delta_time  -- logical, optional: If False, ignore any finite ele%value(delta_time$). 
!                            Default is True.
!   set_trombone        -- logical, optional: Default is False. If True, set the beginning and ending Twiss
!                           values in the element to create a phase trombone.
!
! Output:
!   vec0(6)             -- Real(rp): 0th order part of the transfer map.
!   mat6(6,6)           -- Real(rp): Transfer matrix (1st order part of xfer map).
!   err_flag            -- Logical: Set true if there is an error. False otherwise.
!-

subroutine match_ele_to_mat6 (ele, start_orb, mat6, vec0, err_flag, twiss_ele, include_delta_time, set_trombone)

use equal_mod, dummy => match_ele_to_mat6
use taylor_mod, only: kill_taylor

implicit none

type (ele_struct), target :: ele, ele0, ele1
type (coord_struct) start_orb
type (ele_struct), optional, target :: twiss_ele
type (ele_struct), pointer :: t_ele

real(rp), pointer :: v(:)
real(rp) orb0(6), orb1(6), p0c, mat6(6,6), vec0(6)
real(rp) old_mat6(6,6), old_vec0(6), E, m56, mass

integer species

logical :: err_flag, setit
logical, optional :: include_delta_time, set_trombone

!

v => ele%value

if (present(twiss_ele)) then
  t_ele => twiss_ele
else
  t_ele => pointer_to_next_ele (ele, -1)
endif

! Phase_trombone = True and set_trombone = False means that the trombone is not ready to be set.
! In this case, just set the trombone to the unit matrix.

setit = logic_option(.false., set_trombone)

if (is_true(ele%value(phase_trombone$)) .and. .not. setit) then
  call mat_make_unit(mat6)
  vec0 = [v(x1$)-v(x0$), v(px1$)-v(px0$), v(y1$)-v(y0$), v(py1$)-v(py0$), v(z1$)-v(z0$), v(pz1$)-v(pz0$)]
  err_flag = .false.
  return
endif

! If no Twiss has been set then just set mat6 to a rotation matrix if the Twiss parameters are available.
! If not available, set mat6 to the unit matrix

if ((v(beta_a0$) == 0 .and. v(beta_b0$) == 0 .and. v(alpha_a0$) == 0 .and. v(alpha_b0$) == 0 .and. &
      v(beta_a1$) == 0 .and. v(beta_b1$) == 0 .and. v(alpha_a1$) == 0 .and. v(alpha_b1$) == 0) .or. setit) then
  if ((t_ele%a%beta > 0 .and. t_ele%b%beta > 0) .or. setit) then
    v(beta_a0$)  = t_ele%a%beta
    v(beta_b0$)  = t_ele%b%beta
    v(alpha_a0$) = t_ele%a%alpha
    v(alpha_b0$) = t_ele%b%alpha
    v(eta_x0$)   = t_ele%x%eta
    v(eta_y0$)   = t_ele%y%eta
    v(etap_x0$)  = t_ele%x%etap
    v(etap_y0$)  = t_ele%y%etap
    v(c11_mat0$:c22_mat0$) = [t_ele%c_mat(1,1), t_ele%c_mat(1,2), t_ele%c_mat(2,1), t_ele%c_mat(2,2)]
    v(mode_flip0$) = int_logic(t_ele%mode_flip)
    v(beta_a1$)  = t_ele%a%beta
    v(beta_b1$)  = t_ele%b%beta
    v(alpha_a1$) = t_ele%a%alpha
    v(alpha_b1$) = t_ele%b%alpha
    v(eta_x1$)   = t_ele%x%eta
    v(eta_y1$)   = t_ele%y%eta
    v(etap_x1$)  = t_ele%x%etap
    v(etap_y1$)  = t_ele%y%etap
    v(c11_mat1$:c22_mat1$) = [t_ele%c_mat(1,1), t_ele%c_mat(1,2), t_ele%c_mat(2,1), t_ele%c_mat(2,2)]
    v(mode_flip1$) = int_logic(t_ele%mode_flip)
  else
    call mat_make_unit(mat6)
    vec0 = [v(x1$)-v(x0$), v(px1$)-v(px0$), v(y1$)-v(y0$), v(py1$)-v(py0$), v(z1$)-v(z0$), v(pz1$)-v(pz0$)]
    err_flag = .false.
    return
  endif
endif

! Error Check. Negative beta can be caused via twiss_propagate1 with a non-symplectic transfer matrix.

old_mat6 = mat6
old_vec0 = vec0

if (v(beta_a1$) <= 0 .or. v(beta_b1$) <= 0) then
  mat6 = 0
  vec0 = 0
  err_flag = .true.
  call kill_taylor (ele%taylor)
  return
endif

! Match_end

if (is_true(v(match_end$))) then
  v(beta_a0$)  = t_ele%a%beta
  v(beta_b0$)  = t_ele%b%beta
  v(alpha_a0$) = t_ele%a%alpha
  v(alpha_b0$) = t_ele%b%alpha
  v(eta_x0$)   = t_ele%x%eta
  v(eta_y0$)   = t_ele%y%eta
  v(etap_x0$)  = t_ele%x%etap
  v(etap_y0$)  = t_ele%y%etap
  v(c11_mat0$:c22_mat0$) = [t_ele%c_mat(1,1), t_ele%c_mat(1,2), t_ele%c_mat(2,1), t_ele%c_mat(2,2)]
  v(mode_flip0$) = int_logic(t_ele%mode_flip)
endif

! Special case where match_end is set but there is no beginning beta value yet.
! In this case, just return the unit matrix. 
! This is not an error since it is important for lat_make_mat6 to keep on computing matrices.

if (is_true(v(match_end$)) .and. (v(beta_a0$) == 0 .or. v(beta_b0$) == 0)) then
  call mat_make_unit (mat6)
  vec0 = 0
  err_flag = .false.
  call kill_taylor (ele%taylor)
  return
endif

!

if (v(beta_a0$) <= 0 .or. v(beta_b0$) <= 0) then
  mat6 = 0
  err_flag = .true.
  call kill_taylor (ele%taylor)
  return
endif

!

err_flag = .false.

ele0%a%beta   = v(beta_a0$)
ele0%a%alpha  = v(alpha_a0$)
ele0%a%phi    = 0
ele0%x%eta    = v(eta_x0$)
ele0%x%etap   = v(etap_x0$)

ele0%b%beta   = v(beta_b0$)
ele0%b%alpha  = v(alpha_b0$)
ele0%b%phi    = 0
ele0%y%eta    = v(eta_y0$)
ele0%y%etap   = v(etap_y0$)

ele1%a%beta   = v(beta_a1$)
ele1%a%alpha  = v(alpha_a1$)
ele1%a%phi    = v(dphi_a$)
ele1%x%eta    = v(eta_x1$)
ele1%x%etap   = v(etap_x1$)

ele1%b%beta   = v(beta_b1$)
ele1%b%alpha  = v(alpha_b1$)
ele1%b%phi    = v(dphi_b$)
ele1%y%eta    = v(eta_y1$)
ele1%y%etap   = v(etap_y1$)

ele0%c_mat(1,:) = [v(c11_mat0$), v(c12_mat0$)]
ele0%c_mat(2,:) = [v(c21_mat0$), v(c22_mat0$)]
ele0%gamma_c    = sqrt(1 - determinant(ele0%c_mat))
ele0%mode_flip  = is_true(ele%value(mode_flip0$))

ele1%c_mat(1,:) = [v(c11_mat1$), v(c12_mat1$)]
ele1%c_mat(2,:) = [v(c21_mat1$), v(c22_mat1$)]
ele1%gamma_c    = sqrt(1 - determinant(ele1%c_mat))
ele1%mode_flip  = is_true(ele%value(mode_flip1$))

ele0%name = ele%name
ele1%name = ele%name

orb0 = [v(x0$), v(px0$), v(y0$), v(py0$), v(z0$), v(pz0$)]
orb1 = [v(x1$), v(px1$), v(y1$), v(py1$), v(z1$), v(pz1$)]

call transfer_mat_from_twiss (ele0, ele1, orb0, orb1, mat6)

! Finite delta_time introduces a mat(5,6) component.

if (logic_option(.true., include_delta_time) .and. v(delta_time$) /= 0) then
  species = ele%branch%param%particle
  p0c = v(p0c$)
  mass = mass_of(species)
  E = (1 + start_orb%vec(6)) * v(p0c$) / start_orb%beta
  m56 = -c_light * v(delta_time$) * mass**2 *p0c / E**3
  mat6(:,6) = mat6(:,6) + mat6(:,5) * m56 
endif

! Kick part

vec0 = orb1 - matmul (mat6, orb0)

if (any(mat6 /= old_mat6) .or. any(vec0 /= old_vec0)) call kill_taylor (ele%taylor)

end subroutine match_ele_to_mat6
