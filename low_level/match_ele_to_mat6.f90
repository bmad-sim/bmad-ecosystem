!+
! Subroutine match_ele_to_mat6 (ele, start_orb, mat6, vec0, err_flag, twiss_ele, include_delta_time)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of a match element. 
!
! Note: ele%taylor%term will be deallocated if the xfer map has changed.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele                -- ele_struct: Match element.
!   start_orb          -- coord_struct: Starting orbit.
!   twiss_ele          -- ele_struct, optional: If present and ele%value(match_end$) is True, 
!                           use the Twiss parameters in this element as the upstream end Twiss.
!   include_delta_time -- logical, optional: If False, ignore any finite ele%value(delta_time$). 
!                           Default is True.
!
! Output:
!   vec0(6)            -- Real(rp): 0th order part of the transfer map.
!   mat6(6,6)          -- Real(rp): Transfer matrix (1st order part of xfer map).
!   err_flag           -- Logical: Set true if there is an error. False otherwise.
!-

subroutine match_ele_to_mat6 (ele, start_orb, mat6, vec0, err_flag, twiss_ele, include_delta_time)

use bmad_routine_interface, dummy => match_ele_to_mat6

implicit none

type (ele_struct), target :: ele, ele0, ele1
type (coord_struct) start_orb
type (ele_struct), optional, target :: twiss_ele
type (ele_struct), pointer :: t_ele

real(rp), pointer :: v(:)
real(rp) orb0(6), orb1(6), p0c, mat6(6,6), vec0(6)
real(rp) old_mat6(6,6), old_vec0(6), E, m56, mass

integer species

logical :: err_flag
logical, optional :: include_delta_time

! Error Check

old_mat6 = mat6
old_vec0 = vec0

if (ele%value(beta_a1$) == 0 .or. ele%value(beta_b1$) == 0) then
  mat6 = 0
  vec0 = 0
  err_flag = .true.
  call kill_taylor (ele%taylor)
  return
endif

! Match_end

if (is_true(ele%value(match_end$))) then
  if (present(twiss_ele)) then
    t_ele => twiss_ele
  else
    t_ele => pointer_to_next_ele (ele, -1)
  endif
  ele%value(beta_a0$)    = t_ele%a%beta
  ele%value(beta_b0$)    = t_ele%b%beta
  ele%value(alpha_a0$)   = t_ele%a%alpha
  ele%value(alpha_b0$)   = t_ele%b%alpha
  ele%value(eta_x0$)     = t_ele%x%eta
  ele%value(eta_y0$)     = t_ele%y%eta
  ele%value(etap_x0$)    = t_ele%x%etap
  ele%value(etap_y0$)    = t_ele%y%etap
  ele%value(gamma_c$)    = t_ele%gamma_c
  ele%value(c_11$:c_22$) = [t_ele%c_mat(1,1), t_ele%c_mat(1,2), t_ele%c_mat(2,1), t_ele%c_mat(2,2)]
endif

! Special case where match_end is set but there is no beginning beta value yet.
! In this case, just return the unit matrix. 
! This is not an error since it is important for lat_make_mat6 to keep on computing matrices.

if (is_true(ele%value(match_end$)) .and. (ele%value(beta_a0$) == 0 .or. ele%value(beta_b0$) == 0)) then
  call mat_make_unit (mat6)
  vec0 = 0
  err_flag = .false.
  call kill_taylor (ele%taylor)
  return
endif

!

if (ele%value(beta_a0$) == 0 .or. ele%value(beta_b0$) == 0) then
  mat6 = 0
  err_flag = .true.
  call kill_taylor (ele%taylor)
  return
endif

!

err_flag = .false.

v => ele%value

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

ele0%c_mat(1,:) = [v(c_11$), v(c_12$)]
ele0%c_mat(2,:) = [v(c_21$), v(c_22$)]
ele0%gamma_c    = v(gamma_c$)

ele1%c_mat = 0 
ele1%gamma_c = 1

ele0%name = ele%name
ele1%name = ele%name

orb0 = [v(x0$), v(px0$), v(y0$), v(py0$), v(z0$), v(pz0$)]
orb1 = [v(x1$), v(px1$), v(y1$), v(py1$), v(z1$), v(pz1$)]

call transfer_mat_from_twiss (ele0, ele1, orb0, orb1, mat6)

! Finite delta_time introduces a mat(5,6) component.

if (logic_option(.true., include_delta_time) .and. ele%value(delta_time$) /= 0) then
  species = ele%branch%param%particle
  p0c = ele%value(p0c$)
  mass = mass_of(species)
  E = (1 + start_orb%vec(6)) * ele%value(p0c$) / start_orb%beta
  m56 = -c_light * ele%value(delta_time$) * mass**2 *p0c / E**3
  mat6(:,6) = mat6(:,6) + mat6(:,5) * m56 
endif

! Kick part

vec0 = orb1 - matmul (mat6, orb0)

if (any(mat6 /= old_mat6) .or. any(vec0 /= old_vec0)) call kill_taylor (ele%taylor)

end subroutine match_ele_to_mat6
