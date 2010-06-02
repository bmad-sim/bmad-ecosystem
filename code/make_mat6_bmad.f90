!+
! Subroutine make_mat6_bmad (ele, param, c0, c1, end_in, err)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!   end_in -- Logical, optional: If present and True then the end coords c1
!               will be taken as input. not output as normal.
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %vec0  -- 0th order map component
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!   err    -- Logical, optional: Set True if there is an error. False otherwise.
!-

#include "CESR_platform.inc"

subroutine make_mat6_bmad (ele, param, c0, c1, end_in, err)

use bmad, except_dummy => make_mat6_bmad

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: c0, c1
type (coord_struct) :: c00, c11
type (coord_struct) orb
type (lat_param_struct)  param

real(rp), pointer :: mat6(:,:)

real(rp) mat6_m(6,6), mat4(4,4), kmat1(4,4), kmat2(4,4)
real(rp) angle, k1, ks, length, e2, g, g_err
real(rp) k2l, k3l, c2, s2, cs, ks2, del_l
real(rp) factor, kmat6(6,6), drift(6,6)
real(rp) s_pos, s_pos_old, z_slice(100)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp) c_e, c_m, gamma_old, gamma_new
real(rp) arg, rel_p, rel_p2, r11, r12, r21, r22
real(rp) cy, sy, k2, s_off, x_pitch, y_pitch, y_ave, k_z
real(rp) dz_x(3), dz_y(3), ddz_x(3), ddz_y(3), xp_start, yp_start
real(rp) t5_11, t5_14, t5_22, t5_23, t5_33, t5_34, t5_44
real(rp) t1_16, t1_26, t1_36, t1_46, t2_16, t2_26, t2_36, t2_46
real(rp) t3_16, t3_26, t3_36, t3_46, t4_16, t4_26, t4_36, t4_46
real(rp) lcs, lc2s2, k, L
real(rp) cos_phi, gradient, e_start, e_end, e_ratio
real(rp) alpha, sin_a, cos_a, f, phase, E, pxy2, dE0
real(rp) g_tot, b1, rho, ct, st, x, px, y, py, z, pz, Dxy, Dy, px_t
real(rp) Dxy_t, dpx_t, df_dpy, df_dp, kx_1, ky_1, kx_2, ky_2
real(rp) mc2, pc_start, pc_end, p0c_start, p0c_end
real(rp) dcP2_dz1, dbeta1_dpz1, dbeta2_dpz2, beta_start, beta_end
real(rp) dp_coupler, dp_x_coupler, dp_y_coupler

integer i, n_slice, key

logical, optional :: end_in, err
logical err_flag
character(16) :: r_name = 'make_mat6_bmad'

!--------------------------------------------------------
! init

if (present(err)) err = .false.

mat6 => ele%mat6
call mat_make_unit (mat6)
ele%vec0 = 0

length = ele%value(l$)
rel_p = 1 + c0%vec(6)  ! E/E_0
key = ele%key

if (.not. logic_option (.false., end_in)) then
  if (ele%tracking_method == linear$) then
    param%lost = .false.
    call track1_bmad (c0, ele, param, c1)
  else
    call track1 (c0, ele, param, c1)
  endif
  if (param%lost) then
    mat6 = 0
    if (present(err)) err = .true.
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING AT: ' // ele%name)
    return
  endif
endif

c00 = c0
c11 = c1

!--------------------------------------------------------
! drift or element is off or
! Electric Separator or Kicker.

if (.not. ele%is_on .and. key /= lcavity$) key = drift$

if (any (key == (/ drift$, elseparator$, kicker$, rcollimator$, &
        ecollimator$, monitor$, instrument$, hkicker$, vkicker$, pipe$ /) )) then
  call drift_mat6_calc (mat6, length, c0%vec, c1%vec)
  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)
  return
endif

!--------------------------------------------------------
! selection

if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$

select case (key)

!--------------------------------------------------------
! beam-beam interaction

case (beambeam$)

 call offset_particle (ele, param, c00, set$)
 call offset_particle (ele, param, c11, set$, s_pos = length)

  n_slice = nint(ele%value(n_slice$))
  if (n_slice < 1) then
    print *, 'ERROR IN MAKE_MAT6_BMAD: N_SLICE FOR BEAMBEAM ELEMENT IS NEGATIVE'
    call type_ele (ele, .true., 0, .false., 0, .false.)
    stop
  endif

  if (ele%value(charge$) == 0 .or. param%n_part == 0) return

  ! factor of 2 in orb%vec(5) since relative motion of the two beams is 2*c_light

  if (n_slice == 1) then
    call bbi_kick_matrix (ele, param, c00, 0.0_rp, mat6)
  else
    call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)

    s_pos = 0          ! start at IP
    orb = c00
    orb%vec(2) = c00%vec(2) - ele%value(x_pitch_tot$)
    orb%vec(4) = c00%vec(4) - ele%value(y_pitch_tot$)
    call mat_make_unit (mat4)

    do i = 1, n_slice + 1
      s_pos_old = s_pos  ! current position
      s_pos = (z_slice(i) + c00%vec(5)) / 2 ! position of slice relative to IP
      del_l = s_pos - s_pos_old
      mat4(1,1:4) = mat4(1,1:4) + del_l * mat4(2,1:4)
      mat4(3,1:4) = mat4(3,1:4) + del_l * mat4(4,1:4)
      if (i == n_slice + 1) exit
      orb%vec(1) = c00%vec(1) + s_pos * orb%vec(2)
      orb%vec(3) = c00%vec(3) + s_pos * orb%vec(4)
      call bbi_kick_matrix (ele, param, orb, s_pos, kmat6)
      mat4(2,1:4) = mat4(2,1:4) + kmat6(2,1) * mat4(1,1:4) + &
                                  kmat6(2,3) * mat4(3,1:4)
      mat4(4,1:4) = mat4(4,1:4) + kmat6(4,1) * mat4(1,1:4) + &
                                  kmat6(4,3) * mat4(3,1:4)
    enddo

    mat6(1:4,1:4) = mat4

  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Custom

case (custom$)

  print *, 'ERROR IN MAKE_MAT6_BMAD: MAT6_CALC_METHOD = BMAD_STANDARD IS NOT'
  print *, '      ALLOWED FOR A CUSTOM ELEMENT: ', ele%name
  call err_exit

!--------------------------------------------------------
! LCavity: Linac rf cavity
! Ultra-relativistic formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1. The extra factors of beta are included to make the 
! transverse determinant (beta_i*gamma_i)/(beta_f*gamma_f) which it should
! be at low energies.
!
! One must keep in mind that we are NOT using good canonical coordinates since
!   the energy of the reference particle is changing.
! This means that the resulting matrix will NOT be symplectic.
! Since things are very complicated we simplify things by ignoring the
!   off-axis corrections to mat6.
!
! bmad_com%grad_loss_sr_wake is an internal variable used with macroparticles.
!   It should be zero otherwise.

case (lcavity$)

  if (length == 0) return

  f = twopi * ele%value(rf_frequency$) / c_light
  phase = twopi * (ele%value(phi0$) + ele%value(dphi0$) + &
                            ele%value(phi0_err$)) - f * c0%vec(5)

  cos_phi = cos(phase)
  gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos_phi 
  if (.not. ele%is_on) gradient = 0

  if (bmad_com%sr_wakes_on) then
    if (bmad_com%grad_loss_sr_wake /= 0) then  
      ! use grad_loss_sr_wake and ignore e_loss
      gradient = gradient - bmad_com%grad_loss_sr_wake
    else
      gradient = gradient - ele%value(e_loss$) * param%n_part * &
                                          abs(charge_of(param%particle)) / length
    endif
  endif

  if (gradient == 0) then
    call drift_mat6_calc (mat6, length, c0%vec, c1%vec)
    call add_multipoles_and_s_offset
    return
  endif

  p0c_start = ele%value(p0c_start$) 
  pc_start = p0c_start * (1 + c0%vec(6))
  call convert_pc_to (pc_start, param%particle, &
                                    E_tot = e_start, beta = beta_start)
  e_end = e_start + gradient * length
  if (e_end <= 0) then
    call out_io (s_error$, r_name, 'END ENERGY IS NEGATIVE AT ELEMENT: ' // ele%name)
    mat6 = 0   ! garbage.
    err = .true.
    return 
  endif
  call convert_total_energy_to (e_end, param%particle, &
                                           pc = pc_end, beta = beta_end)
  p0c_end = ele%value(p0c$)
  e_ratio = e_end / e_start

  mat6(6,5) = (e_end / pc_end) * ele%value(gradient$) * length * &
                                            f * sin(phase) / ele%value(p0c$) 
  mat6(6,6) = e_end * pc_start * p0c_start / (e_start * pc_end * p0c_end)

  mc2 = mass_of(param%particle)
  dcP2_dz1 = mat6(6,5) * p0c_end
  dbeta1_dpz1 = p0c_start * mc2**2 / e_start**3
  dbeta2_dpz2 = p0c_end * mc2**2 / e_end**3

  mat6(5,5) = beta_end / beta_start + &
              mat6(6,5) * dbeta2_dpz2 * c1%vec(5) / beta_end - &
              beta_end * dcP2_dz1 / gradient + &
              dcP2_dz1 * beta_end * (pc_end - pc_start) * pc_end / &
                                          (length * gradient**2 * e_end)
  mat6(5,6) = -beta_end * c0%vec(5) * dbeta1_dpz1 / beta_start**2 + &
              c1%vec(5) * mat6(6,6) * dbeta2_dpz2 / beta_end - &
              beta_end * (p0c_end * mat6(6,6) - p0c_start) / gradient

  ! entrence kick

  k1 = -gradient / (2 * pc_start)

  ! body...
  ! Note: dimad/liar formulas are:
  !  r11 = 1
  !  r12 = e_start * log (e_ratio) / gradient
  !  r21 = 0
  !  r22 = 1 / e_ratio

  alpha = log(e_ratio) / (2 * sqrt_2 * cos_phi)
  cos_a = cos(alpha)
  sin_a = sin(alpha)
  f = gradient / (2 * sqrt_2 * cos_phi)   ! body matrix
  r11 =  cos_a
  r12 =  sin_a * beta_start * e_start / f
  r21 = -sin_a * f / (e_end * beta_end)
  r22 =  cos_a * beta_start * e_start / (e_end * beta_end)

  ! exit kick

  k2 = +gradient / (2 * pc_end)

  ! put everything together

  mat6(1,1) = r11 + r12*k1
  mat6(1,2) = r12 
  mat6(2,1) = r21 + k2*r11 + k2*r12*k1 + r22*k1
  mat6(2,2) = r22 + k2*r12

  mat6(3:4,3:4) = mat6(1:2,1:2)

  ! off-energy corrections

  mat6(:,2) = mat6(:,2) / (1 + c0%vec(6))
  mat6(:,4) = mat6(:,4) / (1 + c0%vec(6))
  mat6(2,:) = (1 + c1%vec(6)) * mat6(2,:) 
  mat6(4,:) = (1 + c1%vec(6)) * mat6(4,:)

  ! coupler kicks
  ! dp_x_coupler ~ (2,5) matrix term
  ! dp_y_coupler ~ (4,5) matrix term

  if (ele%value(coupler_strength$) /= 0) then

    f = twopi * ele%value(rf_frequency$) / c_light
    dp_coupler = ele%value(gradient$) * ele%value(coupler_strength$) * &
                            f * sin(phase + twopi * ele%value(coupler_phase$))
    dp_x_coupler = dp_coupler * cos (twopi * ele%value(coupler_angle$))
    dp_y_coupler = dp_coupler * sin (twopi * ele%value(coupler_angle$))

    if (nint(ele%value(coupler_at$)) == both_ends$) then
      dp_x_coupler = dp_x_coupler / 2
      dp_y_coupler = dp_y_coupler / 2
    endif

    if (nint(ele%value(coupler_at$)) == entrance_end$ .or. &
        nint(ele%value(coupler_at$)) == both_ends$) then
      mat6(:,5) = mat6(:,5) + &
          (mat6(:,2) * dp_x_coupler + mat6(:,4) * dp_y_coupler) / p0c_start
    endif

    if (nint(ele%value(coupler_at$)) == exit_end$ .or. &
        nint(ele%value(coupler_at$)) == both_ends$) then
      mat6(2,:) = mat6(2,:) + dp_x_coupler * mat6(5,:) / p0c_end
      mat6(4,:) = mat6(4,:) + dp_y_coupler * mat6(5,:) / p0c_end
    endif

  endif

  ! multipoles and s_offset

  if (ele%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Marker, branch, photon_branch

case (marker$, branch$, photon_branch$) 
  return

!--------------------------------------------------------
! Match

case (match$)
  call match_ele_to_mat6 (ele, ele%vec0, ele%mat6, err_flag)
  if (present(err)) err = err_flag

!--------------------------------------------------------
! Mirror

case (mirror$)

  mat6(1, 1) = -1
  mat6(2, 1) =  2 * ele%value(g_graze$) / sin(ele%value(graze_angle$))
  mat6(2, 2) = -1
  mat6(4, 3) = -2 * ele%value(g_trans$)  

  call offset_photon_mat6(mat6, ele)
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Multipole, AB_Multipole

case (multipole$, ab_multipole$)

  if (.not. ele%multipoles_on) return

  call multipole_ele_to_kt (ele, param%particle, knl, tilt, .true.)
  call mat6_multipole (knl, tilt, c0%vec, 1.0_rp, ele%mat6)

  ! if knl(0) is non-zero then the reference orbit itself is bent
  ! and we need to account for this.

  if (knl(0) /= 0) then
    ele%mat6(2,6) = knl(0) * cos(tilt(0))
    ele%mat6(4,6) = knl(0) * sin(tilt(0))
    ele%mat6(5,1) = -ele%mat6(2,6)
    ele%mat6(5,3) = -ele%mat6(4,6)
  endif

  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Octupole
! the octupole is modeled as kick-drift-kick

case (octupole$)

  call offset_particle (ele, param, c00, set$, set_canonical = .false.)
  call offset_particle (ele, param, c11, set$, set_canonical = .false., s_pos = length)

  k3l = ele%value(k3$) * length 
  call mat4_multipole (k3l/2, 0.0_rp, 3, c00%vec, kmat1)
  call mat4_multipole (k3l/2, 0.0_rp, 3, c11%vec, kmat2)

  c00%vec(1:4) = matmul(kmat1, c00%vec(1:4))
  call drift_mat6_calc (drift, length, c00%vec)

  mat6 = drift
  mat6(1:4,1:4) = matmul(kmat2, matmul(drift(1:4,1:4), kmat1))
  mat6(5,1) = drift(5,2) * kmat1(2,1) + drift(5,4) * kmat1(4,1)
  mat6(5,3) = drift(5,2) * kmat1(2,3) + drift(5,4) * kmat1(4,3)
  mat6(2,6) = kmat2(2,1) * drift(1,6) + kmat2(2,3) * drift(3,6)
  mat6(4,6) = kmat2(4,1) * drift(1,6) + kmat2(4,3) * drift(3,6)

  if (ele%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Patch
! If patch_end is set then the xfer map is just the unit matrix

case (patch$) 

  if (ele%value(patch_end$) /= 0) return

  mat6(2,6) = -ele%value(x_pitch_tot$)
  mat6(4,6) = -ele%value(y_pitch_tot$)
  mat6(5,1) =  ele%value(x_pitch_tot$)
  mat6(5,3) =  ele%value(y_pitch_tot$)

  if (ele%value(tilt_tot$) /= 0) then
    cos_a = cos(ele%value(tilt_tot$)) ; sin_a = sin(ele%value(tilt_tot$))
    mat6(1,1) =  cos_a ; mat6(2,2) =  cos_a
    mat6(1,3) =  sin_a ; mat6(2,4) =  sin_a
    mat6(3,1) = -sin_a ; mat6(4,2) = -sin_a
    mat6(3,3) =  cos_a ; mat6(4,4) =  cos_a
  endif

  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! quadrupole

case (quadrupole$)

  call offset_particle (ele, param, c00, set$)
  call offset_particle (ele, param, c11, set$, s_pos = length)

  k1 = ele%value(k1$) / rel_p

  call quad_mat2_calc (-k1, length, mat6(1:2,1:2), dz_x, ddz_x)
  call quad_mat2_calc ( k1, length, mat6(3:4,3:4), dz_y, ddz_y)

  mat6(1,2) = mat6(1,2) / rel_p
  mat6(2,1) = mat6(2,1) * rel_p

  mat6(3,4) = mat6(3,4) / rel_p
  mat6(4,3) = mat6(4,3) * rel_p

  ! The mat6(i,6) terms are constructed so that mat6 is sympelctic

  if (any(c00%vec(1:4) /= 0)) then
    mat6(5,1) = 2 * c00%vec(1) * dz_x(1) +     c00%vec(2) * dz_x(2)
    mat6(5,2) =    (c00%vec(1) * dz_x(2) + 2 * c00%vec(2) * dz_x(3)) / rel_p
    mat6(5,3) = 2 * c00%vec(3) * dz_y(1) +     c00%vec(4) * dz_y(2)
    mat6(5,4) =    (c00%vec(3) * dz_y(2) + 2 * c00%vec(4) * dz_y(3)) / rel_p
    mat6(5,6) = c00%vec(1)**2 * ddz_x(1) / rel_p + &
                c00%vec(1)*c1%vec(2) * ddz_x(2) / rel_p + &
                c00%vec(2)**2 * ddz_x(3) / rel_p + &
                c00%vec(3)**2 * ddz_y(1) / rel_p + &
                c00%vec(3)*c1%vec(4) * ddz_y(2) / rel_p + &
                c00%vec(4)**2 * ddz_y(3) / rel_p 
  endif

  if (any(mat6(5,1:4) /= 0)) then
    mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2)
    mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2)
    mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4)
    mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4)
  endif

  ! tilt and multipoles

  if (ele%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! rbends are not allowed internally

case (rbend$)

  print *, 'ERROR IN MAKE_MAT6_BMAD: RBEND ELEMENTS NOT ALLOWED INTERNALLY!'
  call err_exit

!--------------------------------------------------------
! rf cavity
! Calculation Uses a 3rd order map assuming a linearized rf voltage vs time.

case (rfcavity$)

  if (ele%value(voltage$) == 0) then
    phase = 0
    k = 0
  else
    if (ele%value(RF_frequency$) == 0) then
      print *, 'ERROR IN MAKE_MAT6_BMAD: ', &
                 '"RF_FREQUENCY" ATTRIBUTE NOT SET FOR RF: ', trim(ele%name)
      print *, '      YOU NEED TO SET THIS OR THE "HARMON" ATTRIBUTE.'
      call err_exit
    endif
    factor = twopi * ele%value(rf_frequency$) / c_light
    phase = twopi * (ele%value(phi0$)+ele%value(dphi0$)) + factor * c0%vec(5) 
    k  =  factor * ele%value(voltage$) * cos(phase) / ele%value(p0c$)
  endif

  px = c0%vec(2)
  py = c0%vec(4)
  pz = c0%vec(6)

  dE0 =  ele%value(voltage$) * sin(phase) / ele%value(p0c$)
  L = ele%value(l$)
  E = 1 + pz
  E2 = E**2
  pxy2 = px**2 + py**2

  !

  mat6(1,1) = 1
  mat6(1,2) = L * (1/E - dE0/2 + L*(3*px**2 + py**2)/12 + pz*dE0 + dE0*dE0/3)
  mat6(1,4) = px*py*L**2/6
  mat6(1,5) = L*px * (-k/2 + pz*k + 2*dE0*k/3)
  mat6(1,6) = L*px * (-1/E2 + dE0)
  mat6(2,2) = 1
  mat6(3,2) = px*py*L**2/6
  mat6(3,3) = 1
  mat6(3,4) = L * (1/E - dE0/2 + L*(3*py**2 + px**2)/12 + pz*dE0 + dE0*dE0/3)
  mat6(3,5) = L*py * (-k/2 + pz*k + 2*dE0*k/3)
  mat6(3,6) = L*py * (-1/E2 + dE0)
  mat6(4,4) = 1
  mat6(5,2) = px*L * (-1/E2 + dE0)
  mat6(5,4) = py*L * (-1/E2 + dE0)
  mat6(5,5) = 1 + pxy2*k*L/2
  mat6(5,6) = pxy2 * L / (E2*E)
  mat6(6,2) = k*px*L * (-1/(2*E2) + dE0/3)
  mat6(6,4) = k*py*L * (-1/(2*E2) + dE0/3)
  mat6(6,5) = k * (1 + pxy2*L*k/6)
  mat6(6,6) = 1 + pxy2*k*L/(2*E2*E)

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! sbend

case (sbend$)

  k1 = ele%value(k1$)
  g = ele%value(g$)
  g_err = ele%value(g_err$)
  g_tot = g + g_err

  if (g_tot == 0) then
    call drift_mat6_calc (mat6, length, c0%vec, c1%vec)
    return
  endif

  call offset_particle (ele, param, c00, set$, set_canonical = .false.)
  call offset_particle (ele, param, c11, set$, set_canonical = .false., s_pos = length)
  call track_bend_edge (c00, ele, .true., .false., kx_1, ky_1)
  call track_bend_edge (c11, ele, .false., .true., kx_2, ky_2)  ! track backwards

  ! Body

  if (k1 /= 0) then
    call sbend_body_with_k1_map (g, g_err, length, k1, c00%vec, mat6 = mat6)

  elseif (length /= 0) then

    ! Used: Eqs (12.18) from Etienne Forest: Beam Dynamics.

    x  = c00%vec(1)
    px = c00%vec(2)
    y  = c00%vec(3)
    py = c00%vec(4)
    z  = c00%vec(5)
    pz = c00%vec(6)
 
    b1 = g_tot
    angle = g * length
    rel_p  = 1 + pz
    rel_p2 = rel_p**2

    ct = cos(angle)
    st = sin(angle)

    pxy2 = px**2 + py**2
    Dxy = sqrt(rel_p2 - pxy2)
    Dy  = sqrt(rel_p2 - py**2)

    if (ele%value(g$) == 0) then
      px_t = px - b1 * length
      dpx_t = -px*length - b1
    else
      rho = 1 / ele%value(g$)
      if (pxy2 < 1e-5) then
        f = pxy2 / (2 * rel_p)
        f = pz - f - f*f/2 - g_err*rho - b1*x
      else
        f = sqrt(rel_p2 - pxy2) - 1 - g_err*rho - b1*x
      endif
      px_t = px*ct + f*st
      dpx_t = -px*st/rho + f*ct/rho
    endif

    Dxy_t = sqrt(rel_p2 - px_t**2 - py**2)
    factor = (angle + asin(px/Dy) - asin(px_t/Dy)) / b1
    df_dpy = px/(Dxy*Dy**2) - px_t/(Dxy_t*Dy**2) + st/(Dxy*Dxy_t)
    df_dp = rel_p * (-px/(Dxy*Dy**2) - st/(Dxy*Dxy_t) + px_t/(Dxy_t*Dy**2))

    mat6(1,1) = px_t * st / Dxy_t + ct
    mat6(1,2) = -px_t * (ct - px*st/Dxy) / (b1 * Dxy_t) + st/b1 + px*ct/(b1*Dxy)
    mat6(1,4) = (-py + px_t*py*st/Dxy) / (b1 * Dxy_t) + py*ct/(b1*Dxy)
    mat6(1,6) = (rel_p - px_t*rel_p*st/Dxy) / (b1 * Dxy_t) - rel_p*ct/(b1*Dxy)
    mat6(2,1) = -b1 * st
    mat6(2,2) = ct - px * st / Dxy
    mat6(2,4) = -py * st / Dxy
    mat6(2,6) = rel_p * st / Dxy
    mat6(3,1) = py * st / Dxy_t
    mat6(3,2) = py * (1/ Dxy - ct/Dxy_t + px*st/(Dxy*Dxy_t)) / b1
    mat6(3,3) = 1
    mat6(3,4) = factor + py**2 * df_dpy / b1
    mat6(3,6) = py * df_dp / b1
    mat6(5,1) = -rel_p * st / Dxy_t
    mat6(5,2) = -rel_p * (1/Dxy - ct/Dxy_t + px*st/(Dxy*Dxy_t)) / b1
    mat6(5,4) = -rel_p * py * df_dpy / b1
    mat6(5,6) = -factor - rel_p * df_dp / b1

  endif

  ! edge focusing terms

  mat6(1:6,1) = mat6(1:6,1) + mat6(1:6,2) * kx_1
  mat6(1:6,3) = mat6(1:6,3) + mat6(1:6,4) * ky_1

  mat6(2,1:6) = mat6(2,1:6) + kx_2 * mat6(1,1:6) 
  mat6(4,1:6) = mat6(4,1:6) + ky_2 * mat6(3,1:6)

  if (ele%value(k2$) /= 0) then
    k2l = ele%value(k2$) * ele%value(l$) / 2
    call mat4_multipole (k2l, 0.0_rp, 2, c00%vec, mat6_m(1:4,1:4))
    mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
    mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
    call mat4_multipole (k2l, 0.0_rp, 2, c11%vec, mat6_m(1:4,1:4))
    mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)
  endif

  if (ele%value(tilt_tot$)+ele%value(roll$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$)+ele%value(roll$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Sextupole.
! the sextupole is modeled as kick-drift-kick

case (sextupole$)


  call offset_particle (ele, param, c00, set$, set_canonical = .false.)
  call offset_particle (ele, param, c11, set$, set_canonical = .false., s_pos = length)

  k2l = ele%value(k2$) * length 
  call mat4_multipole (k2l/2, 0.0_rp, 2, c00%vec, kmat1)
  call mat4_multipole (k2l/2, 0.0_rp, 2, c11%vec, kmat2)

  c00%vec(1:4) = matmul(kmat1, c00%vec(1:4))
  call drift_mat6_calc (drift, length, c00%vec)

  mat6 = drift
  mat6(1:4,1:4) = matmul(kmat2, matmul(drift(1:4,1:4), kmat1))
  mat6(5,1) = drift(5,2) * kmat1(2,1) + drift(5,4) * kmat1(4,1)
  mat6(5,3) = drift(5,2) * kmat1(2,3) + drift(5,4) * kmat1(4,3)
  mat6(2,6) = kmat2(2,1) * drift(1,6) + kmat2(2,3) * drift(3,6)
  mat6(4,6) = kmat2(4,1) * drift(1,6) + kmat2(4,3) * drift(3,6)

  if (ele%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! solenoid

case (solenoid$)

  call offset_particle (ele, param, c00, set$)
  call offset_particle (ele, param, c11, set$, s_pos = length)

  ks = ele%value(ks$) / rel_p

  call solenoid_mat_calc (ks, length, mat6(1:4,1:4))

  mat6(1,2) = mat6(1,2) / rel_p
  mat6(1,4) = mat6(1,4) / rel_p

  mat6(2,1) = mat6(2,1) * rel_p
  mat6(2,3) = mat6(2,3) * rel_p

  mat6(3,2) = mat6(3,2) / rel_p
  mat6(3,4) = mat6(3,4) / rel_p

  mat6(4,1) = mat6(4,1) * rel_p
  mat6(4,3) = mat6(4,3) * rel_p


  c2 = mat6(1,1)
  s2 = mat6(1,4) * ks / 2
  cs = mat6(1,3)

  lcs = length * cs
  lc2s2 = length * (c2 - s2) / 2

  t1_16 =  lcs * ks
  t1_26 = -lc2s2 * 2
  t1_36 = -lc2s2 * ks
  t1_46 = -lcs * 2

  t2_16 =  lc2s2 * ks**2 / 2
  t2_26 =  lcs * ks
  t2_36 =  lcs * ks**2 / 2
  t2_46 = -lc2s2 * ks

  t3_16 =  lc2s2 * ks
  t3_26 =  lcs * 2
  t3_36 =  lcs * ks
  t3_46 = -lc2s2 * 2

  t4_16 = -lcs * ks**2 / 2
  t4_26 =  lc2s2 * ks
  t4_36 =  t2_16
  t4_46 =  lcs * ks

  arg = length / 2
  t5_11 = -arg * (ks/2)**2
  t5_14 =  arg * ks
  t5_22 = -arg
  t5_23 = -arg * ks
  t5_33 = -arg * (ks/2)**2
  t5_44 = -arg

  ! the mat6(i,6) terms are constructed so that mat6 is sympelctic

  mat6(5,1) =  2 * c00%vec(1) * t5_11 + c00%vec(4) * t5_14
  mat6(5,2) = (2 * c00%vec(2) * t5_22 + c00%vec(3) * t5_23) / rel_p
  mat6(5,3) =  2 * c00%vec(3) * t5_33 + c00%vec(2) * t5_23
  mat6(5,4) = (2 * c00%vec(4) * t5_44 + c00%vec(1) * t5_14) / rel_p

  mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2) + &
                  mat6(5,4) * mat6(1,3) - mat6(5,3) * mat6(1,4)
  mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2) + &
                  mat6(5,4) * mat6(2,3) - mat6(5,3) * mat6(2,4)
  mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4) + &
                  mat6(5,2) * mat6(3,1) - mat6(5,1) * mat6(3,2)
  mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4) + &
                  mat6(5,2) * mat6(4,1) - mat6(5,1) * mat6(4,2)

  ! mat6(5,6) 

  ks2 = ele%value(ks$) / 2
  xp_start = (c0%vec(2) + ks2 * c0%vec(3)) 
  yp_start = (c0%vec(4) - ks2 * c0%vec(1)) 
  mat6(5,6) = length * (xp_start**2 + yp_start**2 ) / rel_p**3

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! solenoid/quad

case (sol_quad$)

  call offset_particle (ele, param, c00, set$)
  call offset_particle (ele, param, c11, set$, s_pos = length)

  call sol_quad_mat6_calc (ele%value(ks$), ele%value(k1$), length, mat6, c00%vec)

  if (ele%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! taylor

case (taylor$)

  call make_mat6_taylor (ele, param, c0)

!--------------------------------------------------------
! wiggler

case (wiggler$)

  call offset_particle (ele, param, c00, set$)
  call offset_particle (ele, param, c11, set$, s_pos = length)

  call mat_make_unit (mat6)     ! make a unit matrix

  if (length == 0) then
    call add_multipoles_and_s_offset
    return
  endif

  k1 = -0.5 * (c_light * ele%value(b_max$) / &
                  (ele%value(p0c$) * rel_p))**2

  ! octuple correction to k1

  y_ave = (c00%vec(3) + c11%vec(3)) / 2
  if (ele%value(l_pole$) == 0) then
    k_z = 0
  else
    k_z = pi / ele%value(l_pole$)
  endif
  k1 = k1 * (1 + 2 * (k_z * y_ave)**2)

  !

  mat6(1, 1) = 1
  mat6(1, 2) = length
  mat6(2, 1) = 0
  mat6(2, 2) = 1

  call quad_mat2_calc (k1, length, mat6(3:4,3:4))

  cy = mat6(3, 3)
  sy = mat6(3, 4)

  t5_22 = -length / 2
  t5_33 =  k1 * (length - sy*cy) / 4
  t5_34 = -k1 * sy**2 / 2
  t5_44 = -(length + sy*cy) / 4

  ! the mat6(i,6) terms are constructed so that mat6 is sympelctic

  mat6(5,2) = 2 * c00%vec(2) * t5_22
  mat6(5,3) = 2 * c00%vec(3) * t5_33 +     c00%vec(4) * t5_34
  mat6(5,4) =     c00%vec(3) * t5_34 + 2 * c00%vec(4) * t5_44

  mat6(1,6) = mat6(5,2) * mat6(1,1)
  mat6(2,6) = mat6(5,2) * mat6(2,1)
  mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4)
  mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4)

  if (ele%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  call add_multipoles_and_s_offset
  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! Accelerating solenoid with steerings
! WARNING: This 6x6 matrix may produce bad results at low energies!

case (accel_sol$)

!  if ((ele%value(s_st1$) < 0.) .or.  &
!      (ele%value(s_st1$) + ele%value(l_st1$) > ele%value(s_st2$)) .or.  &
!      (ele%value(s_st2$) + ele%value(l_st2$) > length)) then
!    print *, 'ERROR IN MAKE_MAT6_BMAD: STEERINGS MUST NOT OVERLAP AND MUST BE',  &
!      ' CONTAINED WITHIN'
!    print *, 'THE ACCEL_SOL ELEMENT!'
!    call type_ele(ele, .true., 0, .false., 0, .false.)
!    stop
!  endif
!
!  call mat_make_unit (mat6)     ! make a unit matrix
!  if (ele%value(voltage$) /= 0) then
!    if (ele%value(rf_wavelength$) == 0) then
!      print *, 'ERROR IN MAKE_MAT6_BMAD: RF IS ON BUT "RF_WAVELENGTH" NOT SET',  &
!            ' IN ACCEL_SOL!'
!      call err_exit
!    else
!      phase = twopi * (ele%value(phi0$)+ele%value(dphi0$)) 
!      mat6(6,5) = ele%value(voltage$) * cos(phase) *  &
!                    twopi / ele%value(rf_wavelength$) /ele%value(p0c$)
!      c_e = ele%value(voltage$) * sin(phase) / (mass_of(param%particle) * length)
!    endif
!  else
!    c_e = 0.0
!  endif
!  c_m = param%particle * c_light * ele%value(b_z$) / mass_of(param%particle)
!  gamma_old = ele%value(p0c$) * rel_p / mass_of(param%particle)
!  gamma_new = gamma_old + c_e * length
!    call accel_sol_mat_calc (length, c_m, c_e, gamma_old, gamma_new, &
!                                    0.0_rp, 0.0_rp, c00%vec, mat4, vec_st)
!  mat4 = mat6(1:4,1:4)
!
!  call add_multipoles_and_s_offset
!  ele%vec0 = c1%vec - matmul(mat6, c0%vec)

!--------------------------------------------------------
! unrecognized element

case default

  print *, 'ERROR IN MAKE_MAT6_BMAD: UNKNOWN ELEMENT KEY:', ele%key
  print *, '      FOR ELEMENT: ', ele%name
  call err_exit

end select

!--------------------------------------------------------
! put in multipole components

contains

subroutine add_multipoles_and_s_offset

if (associated(ele%a_pole) .and. key /= multipole$ .and. key /= ab_multipole$) then
  mat6_m = 0
  call multipole_ele_to_kt (ele, param%particle, knl, tilt, .true.)
  call mat6_multipole (knl, tilt, c0%vec, 0.5_rp, mat6_m)
  mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
  mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
  call mat6_multipole (knl, tilt, c1%vec, 0.5_rp, mat6_m)
  mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
  mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)
endif

if (ele%value(s_offset_tot$) /= 0) then
  s_off = ele%value(s_offset_tot$)
  mat6(1,:) = mat6(1,:) - s_off * mat6(2,:)
  mat6(3,:) = mat6(3,:) - s_off * mat6(4,:)
  mat6(:,2) = mat6(:,2) + mat6(:,1) * s_off
  mat6(:,4) = mat6(:,4) + mat6(:,3) * s_off
endif

! pitch corrections

call mat6_add_pitch (ele, ele%mat6)

end subroutine

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine bbi_kick_matrix (ele, param, orb, s_pos, mat6)

use bmad_struct
use bmad_interface, except_dummy => bbi_kick_matrix

implicit none

type (ele_struct)  ele
type (coord_struct)  orb
type (lat_param_struct) param

real(rp) x_pos, y_pos, del, sig_x, sig_y, coef, garbage, s_pos
real(rp) ratio, k0_x, k1_x, k0_y, k1_y, mat6(6,6), beta, bbi_const

!

call mat_make_unit (mat6)

sig_x = ele%value(sig_x$)
sig_y = ele%value(sig_y$)

if (sig_x == 0 .or. sig_y == 0) return

if (s_pos /= 0 .and. ele%a%beta /= 0) then
  beta = ele%a%beta - 2 * ele%a%alpha * s_pos + ele%a%gamma * s_pos**2
  sig_x = sig_x * sqrt(beta / ele%a%beta)
  beta = ele%b%beta - 2 * ele%b%alpha * s_pos + ele%b%gamma * s_pos**2
  sig_y = sig_y * sqrt(beta / ele%b%beta)
endif

x_pos = orb%vec(1) / sig_x  ! this has offset in it
y_pos = orb%vec(3) / sig_y

del = 0.001

ratio = sig_y / sig_x
call bbi_kick (x_pos, y_pos, ratio, k0_x, k0_y)
call bbi_kick (x_pos+del, y_pos, ratio, k1_x, garbage)
call bbi_kick (x_pos, y_pos+del, ratio, garbage, k1_y)

bbi_const = -param%n_part * ele%value(charge$) * classical_radius_factor /  &
                    (2 * pi * ele%value(p0c$) * (sig_x + sig_y))

coef = bbi_const / (ele%value(n_slice$) * del * (1 + orb%vec(6)))

mat6(2,1) = coef * (k1_x - k0_x) / sig_x
mat6(4,3) = coef * (k1_y - k0_y) / sig_y

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+      
! Subroutine offset_photon_mat6 (mat6, start, ele)
!
! Subroutine to transform a 6x6 transfer matrix to a new reference frame
! with the given offsets, pitches and tilts of the given element.
!
! Modules needed:
!   use bmad
!
! Input:
!   mat6(6,6) -- Real(rp): Untilted matrix.
!   ele       -- Ele_struct: Mirror or equivalent.
!   start     -- Coord_struct: Starting coords
!
! Output:
!   mat6(6,6) -- Real(rp): Tilted matrix.
!-

subroutine offset_photon_mat6 (mat6, ele)

use bmad_struct
use bmad_interface

implicit none

type (ele_struct), target :: ele

real(rp) mat6(6,6), mm(6,6)
real(rp), pointer :: p(:)
real(rp), save :: old_tilt = 0, ct, st
real(rp) c2g, s2g, offset(6), tilt, graze
real(rp) off(3), rot(3), project_x(3), project_y(3), project_s(3)

!

p => ele%value

! Set: Work backward form element mat6 matrix...

! Set: Graze angle error

if (p(graze_angle_err$) /= 0) then
  mat6(:,1) = mat6(:,1) + mat6(:,1) * p(graze_angle_err$) 
endif

! Set: Tilt

tilt = p(tilt_tot$) + p(tilt_err$)

if (tilt /= 0) then

  if (tilt /= old_tilt) then
    ct = cos(tilt)
    st = sin(tilt)
    old_tilt = tilt
  endif

  mm(:,1) = mat6(:,1) * ct - mat6(:,3) * st
  mm(:,2) = mat6(:,2) * ct - mat6(:,4) * st
  mm(:,3) = mat6(:,3) * ct + mat6(:,1) * st
  mm(:,4) = mat6(:,4) * ct + mat6(:,2) * st
  mm(:,5) = mat6(:,5)
  mm(:,6) = mat6(:,6)

else
  mm = mat6
endif

! Set: transverse offsets and pitches

mm(:,1) = mm(:,1) + mm(:,5) * p(x_pitch_tot$) 
mm(:,3) = mm(:,3) + mm(:,5) + p(y_pitch_tot$)

! Set: s_offset

mm(:,2) = mm(:,2) + mm(:,1) * p(s_offset_tot$)
mm(:,4) = mm(:,4) + mm(:,3) * p(s_offset_tot$)

!------------------------------------------------------
! Unset: 

c2g = cos(2*p(graze_angle$)) 
s2g = sin(2*p(graze_angle$))

if (p(tilt_err$) /= 0) then
  ct = cos(p(tilt$)) 
  st = sin(p(tilt$))
  old_tilt = p(tilt$)
endif

project_x = (/ c2g * ct**2 + st**2, -ct * st + c2g * ct * st, -ct * s2g /)
project_y = (/ -ct * st + c2g * ct * st, ct**2 + c2g * st**2, -s2g * st /) 
project_s = (/ ct * s2g, s2g * st, c2g /)

! Unset: graze_angle_error

if (p(graze_angle_err$) /= 0) then
  mm(5,:) = mm(5,:) + p(graze_angle_err$) * mm(1,:)
endif

! Unset tilt

if (p(tilt$) /= 0) then
  mat6(1,:) = ct * mm(1,:) - st * mm(3,:)
  mat6(2,:) = ct * mm(2,:) - st * mm(4,:)
  mat6(3,:) = ct * mm(3,:) + st * mm(1,:)
  mat6(4,:) = ct * mm(4,:) + st * mm(2,:)
  mat6(5,:) =     mm(5,:)
  mat6(6,:) =     mm(6,:)
else
  mat6 = mm
endif

! Unset: tilt_err

if (p(tilt_err$) /= 0) then
  rot = project_s * p(tilt_err$)

  ct = cos(rot(3)) 
  st = sin(rot(3))
  old_tilt = rot(3)

  mm(1,:) = ct * mat6(1,:) - st * mat6(3,:)
  mm(2,:) = ct * mat6(2,:) - st * mat6(4,:)
  mm(3,:) = ct * mat6(3,:) + st * mat6(1,:)
  mm(4,:) = ct * mat6(4,:) + st * mat6(2,:)
  mm(5,:) =     mat6(5,:)
  mm(6,:) =     mat6(6,:)

  mm(5,:) = mm(5,:) - rot(2) * mm(2,:) + rot(1) * mm(3,:)

  mat6 = mm

endif

! Unset pitch

rot = project_x * p(y_pitch_tot$) - project_y * p(x_pitch_tot$)
mat6(5,:) = mat6(5,:) + rot(2) * mat6(2,:) - rot(1) * mat6(3,:)

! Unset: offset

off = project_x * p(x_offset_tot$) + project_y * p(y_offset_tot$) + project_s * p(s_offset_tot$)

mat6(1,:) = mat6(1,:) - off(3) * mat6(2,:)
mat6(3,:) = mat6(3,:) - off(3) * mat6(4,:)

end subroutine

