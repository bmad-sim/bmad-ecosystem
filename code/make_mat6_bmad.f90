!+
! Subroutine make_mat6_bmad (ele, param, c0, c1, end_in)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- Param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!   end_in -- Logical, optional: If present and True then the end coords c1
!               will be taken as input. not output as normal.
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!-

#include "CESR_platform.inc"

subroutine make_mat6_bmad (ele, param, c0, c1, end_in)

  use bmad

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0, c1
  type (coord_struct) :: c00, c11, c0t, c1t
  type (coord_struct) orb
  type (param_struct)  param

  real(rp), pointer :: mat6(:,:)

  real(rp) mat6_m(6,6), mat2(2,2), mat4(4,4), kmat1(4,4), kmat2(4,4)
  real(rp) e1, e2, angle, g, cos_angle, sin_angle, k1, ks, length, kc
  real(rp) phi, k2l, k3l, c2, s2, cs, ks2, del_l, g_bend, l_period, l_bend
  real(rp) factor, l_drift, dx, kmat6(6,6), drift(6,6)
  real(rp) s_pos, s_pos_old, z_slice(100)
  real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rp) r, c_e, c_m, gamma_old, gamma_new, vec_st(4)
  real(rp) sqrt_k, arg, kick2, rel_E, dE, r11, r12, r21, r22
  real(rp) cx, sx, cy, sy, k2l_2, k2l_3, k2l_4, k2
  real(rp) x_off, y_off, s_off, x_pit, y_pit, y_ave, k_z, del_x, del_y
  real(rp) t5_11, t5_12, t5_22, t5_33, t5_34, t5_44, t5_14, t5_23
  real(rp) t1_16, t1_26, t1_36, t1_46, t2_16, t2_26, t2_36, t2_46
  real(rp) t3_16, t3_26, t3_36, t3_46, t4_16, t4_26, t4_36, t4_46
  real(rp) lcs, lc2s2, error, rho, px, py, pz, k, L, k_gradient
  real(rp) cos_phi, gradient, e_start, e_end, e_ratio
  real(rp) alpha, sin_a, cos_a, f, phase, E, pxy2, dE0

  integer i, n, n_slice, n_pole, key

  logical, optional :: end_in

!--------------------------------------------------------
! init

  length = ele%value(l$)
  mat6 => ele%mat6
  rel_E = 1 + c0%vec(6)  ! E/E_0

  if (present(end_in)) then
    if (.not. end_in) call track1 (c0, ele, param, c1)
  else
    call track1 (c0, ele, param, c1)
  endif

  c00 = c0
  c11 = c1

!--------------------------------------------------------
! drift or element is off or
! Electric Separator or Kicker.

  key = ele%key
  if (.not. ele%is_on) key = drift$
  if (any (key == (/ drift$, elseparator$, kicker$, rcollimator$, &
          ecollimator$, monitor$, instrument$, hkicker$, vkicker$ /) )) then

    call drift_mat6_calc (mat6, length, c0%vec, c1%vec)

    goto 8000   ! put in multipole ends if needed
  else
    call mat_make_unit (mat6)
  endif

!--------------------------------------------------------
! marker

  if (ele%key == marker$) return

!--------------------------------------------------------
! Put in offsets, etc.
! Note: c0t and c1t are the coords in the frame of reference where the element
!       is upright (element frame of reference).

  call offset_particle (ele, param, c00, set$, set_canonical = .false.)
  call offset_particle (ele, param, c11, set$, &
                                set_canonical = .false., s_pos = length)

  c0t = c00
  c1t = c11

  if (ele%value(tilt$) /= 0) then
    if (ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
      call tilt_coords (ele%value(tilt$), c0t%vec, .true.)
      call tilt_coords (ele%value(tilt$), c1t%vec, .true.)
    endif
  endif

!--------------------------------------------------------
! selection

  key = ele%key
  if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$

  select case (key)

!--------------------------------------------------------
! Patch

  case (patch$) 

    mat6(2,6) = -ele%value(x_pitch$)
    mat6(4,6) = -ele%value(y_pitch$)
    mat6(5,1) =  ele%value(x_pitch$)
    mat6(5,3) =  ele%value(y_pitch$)

    if (ele%value(tilt$) /= 0) then
      cos_a = cos(ele%value(tilt$)) ; sin_a = sin(ele%value(tilt$))
      mat6(1,1) =  cos_a ; mat6(2,2) =  cos_a
      mat6(1,3) =  sin_a ; mat6(2,4) =  sin_a
      mat6(3,1) = -sin_a ; mat6(4,2) = -sin_a
      mat6(3,3) =  cos_a ; mat6(4,4) =  cos_a
    endif

!--------------------------------------------------------
! sbend

  case (sbend$)

    e1 = ele%value(e1$)
    e2 = ele%value(e2$)
    k1 = ele%value(k1$)
    g = ele%value(g$) + ele%value(delta_g$)

    if (length == 0) return

    if (g == 0) then
      call drift_mat6_calc (mat6, length, c0%vec, c1%vec)
      return
    endif

    angle = length * g
    e1 = e1 + c0t%vec(2) / rel_E
    e2 = e2 - c1t%vec(2) / rel_E
    angle = angle + (c0t%vec(2) - c1t%vec(2)) / rel_E
    g = g / rel_E
    length = angle / g
    k1 = k1 / rel_E

    kc = g**2 + k1

    call quad_mat_calc (-kc, length, mat6(1:2,1:2))
    call quad_mat_calc (k1, length, mat6(3:4,3:4))

    phi = sqrt(abs(kc)) * length
    if (kc < 0.0) then
      mat6(1,6) = (1 - cosh(phi)) * g / kc
      mat6(2,6) = sinh(phi) * g / sqrt(-kc)
      mat6(5,6) = (phi - sinh(phi)) * g**2 / abs(kc)**1.5
    else
      mat6(1,6) = (1 - cos(phi)) * g / kc
      mat6(2,6) = sin(phi) * g/ sqrt(kc)
      mat6(5,6) = (sin(phi) - phi) * g**2 / kc**1.5
    endif

! Above correct for (x, x', y, y') coord system. 
! Actually we have (x, p_x, y, p_y) so transform.

    mat6(1,6) = mat6(1,6) / rel_E - mat6(1,2) * c0t%vec(2) / rel_E**2
    mat6(2,6) = rel_E * mat6(2,6) 
    mat6(3,6) = -c0t%vec(4) * (mat6(5,6) + length)

    mat6(1,2) = mat6(1,2) / rel_E
    mat6(2,1) = mat6(2,1) * rel_E

    mat6(3,4) = mat6(3,4) / rel_E
    mat6(4,3) = mat6(4,3) * rel_E

!

    mat6(5,1) = mat6(1,6) * mat6(2,1) - mat6(2,6) * mat6(1,1) 
    mat6(5,2) = mat6(1,6) * mat6(2,2) - mat6(2,6) * mat6(1,2)
  
    mat6(5,4) = -c0t%vec(4) * (mat6(5,6) + length)

    if (e1 /= 0) then
      arg = tan(e1) * g
      mat6(1:6,1) = mat6(1:6,1) + mat6(1:6,2) * arg
      mat6(1:6,3) = mat6(1:6,3) - mat6(1:6,4) * arg
    endif

    if (e2 /= 0) then
      arg = tan(e2) * g
      mat6(2,1:6) = mat6(2,1:6) + mat6(1,1:6) * arg
      mat6(4,1:6) = mat6(4,1:6) - mat6(3,1:6) * arg
    endif

    if (ele%value(tilt$)+ele%value(roll$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$)+ele%value(roll$))
    endif

!--------------------------------------------------------
! quadrupole

  case (quadrupole$)


    k1 = ele%value(k1$) / rel_E

    if (k1 == 0) then
      call drift_mat6_calc (mat6, length, c0%vec, c1%vec)
      goto 8000
    endif

    call quad_mat_calc (-k1, length, mat6(1:2,1:2))
    call quad_mat_calc ( k1, length, mat6(3:4,3:4))

    mat6(1,2) = mat6(1,2) / rel_E
    mat6(2,1) = mat6(2,1) * rel_E

    mat6(3,4) = mat6(3,4) / rel_E
    mat6(4,3) = mat6(4,3) * rel_E

! The mat6(i,6) terms are constructed so that mat6 is sympelctic

    if (any(c0t%vec(1:4) /= 0)) then

      cx = mat6(1, 1)
      sx = mat6(1, 2)
      cy = mat6(3, 3)
      sy = mat6(3, 4)

      t5_11 =  k1 * (cx*sx - length) / 4
      t5_12 =  k1 * sx**2 / 2
      t5_22 = -(length + sx*cx) / 4
      t5_33 =  k1 * (length - sy*cy) / 4
      t5_34 = -k1 * sy**2 / 2
      t5_44 = -(length + sy*cy) / 4

      mat6(5,1) = 2 * c0t%vec(1) * t5_11 +     c0t%vec(2) * t5_12
      mat6(5,2) =     c0t%vec(1) * t5_12 + 2 * c0t%vec(2) * t5_22
      mat6(5,3) = 2 * c0t%vec(3) * t5_33 +     c0t%vec(4) * t5_34
      mat6(5,4) =     c0t%vec(3) * t5_34 + 2 * c0t%vec(4) * t5_44

      mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2)
      mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2)
      mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4)
      mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4)

    endif

! tilt

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! Sextupole.
! the sextupole is modeled as kick-drift-kick

  case (sextupole$)

    k2l = ele%value(k2$) * length 
    call mat4_multipole (k2l/2, 0.0_rp, 2, c0t%vec, kmat1)
    call mat4_multipole (k2l/2, 0.0_rp, 2, c1t%vec, kmat2)

    c0t%vec(1:4) = matmul(kmat1, c0t%vec(1:4))
    call drift_mat6_calc (drift, length, c0t%vec)

    mat6 = drift
    mat6(1:4,1:4) = matmul(kmat2, matmul(drift(1:4,1:4), kmat1))
    mat6(5,1) = drift(5,2) * kmat1(2,1) + drift(5,4) * kmat1(4,1)
    mat6(5,3) = drift(5,2) * kmat1(2,3) + drift(5,4) * kmat1(4,3)
    mat6(2,6) = kmat2(2,1) * drift(1,6) + kmat2(2,3) * drift(3,6)
    mat6(4,6) = kmat2(4,1) * drift(1,6) + kmat2(4,3) * drift(3,6)

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! octupole
! the octupole is modeled as kick-drift-kick

  case (octupole$)

    k3l = ele%value(k3$) * length 
    call mat4_multipole (k3l/2, 0.0_rp, 3, c0t%vec, kmat1)
    call mat4_multipole (k3l/2, 0.0_rp, 3, c1t%vec, kmat2)

    c0t%vec(1:4) = matmul(kmat1, c0t%vec(1:4))
    call drift_mat6_calc (drift, length, c0t%vec)

    mat6 = drift
    mat6(1:4,1:4) = matmul(kmat2, matmul(drift(1:4,1:4), kmat1))
    mat6(5,1) = drift(5,2) * kmat1(2,1) + drift(5,4) * kmat1(4,1)
    mat6(5,3) = drift(5,2) * kmat1(2,3) + drift(5,4) * kmat1(4,3)
    mat6(2,6) = kmat2(2,1) * drift(1,6) + kmat2(2,3) * drift(3,6)
    mat6(4,6) = kmat2(4,1) * drift(1,6) + kmat2(4,3) * drift(3,6)

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! solenoid

  case (solenoid$)

    ks = ele%value(ks$) / rel_E

    call solenoid_mat_calc (ks, length, mat6(1:4,1:4))

    mat6(1,2) = mat6(1,2) / rel_E
    mat6(1,4) = mat6(1,4) / rel_E

    mat6(2,1) = mat6(2,1) * rel_E
    mat6(2,3) = mat6(2,3) * rel_E

    mat6(3,2) = mat6(3,2) / rel_E
    mat6(3,4) = mat6(3,4) / rel_E

    mat6(4,1) = mat6(4,1) * rel_E
    mat6(4,3) = mat6(4,3) * rel_E


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

    mat6(5,1) = 2 * c0t%vec(1) * t5_11 + c0t%vec(4) * t5_14
    mat6(5,2) = 2 * c0t%vec(2) * t5_22 + c0t%vec(3) * t5_23
    mat6(5,3) = 2 * c0t%vec(3) * t5_33 + c0t%vec(2) * t5_23
    mat6(5,4) = 2 * c0t%vec(4) * t5_44 + c0t%vec(1) * t5_14

    mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2) + &
                    mat6(5,4) * mat6(1,3) - mat6(5,3) * mat6(1,4)
    mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2) + &
                    mat6(5,4) * mat6(2,3) - mat6(5,3) * mat6(2,4)
    mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4) + &
                    mat6(5,2) * mat6(3,1) - mat6(5,1) * mat6(3,2)
    mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4) + &
                    mat6(5,2) * mat6(4,1) - mat6(5,1) * mat6(4,2)

!--------------------------------------------------------
! linac rf cavity
!
! One must keep in mind that we are NOT using good canonical coordinates since
!   the energy of the reference particle is changing.
! This means that the resulting matrix will NOT be symplectic.
! Since things are very complicated we simplify things by ignoring the
!   off-axis corrections to mat6.

  case (lcavity$)

    f = twopi * ele%value(rf_frequency$) / c_light
    phase = f * c0%vec(5) + twopi * ele%value(phi0$)
    mat6(6,5) = -ele%value(gradient$) * ele%value(l$) * f * sin(phase)

    cos_phi = cos(phase)
    gradient = ele%value(gradient$) * cos_phi
    k_gradient = -ele%value(k_loss$) * abs(param%charge) 

    if (bmad_com%emulate_liar_bug) then
      gradient = gradient + k_gradient
      k_gradient = 0
    endif

    if (gradient == 0 .and. k_gradient == 0) then
      call drift_mat6_calc (mat6, length, c0%vec, c1%vec)
      goto 8000  ! put in mulipole ends if needed
    endif

! entrence kick

    e_start = ele%value(energy_start$) * (1 + c0%vec(6))
    k1 = -gradient / (2 * e_start)

! body 

    e_start = e_start + (k_gradient * ele%value(l$))/2
    e_end = e_start + gradient * ele%value(l$)
    e_ratio = e_end / e_start

    if (bmad_com%use_dimad_lcavity) then  ! use dimad formula
      r11 = 1
      r12 = e_start * log (e_ratio) / gradient
      r21 = 0
      r22 = 1 / e_ratio
    else
      alpha = log(e_ratio) / (2 * sqrt_2 * cos_phi)
      cos_a = cos(alpha)
      sin_a = sin(alpha)
      f = gradient / (2 * sqrt_2 * cos_phi)   ! body matrix
      r11 =  cos_a
      r12 =  sin_a * e_start / f
      r21 = -sin_a * f / e_end
      r22 =  cos_a * e_start / e_end
    endif

! exit kick

    k2 = +gradient / (2 * e_end)

! put everything together

    mat6(1,1) = r11 + r12*k1
    mat6(1,2) = r12 
    mat6(2,1) = r21 + k2 + k2*r12*k1 + r22*k1
    mat6(2,2) = r22 + k2*r12

    mat6(3:4,3:4) = mat6(1:2,1:2)

! off-energy corrections

    mat6(:,2) = mat6(:,2) / (1 + c0%vec(6))
    mat6(:,4) = mat6(:,4) / (1 + c0%vec(6))
    mat6(2,:) = (1 + c1%vec(6)) * mat6(2,:) 
    mat6(4,:) = (1 + c1%vec(6)) * mat6(4,:)

!--------------------------------------------------------
! rf cavity
! Calculation Uses a 3rd order map assuming a linearized rf voltage vs time.

  case (rfcavity$)

    if (ele%value(volt$) == 0) then
      phase = 0
      k = 0
    else
      if (ele%value(RF_wavelength$) == 0) then
        print *, 'ERROR IN MAKE_MAT6_BMAD: ', &
                   '"RF_WAVELENGTH" ATTRIBUTE NOT SET FOR RF: ', trim(ele%name)
        print *, '      YOU NEED TO SET THIS OR THE "HARMON" ATTRIBUTE.'
        call err_exit
      endif
      phase = twopi * (ele%value(phi0$) + c0%vec(5) / ele%value(rf_wavelength$))
      k  =  twopi * ele%value(volt$) * cos(phase) / &
                              (param%beam_energy * ele%value(rf_wavelength$))
    endif

    px = c0%vec(2)
    py = c0%vec(4)
    pz = c0%vec(6)

    dE0 =  ele%value(volt$) * sin(phase) / param%beam_energy
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

!--------------------------------------------------------
! beam-beam interaction

  case (beambeam$)

    n_slice = nint(ele%value(n_slice$))
    if (n_slice < 1) then
      print *, 'ERROR IN MAKE_MAT6_BMAD: N_SLICE FOR BEAMBEAM ELEMENT IS NEGATIVE'
      call type_ele (ele, .true., 0, .false., 0, .false.)
      stop
    endif

    if (ele%value(charge$) == 0 .or. param%n_part == 0) return

! factor of 2 in orb%vec(5) since relative motion of the two beams is 2*c_light

    if (n_slice == 1) then
      call bbi_kick_matrix (ele, c0t, 0.0_rp, mat6)
    else
      call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)

      s_pos = 0          ! start at IP
      orb%vec(2) = c0t%vec(2) - ele%value(x_pitch$)
      orb%vec(4) = c0t%vec(4) - ele%value(y_pitch$)
      call mat_make_unit (mat4)

      do i = 1, n_slice + 1
        s_pos_old = s_pos  ! current position
        s_pos = (z_slice(i) + c0t%vec(5)) / 2 ! position of slice relative to IP
        del_l = s_pos - s_pos_old
        mat4(1,1:4) = mat4(1,1:4) + del_l * mat4(2,1:4)
        mat4(3,1:4) = mat4(3,1:4) + del_l * mat4(4,1:4)
        if (i == n_slice + 1) exit
        orb%vec(1) = c0t%vec(1) + s_pos * orb%vec(2)
        orb%vec(3) = c0t%vec(3) + s_pos * orb%vec(4)
        call bbi_kick_matrix (ele, orb, s_pos, kmat6)
        mat4(2,1:4) = mat4(2,1:4) + kmat6(2,1) * mat4(1,1:4) + &
                                    kmat6(2,3) * mat4(3,1:4)
        mat4(4,1:4) = mat4(4,1:4) + kmat6(4,1) * mat4(1,1:4) + &
                                    kmat6(4,3) * mat4(3,1:4)
      enddo

      mat6(1:4,1:4) = mat4

    endif

!--------------------------------------------------------
! wiggler

  case (wiggler$)

    call mat_make_unit (mat6)     ! make a unit matrix
    if (length == 0) goto 8000

    k1 = ele%value(k1$) / rel_E**2

! octuple correction to k1

    y_ave = (c0t%vec(3) + c1t%vec(3)) / 2
    k_z = pi * ele%value(n_pole$) / length
    k1 = k1 * (1 + 2 * (k_z * y_ave)**2)

! correction for fact that wigglers with odd number of poles have end
! poles with modified bending radius

!    n_pole = nint(ele%value(n_pole$))
!    if (mod(n_pole, 2) == 1) then
!      rho = 4 * ele%value(rho$) / pi
!      l_period = length / n_pole
!      l_bend = 8 * l_period / pi**2
!      l_drift = l_period - l_bend
!      factor = sqrt(rho**2 - (l_bend/2)**2)
!      dx = 2 * (rho - factor) + l_drift * l_bend / (2 * factor)
!      factor = -2 * rho * dx / (l_bend * (length-l_bend/2))
!      k1 = k1 * (1 + factor / n_pole)
!    endif

!

    mat6(1, 1) = 1
    mat6(1, 2) = length
    mat6(2, 1) = 0
    mat6(2, 2) = 1

    call quad_mat_calc (k1, length, mat6(3:4,3:4))

    cy = mat6(3, 3)
    sy = mat6(3, 4)

    t5_22 = -length / 2
    t5_33 =  k1 * (length - sy*cy) / 4
    t5_34 = -k1 * sy**2 / 2
    t5_44 = -(length + sy*cy) / 4

! the mat6(i,6) terms are constructed so that mat6 is sympelctic

    mat6(5,2) = 2 * c0t%vec(2) * t5_22
    mat6(5,3) = 2 * c0t%vec(3) * t5_33 +     c0t%vec(4) * t5_34
    mat6(5,4) =     c0t%vec(3) * t5_34 + 2 * c0t%vec(4) * t5_44

    mat6(1,6) = mat6(5,2) * mat6(1,1)
    mat6(2,6) = mat6(5,2) * mat6(2,1)
    mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4)
    mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4)

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! solenoid/quad

  case (sol_quad$)

    ks = ele%value(ks$) / rel_E
    k1 = ele%value(k1$) / rel_E

    call sol_quad_mat6_calc (ks, k1, length, mat6, c0t%vec)

    mat6(1,2) = mat6(1,2) / rel_E
    mat6(1,4) = mat6(1,4) / rel_E

    mat6(2,1) = mat6(2,1) * rel_E
    mat6(2,3) = mat6(2,3) * rel_E

    mat6(3,2) = mat6(3,2) / rel_E
    mat6(3,4) = mat6(3,4) / rel_E

    mat6(4,1) = mat6(4,1) * rel_E
    mat6(4,3) = mat6(4,3) * rel_E

    mat6(5,1) = -mat6(2,6)*mat6(1,1) + mat6(1,6)*mat6(2,1) - mat6(4,6)*mat6(3,1) + mat6(3,6)*mat6(4,1)
    mat6(5,2) = -mat6(2,6)*mat6(1,2) + mat6(1,6)*mat6(2,2) - mat6(4,6)*mat6(3,2) + mat6(3,6)*mat6(4,2)
    mat6(5,3) = -mat6(2,6)*mat6(1,3) + mat6(1,6)*mat6(2,3) - mat6(4,6)*mat6(3,3) + mat6(3,6)*mat6(4,3)
    mat6(5,4) = -mat6(2,6)*mat6(1,4) + mat6(1,6)*mat6(2,4) - mat6(4,6)*mat6(3,4) + mat6(3,6)*mat6(4,4)

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! multipole

  case (multipole$, ab_multipole$)
    if (.not. ele%multipoles_on) return
    call multipole_ele_to_kt (ele, param%particle, knl, tilt, .true.)
    call mat6_multipole (knl, tilt, c0t%vec, 1.0_rp, ele%mat6)

!--------------------------------------------------------
! accelerating solenoid with steerings
! WARNING: This 6x6 matrix may produce bad results at low energies!

  case (accel_sol$)

    if ((ele%value(s_st1$) < 0.) .or.  &
        (ele%value(s_st1$) + ele%value(l_st1$) > ele%value(s_st2$)) .or.  &
        (ele%value(s_st2$) + ele%value(l_st2$) > length)) then
      print *, 'ERROR IN MAKE_MAT6_BMAD: STEERINGS MUST NOT OVERLAP AND MUST BE',  &
        ' CONTAINED WITHIN'
      print *, 'THE ACCEL_SOL ELEMENT!'
      call type_ele(ele, .true., 0, .false., 0, .false.)
      stop
    endif

    call mat_make_unit (mat6)     ! make a unit matrix
    if (ele%value(volt$) /= 0) then
      if (ele%value(rf_wavelength$) == 0) then
        print *, 'ERROR IN MAKE_MAT6_BMAD: RF IS ON BUT "RF_WAVELENGTH" NOT SET',  &
              ' IN ACCEL_SOL!'
        call err_exit
      else
        mat6(6,5) = ele%value(volt$) * cos(twopi*ele%value(phi0$)) *  &
                      twopi / ele%value(rf_wavelength$) /param%beam_energy
        c_e = ele%value(volt$) * sin(twopi*ele%value(phi0$))  &
              / (m_electron * length)
      endif
    else
      c_e = 0.0
    endif
    c_m = param%particle * c_light * ele%value(b_z$) / m_electron
    gamma_old = param%beam_energy * rel_E / m_electron
    gamma_new = gamma_old + c_e * length
!!    call accel_sol_mat_calc (length, c_m, c_e, gamma_old, gamma_new, &
!!                                    0.0_rp, 0.0_rp, c0t%vec, mat4, vec_st)
    mat4 = mat6(1:4,1:4)

!--------------------------------------------------------
! rbends are not allowed internally

  case (rbend$)

    print *, 'ERROR IN MAKE_MAT6_BMAD: RBEND ELEMENTS NOT ALLOWED INTERNALLY!'
    call err_exit

!--------------------------------------------------------
! Custom

  case (custom$)

    print *, 'ERROR IN MAKE_MAT6_BMAD: MAT6_CALC_METHOD = BMAD_STANDARD IS NOT'
    print *, '      ALLOWED FOR A CUSTOM ELEMENT: ', ele%name
    call err_exit

!--------------------------------------------------------
! unrecognized element

  case default

    print *, 'ERROR IN MAKE_MAT6_BMAD: UNKNOWN ELEMENT KEY:', ele%key
    print *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

!--------------------------------------------------------
! put in multipole components

8000 continue

  if (associated(ele%a) .and. key /= multipole$ .and. key /= ab_multipole$) then
    mat6_m = 0
    call multipole_ele_to_kt (ele, param%particle, knl, tilt, .true.)
    call mat6_multipole (knl, tilt, c00%vec, 0.5_rp, mat6_m)
    mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
    mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
    call mat6_multipole (knl, tilt, c11%vec, 0.5_rp, mat6_m)
    mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)
  endif

  if (ele%value(s_offset$) /= 0) then
    s_off = ele%value(s_offset$)
    mat6(1,:) = mat6(1,:) - s_off * mat6(2,:)
    mat6(3,:) = mat6(3,:) - s_off * mat6(4,:)
    mat6(:,2) = mat6(:,2) + mat6(:,1) * s_off
    mat6(:,4) = mat6(:,4) + mat6(:,3) * s_off
  endif

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


subroutine bbi_kick_matrix (ele, orb, s_pos, mat6)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct)  ele
  type (coord_struct)  orb

  real(rp) x_pos, y_pos, del, sig_x, sig_y, coef, garbage, s_pos
  real(rp) ratio, k0_x, k1_x, k0_y, k1_y, mat6(6,6), beta

!

  sig_x = ele%value(sig_x$)
  sig_y = ele%value(sig_y$)

  if (s_pos /= 0 .and. ele%x%beta /= 0) then
    beta = ele%x%beta - 2 * ele%x%alpha * s_pos + ele%x%gamma * s_pos**2
    sig_x = sig_x * sqrt(beta / ele%x%beta)
    beta = ele%y%beta - 2 * ele%y%alpha * s_pos + ele%y%gamma * s_pos**2
    sig_y = sig_y * sqrt(beta / ele%y%beta)
  endif


  x_pos = orb%vec(1) / sig_x  ! this has offset in it
  y_pos = orb%vec(3) / sig_y

  del = 0.001

  ratio = sig_y / sig_x
  call bbi_kick (x_pos, y_pos, ratio, k0_x, k0_y)
  call bbi_kick (x_pos+del, y_pos, ratio, k1_x, garbage)
  call bbi_kick (x_pos, y_pos+del, ratio, garbage, k1_y)

  coef = ele%value(bbi_const$) / (1 + orb%vec(6))

  call mat_make_unit (mat6)
  mat6(2,1) = coef * (k1_x - k0_x) / (ele%value(n_slice$) * del * sig_x)
  mat6(4,3) = coef * (k1_y - k0_y) / (ele%value(n_slice$) * del * sig_y)

end subroutine
