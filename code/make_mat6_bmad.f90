!+
! Subroutine make_mat6_bmad (ele, param, c0, c1)
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
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!-

#include "CESR_platform.inc"

subroutine make_mat6_bmad (ele, param, c0, c1)

  use bmad

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0, c1
  type (coord_struct) :: c00, c11
  type (coord_struct) orb
  type (param_struct)  param

  real(rdef), pointer :: mat6(:,:)

  real(rdef) mat6_m(6,6), mat2(2,2), mat4(4,4), kmat1(4,4), kmat2(4,4)
  real(rdef) e1, e2, angle, g, cos_angle, sin_angle, k1, ks, length, kc
  real(rdef) phi, k2l, k3l, c2, s2, cs, ks2, del_l, g_bend, l_period, l_bend
  real(rdef) factor, l_drift, dx, kmat6(6,6)
  real(rdef) s_pos, s_pos_old, z_slice(100)
  real(rdef) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rdef) r, c_e, c_m, gamma_old, gamma_new, vec_st(4)
  real(rdef) sqrt_k, arg, kick2
  real(rdef) cx, sx, cy, sy, k2l_2, k2l_3, k2l_4
  real(rdef) x_off, y_off, s_off, x_pit, y_pit, y_ave, k_z, del_x, del_y
  real(rdef) t5_11, t5_12, t5_22, t5_33, t5_34, t5_44, t5_14, t5_23
  real(rdef) t1_16, t1_26, t1_36, t1_46, t2_16, t2_26, t2_36, t2_46
  real(rdef) t3_16, t3_26, t3_36, t3_46, t4_16, t4_26, t4_36, t4_46
  real(rdef) lcs, lc2s2, error, rho

  integer i, n, n_slice, n_pole, key

!--------------------------------------------------------
! init

  length = ele%value(l$)
  mat6 => ele%mat6

  call track1 (c0, ele, param, c1)

!--------------------------------------------------------
! drift or element is off or
! Electric Separator or Kicker.

  if (ele%key == drift$ .or. (.not. ele%is_on) .or. &
      ele%key == elseparator$ .or. ele%key == kicker$) then

    orb%vec = (c0%vec + c1%vec) / 2
    call drift_mat6_calc (mat6, length, orb%vec)

    goto 8000   ! put in multipole ends if needed
  else
    call mat_make_unit (mat6)
  endif

!--------------------------------------------------------
! marker

  if (ele%key == marker$) return

!--------------------------------------------------------
! Put in offsets, etc.
! Note: c00 and c11 are the coords in the frame of reference where the element
!       is upright (no tilt).

  if (ele%value(x_offset$) /= 0 .or. ele%value(y_offset$) /= 0 .or. &
      ele%value(s_offset$) /= 0 .or. ele%value(x_pitch$) /= 0 .or. &
      ele%value(y_pitch$) /= 0) then
    x_off = ele%value(x_offset$)
    y_off = ele%value(y_offset$)
    x_pit = ele%value(x_pitch$) * length / 2
    y_pit = ele%value(y_pitch$) * length / 2
    s_off = ele%value(s_offset$)
    del_x = s_off * c0%x%vel / (1 + c0%z%vel)
    del_y = s_off * c0%y%vel / (1 + c0%z%vel)

    c00%vec = c0%vec - (/ x_off - x_pit*length/2 + del_x, x_pit, &
             y_off - y_pit*length/2 + del_y, y_pit, 0.0_rdef, 0.0_rdef /) 
    c11%vec = c1%vec - (/ x_off + x_pit*length/2 - del_x, x_pit, &
             y_off + y_pit*length/2 - del_y, y_pit, 0.0_rdef, 0.0_rdef /) 
  else
    c00 = c0
    c11 = c1
  endif

  if (ele%value(tilt$) /= 0) then
    if (ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
      call tilt_coords (ele%value(tilt$), c00%vec, .true.)
      call tilt_coords (ele%value(tilt$), c11%vec, .true.)
    endif
  endif
               
!--------------------------------------------------------
! selection

  key = ele%key
  if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$

  select case (key)

!--------------------------------------------------------
! sbend

  case (sbend$) 

    e1 = ele%value(e1$)
    e2 = ele%value(e2$)
    k1 = ele%value(k1$)

    if (length == 0) return

    if (ele%value(g$) == 0) then
      call drift_mat6_calc (mat6, length, orb%vec)
      return
    endif

    angle = ele%value(l$) * ele%value(g$)
    e1 = e1 + c00%x%vel
    e2 = e2 - c11%x%vel
    angle = angle + c00%x%vel - c11%x%vel
    g = ele%value(g$) / (1 + c00%z%vel)
    length = angle / g
    k1 = k1 / (1 + c00%z%vel)
   
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

    mat6(5,1) = -mat6(2,6)
    mat6(5,2) = -mat6(1,6)

    mat6(5,4) = -c00%vec(4) * (mat6(5,6) + length)
    mat6(3,6) = -c00%vec(4) * (mat6(5,6) + length)

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

    
    k1 = ele%value(k1$) / (1 + c00%z%vel)

    if (k1 == 0) then
      mat6(1,2)  =  length
      mat6(3,4)  =  length
      return
    endif

    call quad_mat_calc (-k1, length, mat6(1:2,1:2))
    call quad_mat_calc ( k1, length, mat6(3:4,3:4))

! The mat6(i,6) terms are constructed so that mat6 is sympelctic

    if (any(c00%vec(1:4) /= 0)) then

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

      mat6(5,1) = 2 * c00%vec(1) * t5_11 +     c00%vec(2) * t5_12
      mat6(5,2) =     c00%vec(1) * t5_12 + 2 * c00%vec(2) * t5_22
      mat6(5,3) = 2 * c00%vec(3) * t5_33 +     c00%vec(4) * t5_34
      mat6(5,4) =     c00%vec(3) * t5_34 + 2 * c00%vec(4) * t5_44
  
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
   
    k2l = ele%value(k2$) * length / (1 + c00%z%vel)
    call mat4_multipole (k2l/2, 0.0_rdef, 2, c00%vec, kmat1)
    mat4 = kmat1
    mat4(1,1:4) = kmat1(1,1:4) + length * kmat1(2,1:4) ! kick * length
    mat4(3,1:4) = kmat1(3,1:4) + length * kmat1(4,1:4)
    call mat4_multipole (k2l/2, 0.0_rdef, 2, c11%vec, kmat2)
    mat6(1:4,1:4) = matmul (kmat2, mat4)

    c00%vec(1:4) = matmul(kmat1, c00%vec(1:4))

    mat6(1,6) = -(kmat2(1,1)*c00%vec(2) + kmat2(1,3)*c00%vec(4)) * length
    mat6(2,6) = -(kmat2(2,1)*c00%vec(2) + kmat2(2,3)*c00%vec(4)) * length
    mat6(3,6) = -(kmat2(3,1)*c00%vec(2) + kmat2(3,3)*c00%vec(4)) * length
    mat6(4,6) = -(kmat2(4,1)*c00%vec(2) + kmat2(4,3)*c00%vec(4)) * length

    mat6(5,1) = -(c00%vec(2)*kmat2(2,1) + c00%vec(4)*kmat2(4,1)) * length
    mat6(5,2) = -(c00%vec(2)*kmat2(2,2) + c00%vec(4)*kmat2(4,2)) * length
    mat6(5,3) = -(c00%vec(2)*kmat2(2,3) + c00%vec(4)*kmat2(4,3)) * length
    mat6(5,4) = -(c00%vec(2)*kmat2(2,4) + c00%vec(4)*kmat2(4,4)) * length

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! octupole
! the octupole is modeled as kick-drift-kick

  case (octupole$) 

    k3l = ele%value(k3$) * length / (1 + c00%z%vel)
    call mat4_multipole (k3l/2, 0.0_rdef, 3, c00%vec, kmat1)
    mat4 = kmat1
    mat4(1,1:4) = kmat1(1,1:4) + length * kmat1(2,1:4) ! kick * length
    mat4(3,1:4) = kmat1(3,1:4) + length * kmat1(4,1:4)
    call mat4_multipole (k3l/2, 0.0_rdef, 3, c11%vec, kmat2)
    mat6(1:4,1:4) = matmul (kmat2, mat4)

    c00%vec(1:4) = matmul(kmat1, c00%vec(1:4))

    mat6(1,6) = -(kmat2(1,1)*c00%vec(2) + kmat2(1,3)*c00%vec(4)) * length
    mat6(2,6) = -(kmat2(2,1)*c00%vec(2) + kmat2(2,3)*c00%vec(4)) * length
    mat6(3,6) = -(kmat2(3,1)*c00%vec(2) + kmat2(3,3)*c00%vec(4)) * length
    mat6(4,6) = -(kmat2(4,1)*c00%vec(2) + kmat2(4,3)*c00%vec(4)) * length

    mat6(5,1) = -(c00%vec(2)*kmat2(2,1) + c00%vec(4)*kmat2(4,1)) * length
    mat6(5,2) = -(c00%vec(2)*kmat2(2,2) + c00%vec(4)*kmat2(4,2)) * length
    mat6(5,3) = -(c00%vec(2)*kmat2(2,3) + c00%vec(4)*kmat2(4,3)) * length
    mat6(5,4) = -(c00%vec(2)*kmat2(2,4) + c00%vec(4)*kmat2(4,4)) * length

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! solenoid

  case (solenoid$) 

    ks = ele%value(ks$) / (1 + c00%z%vel)

    call solenoid_mat_calc (ks, length, mat6(1:4,1:4))

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

    mat6(5,1) = 2 * c00%vec(1) * t5_11 + c00%vec(4) * t5_14
    mat6(5,2) = 2 * c00%vec(2) * t5_22 + c00%vec(3) * t5_23
    mat6(5,3) = 2 * c00%vec(3) * t5_33 + c00%vec(2) * t5_23
    mat6(5,4) = 2 * c00%vec(4) * t5_44 + c00%vec(1) * t5_14

    mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2) + &
                    mat6(5,4) * mat6(1,3) - mat6(5,3) * mat6(1,4)
    mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2) + &
                    mat6(5,4) * mat6(2,3) - mat6(5,3) * mat6(2,4)
    mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4) + &
                    mat6(5,2) * mat6(3,1) - mat6(5,1) * mat6(3,2) 
    mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4) + &
                    mat6(5,2) * mat6(4,1) - mat6(5,1) * mat6(4,2)

!--------------------------------------------------------
! rf cavity
! this is not quite correct.

  case (rfcavity$) 
    if (ele%value(volt$) /= 0) then
      if (ele%value(harmon$) == 0) then
        type *, 'ERROR IN MAKE_MAT6_BMAD: "HARMON" ATTRIBUTE NOT SET FOR RF.'
        type *, '      FOR ELEMENT: ', ele%name
        call err_exit
      else
        mat6(6,5) = ele%value(volt$) * cos(twopi*ele%value(lag$)) *  &
                   twopi / ele%value(rf_wavelength$) /param%energy / 1.e9
      endif
    endif

    c00%vec = (c00%vec + c11%vec) / 2

    mat6(1,2) = length
    mat6(3,4) = length
    mat6(5,2) = - c00%vec(2) * length
    mat6(5,4) = - c00%vec(4) * length

    mat6(1,5) = mat6(5,2) * mat6(6,5)      
    mat6(1,6) = mat6(5,2)
    mat6(3,5) = mat6(5,4) * mat6(6,5)
    mat6(3,6) = mat6(5,4)

!--------------------------------------------------------
! beam-beam interaction

  case (beambeam$)         

    n_slice = nint(ele%value(n_slice$))
    if (n_slice < 1) then
      type *, 'ERROR IN MAKE_MAT6_BMAD: N_SLICE FOR BEAMBEAM ELEMENT IS NEGATIVE'
      call type_ele (ele, .true., 0, .false., 0, .false.)
      call exit
    endif

    if (ele%value(charge$) == 0 .or. param%n_part == 0) return

! factor of 2 in orb.z.pos since relative motion of the two beams is 2*c_light

    if (n_slice == 1) then
      call bbi_kick_matrix (ele, c00, 0.0_rdef, mat6)
    else
      call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)

      s_pos = 0          ! start at IP
      orb%x%vel = c00%x%vel - ele%value(x_pitch$)
      orb%y%vel = c00%y%vel - ele%value(y_pitch$)
      call mat_make_unit (mat4)

      do i = 1, n_slice + 1
        s_pos_old = s_pos  ! current position
        s_pos = (z_slice(i) + c00%z%pos) / 2 ! position of slice relative to IP
        del_l = s_pos - s_pos_old
        mat4(1,1:4) = mat4(1,1:4) + del_l * mat4(2,1:4)
        mat4(3,1:4) = mat4(3,1:4) + del_l * mat4(4,1:4)
        if (i == n_slice + 1) exit
        orb%x%pos = c00%x%pos + s_pos * orb%x%vel
        orb%y%pos = c00%y%pos + s_pos * orb%y%vel
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

    k1 = ele%value(k1$) / (1 + c00%z%vel)**2
  
! octuple correction to k1

    y_ave = (c00%y%pos + c11%y%pos) / 2
    k_z = pi * ele%value(n_pole$) / length
    k1 = k1 * (1 + 2 * (k_z * y_ave)**2)   

! correction for fact that wigglers with odd number of poles have end
! poles with modified bending radius

    n_pole = nint(ele%value(n_pole$))
    if (mod(n_pole, 2) == 1) then
      rho = 4 * ele%value(rho$) / pi
      l_period = length / n_pole
      l_bend = 8 * l_period / pi**2
      l_drift = l_period - l_bend
      factor = sqrt(rho**2 - (l_bend/2)**2)
      dx = 2 * (rho - factor) + l_drift * l_bend / (2 * factor)
      factor = -2 * rho * dx / (l_bend * (length-l_bend/2))
      k1 = k1 * (1 + factor / n_pole)
    endif

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
             
    mat6(5,2) = 2 * c00%vec(2) * t5_22
    mat6(5,3) = 2 * c00%vec(3) * t5_33 +     c00%vec(4) * t5_34
    mat6(5,4) =     c00%vec(3) * t5_34 + 2 * c00%vec(4) * t5_44

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

    ks = ele%value(ks$) / (1 + c00%z%vel)
    k1 = ele%value(k1$) / (1 + c00%z%vel)

    call sol_quad_mat6_calc (ks, k1, length, mat6, c00%vec)

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! multipole

  case (multipole$, ab_multipole$) 
    if (.not. ele%multipoles_on) return
    call multipole_ele_to_kt (ele, param%particle, knl, tilt, .true.)
    call mat6_multipole (knl, tilt, c00%vec, 1.0_rdef, ele%mat6)

!--------------------------------------------------------
! accelerating solenoid with steerings
! WARNING: This 6x6 matrix may produce bad results at low energies!

  case (accel_sol$) 

    if ((ele%value(s_st1$) < 0.) .or.  &
        (ele%value(s_st1$) + ele%value(l_st1$) > ele%value(s_st2$)) .or.  &
        (ele%value(s_st2$) + ele%value(l_st2$) > ele%value(l$))) then
      type *, 'ERROR IN MAKE_MAT6_BMAD: STEERINGS MUST NOT OVERLAP AND MUST BE',  &
        ' CONTAINED WITHIN'
      type *, 'THE ACCEL_SOL ELEMENT!'
      call type_ele(ele, .true., 0, .false., 0, .false.)
      call exit
    endif

    call mat_make_unit (mat6)     ! make a unit matrix
    if (ele%value(volt$) /= 0) then
      if (ele%value(rf_wavelength$) == 0) then
        type *, 'ERROR IN MAKE_MAT6_BMAD: RF IS ON BUT "RF_WAVELENGTH" NOT SET',  &
              ' IN ACCEL_SOL!'
        call err_exit
      else
        mat6(6,5) = ele%value(volt$) * cos(twopi*ele%value(lag$)) *  &
                      twopi / ele%value(rf_wavelength$) /param%energy / 1.e9
        c_e = ele%value(volt$) * sin(twopi * ele%value(lag$))  &
              / (e_mass * 1.e9 * length)
      endif
    else
      c_e = 0.0
    endif
    c_m = param%particle * c_light * ele%value(b_z$) / (e_mass * 1.e9)
    gamma_old = param%energy * (c00%z%vel + 1) / e_mass
    gamma_new = gamma_old + c_e * length
    call accel_sol_mat_calc (length, c_m, c_e, gamma_old, gamma_new, &
                                    0.0_rdef, 0.0_rdef, c00%vec, mat4, vec_st)
    mat4 = mat6(1:4,1:4)

!--------------------------------------------------------
! rbends are not allowed internally

  case (rbend$) 

    type *, 'ERROR IN MAKE_MAT6_BMAD: RBEND ELEMENTS NOT ALLOWED INTERNALLY!'
    call err_exit

!--------------------------------------------------------
! Custom

  case (custom$)

    type *, 'ERROR IN MAKE_MAT6_BMAD: MAT6_CALC_METHOD = BMAD_STANDARD IS NOT'
    type *, '      ALLOWED FOR A CUSTOM ELEMENT: ', ele%name
    call err_exit

!--------------------------------------------------------
! unrecognized element

  case default

    type *, 'ERROR IN MAKE_MAT6_BMAD: UNKNOWN ELEMENT KEY:', ele%key
    type *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

!--------------------------------------------------------
! put in multipole components

8000 continue

  if (associated(ele%a)) then
    mat6_m = 0
    call multipole_ele_to_kt (ele, param%particle, knl, tilt, .true.)
    call mat6_multipole (knl, tilt, c00%vec, 0.5_rdef, mat6_m)
    mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)
    mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
    mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
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

  use bmad

  implicit none

  type (ele_struct)  ele
  type (coord_struct)  orb

  real(rdef) x_pos, y_pos, del, sig_x, sig_y, coef, garbage, s_pos
  real(rdef) ratio, k0_x, k1_x, k0_y, k1_y, mat6(6,6), beta

!

  sig_x = ele%value(sig_x$)
  sig_y = ele%value(sig_y$)

  if (s_pos /= 0 .and. ele%x%beta /= 0) then
    beta = ele%x%beta - 2 * ele%x%alpha * s_pos + ele%x%gamma * s_pos**2
    sig_x = sig_x * sqrt(beta / ele%x%beta)
    beta = ele%y%beta - 2 * ele%y%alpha * s_pos + ele%y%gamma * s_pos**2
    sig_y = sig_y * sqrt(beta / ele%y%beta)
  endif


  x_pos = orb%x%pos / sig_x  ! this has offset in it
  y_pos = orb%y%pos / sig_y

  del = 0.001

  ratio = sig_y / sig_x
  call bbi_kick (x_pos, y_pos, ratio, k0_x, k0_y)
  call bbi_kick (x_pos+del, y_pos, ratio, k1_x, garbage)
  call bbi_kick (x_pos, y_pos+del, ratio, garbage, k1_y)

  coef = ele%value(bbi_const$) / (1 + orb%z%vel)

  call mat_make_unit (mat6)
  mat6(2,1) = coef * (k1_x - k0_x) / (ele%value(n_slice$) * del * sig_x)
  mat6(4,3) = coef * (k1_y - k0_y) / (ele%value(n_slice$) * del * sig_y)

end subroutine

