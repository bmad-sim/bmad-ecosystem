!+
! Subroutine MAKE_MAT6 (ELE, PARAM, C0, C1)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     ELE    -- Ele_struct: Element
!     PARAM  -- Param_struct: Parameters are needed for some elements.
!     C0, C1 -- [Optional] Coord_struct: Coordinates at the beginning and 
!               end of the element. If
!               not present then the closed orbit is assumed to be the origin. 
!
! Output:
!     ELE%MAT6 -- Real: 6x6 transfer matrix.
!-

subroutine make_mat6 (ele, param, c0, c1)

  use bmad_struct
  implicit none

  type (ele_struct), target :: ele
  type (coord_struct), optional :: c0, c1
  type (coord_struct) orb, c00, c11
  type (param_struct)  param

  real, pointer :: mat6(:,:)

  real mat6_m(6,6), mat2(2,2), mat4(4,4), kmat1(4,4), kmat2(4,4)
  real e1, e2, angle, rho, cos_angle, sin_angle, k1, ks, length, kc
  real phi, k2l, k3l, c2, s2, cs, ks2, del_l, rho_bend, l_period, l_bend
  real factor, l_drift, dx
  real s_pos, s_pos_old, z_slice(100)
  real knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real r, c_e, c_m, gamma_old, gamma_new, vec_st(4)
  real sqrt_k, arg, kick2
  real cx, sx, cy, sy, k2l_2, k2l_3, k2l_4
  real x_off, y_off, x_pit, y_pit, y_ave, k_z
  real t5_11, t5_12, t5_22, t5_33, t5_34, t5_44, t5_14, t5_23
  real t1_16, t1_26, t1_36, t1_46, t2_16, t2_26, t2_36, t2_46
  real t3_16, t3_26, t3_36, t3_46, t4_16, t4_26, t4_36, t4_46
  real lcs, lc2s2

  integer i, n, n_slice, n_pole

  logical unit_multipole_matrix

!--------------------------------------------------------
! init

  length = ele%value(l$)
  mat6 => ele%mat6
  mat6 = 0
  forall (i = 1:6) mat6(i,i) = 1

  ele%coupled = .true.  ! initially assume x-y coupling

  if (.not. present(c0)) then
    c00%vec = 0
    c11%vec = 0
  else
    c00 = c0
    c11 = c1
  endif

  if (ele%is_on .and. ele%multipoles_on .and. &
                                     any(ele%value(ix1_m$:ix2_m$) /= 0)) then
    ele%nonzero_multipoles = .true.
  else
    ele%nonzero_multipoles = .false.
  endif

!--------------------------------------------------------
! marker

  if (ele%key == marker$) then
    ele%coupled = .false.
    return
  endif

!--------------------------------------------------------
! drift or element is off or
! Electric Separator or Kicker.

  if (ele%key == drift$ .or. (.not. ele%is_on) .or. &
      ele%key == elseparator$ .or. ele%key == kicker$) then

    c00%vec = (c00%vec + c11%vec) / 2

    mat6(1,2) = length
    mat6(3,4) = length
    mat6(1,6) =        - c00%vec(2) * length
    mat6(3,6) =        - c00%vec(4) * length
    mat6(5,2) =        - c00%vec(2) * length
    mat6(5,4) =        - c00%vec(4) * length

    ele%coupled = .false.
    goto 8000   ! put in multipole ends if needed

  endif

!--------------------------------------------------------
! Put in offsets, etc.
! Note: c00 and c11 are the coords in the frame of reference where the element
!       is upright (no tilt).

  x_off = ele%value(x_offset$); y_off = ele%value(y_offset$)
  x_pit = ele%value(x_pitch$) * length / 2
  y_pit = ele%value(y_pitch$) * length / 2
  c00%vec = c00%vec - (/ x_off - x_pit*length/2, x_pit, &
                           y_off - y_pit*length/2, y_pit, 0.0, 0.0 /) 
  c11%vec = c11%vec - (/ x_off + x_pit*length/2, x_pit, &
                           y_off + y_pit*length/2, y_pit, 0.0, 0.0 /) 

  if (ele%value(tilt$) /= 0) then
    if (ele%key /= multipole$ .and. ele%key /= ab_multipole$) then
      call tilt_coords (ele%value(tilt$), c00%vec, .true.)
      call tilt_coords (ele%value(tilt$), c11%vec, .true.)
    endif
  endif

!--------------------------------------------------------
! sbend

  select case (ele%key)

  case (sbend$) 

    angle = ele%value(angle$)
    e1 = ele%value(e1$)
    e2 = ele%value(e2$)
    k1 = ele%value(k1$)
    if (length == 0) then
      rho = 1  ! so subroutine will not bomb
    else
      rho = length / angle
    endif
    ele%value(rho$) = rho

    e1 = e1 + c00%x%vel
    e2 = e2 - c11%x%vel
    angle = angle + c00%x%vel - c11%x%vel
    rho = rho * (1 + c00%z%vel)
    length = rho * angle
    k1 = k1 / (1 + c00%z%vel)
   
    kc = 1/rho**2 + k1
  
    call quad_mat_calc (-kc, length, mat6(1:2,1:2))
    call quad_mat_calc (k1, length, mat6(3:4,3:4))

    phi = sqrt(abs(kc)) * length
    if (kc < 0.0) then
      mat6(1,6) = (1 - cosh(phi)) / (rho * kc)
      mat6(2,6) = sinh(phi) / (rho * sqrt(-kc))
      mat6(5,6) = (phi - sinh(phi)) / (rho**2 * abs(kc)**1.5)
    else
      mat6(1,6) = (1 - cos(phi)) / (rho * kc)
      mat6(2,6) = sin(phi) / (rho * sqrt(kc))
      mat6(5,6) = (sin(phi) - phi) / (rho**2 * kc**1.5)
    endif

    mat6(5,1) = -mat6(2,6)
    mat6(5,2) = -mat6(1,6)

    mat6(5,4) = -c00%vec(4) * (mat6(5,6) + length)
    mat6(3,6) = -c00%vec(4) * (mat6(5,6) + length)

    if (e1 /= 0) then
      arg = tan(e1) / rho
      mat6(1:6,1) = mat6(1:6,1) + mat6(1:6,2) * arg
      mat6(1:6,3) = mat6(1:6,3) - mat6(1:6,4) * arg
    endif

    if (e2 /= 0) then
      arg = tan(e2) / rho
      mat6(2,1:6) = mat6(2,1:6) + mat6(1,1:6) * arg
      mat6(4,1:6) = mat6(4,1:6) - mat6(3,1:6) * arg
    endif

    if (ele%value(tilt$)+ele%value(roll$) == 0) then
      ele%coupled = .false.
    else
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

    if (ele%value(tilt$) == 0) then
      ele%coupled = .false.
    else
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! Sextupole.
! the sextupole is modeled as kick-drift-kick

  case (sextupole$) 
   
    k2l = ele%value(k2$) * length / (1 + c00%z%vel)
    call mat4_multipole (ele, k2l/2, 0.0, 2, c00, kmat1)
    mat4 = kmat1
    mat4(1,1:4) = kmat1(1,1:4) + length * kmat1(2,1:4) ! kick * length
    mat4(3,1:4) = kmat1(3,1:4) + length * kmat1(4,1:4)
    call mat4_multipole (ele, k2l/2, 0.0, 2, c11, kmat2)
    mat6(1:4,1:4) = matmul (kmat1, mat4)

    c00%vec(1:4) = matmul(kmat1, c00%vec(1:4))

    mat6(1,6) = -(kmat2(1,1)*c00%vec(2) + kmat2(1,3)*c00%vec(4)) * length
    mat6(2,6) = -(kmat2(2,1)*c00%vec(2) + kmat2(2,3)*c00%vec(4)) * length
    mat6(3,6) = -(kmat2(3,1)*c00%vec(2) + kmat2(3,3)*c00%vec(4)) * length
    mat6(4,6) = -(kmat2(4,1)*c00%vec(2) + kmat2(4,3)*c00%vec(4)) * length

    mat6(5,1) = -(c00%vec(2)*kmat2(2,1) + c00%vec(4)*kmat2(4,1)) * length
    mat6(5,2) = -(c00%vec(2)*kmat2(2,2) + c00%vec(4)*kmat2(4,2)) * length
    mat6(5,3) = -(c00%vec(2)*kmat2(2,3) + c00%vec(4)*kmat2(4,3)) * length
    mat6(5,4) = -(c00%vec(2)*kmat2(2,4) + c00%vec(4)*kmat2(4,4)) * length

    if (k2l /= 0) ele%coupled = .true.

    if (ele%value(tilt$) /= 0) then
      call tilt_mat6 (mat6, ele%value(tilt$))
    endif

!--------------------------------------------------------
! octupole
! the octupole is modeled as kick-drift-kick

  case (octupole$) 

    k3l = ele%value(k3$) * length / (1 + c00%z%vel)
    call mat4_multipole (ele, k3l/2, 0.0, 3, c00, kmat1)
    mat4 = kmat1
    mat4(1,1:4) = kmat1(1,1:4) + length * kmat1(2,1:4) ! kick * length
    mat4(3,1:4) = kmat1(3,1:4) + length * kmat1(4,1:4)
    call mat4_multipole (ele, k3l/2, 0.0, 3, c11, kmat2)
    mat6(1:4,1:4) = matmul (kmat1, mat4)

    c00%vec(1:4) = matmul(kmat1, c00%vec(1:4))

    mat6(1,6) = -(kmat2(1,1)*c00%vec(2) + kmat2(1,3)*c00%vec(4)) * length
    mat6(2,6) = -(kmat2(2,1)*c00%vec(2) + kmat2(2,3)*c00%vec(4)) * length
    mat6(3,6) = -(kmat2(3,1)*c00%vec(2) + kmat2(3,3)*c00%vec(4)) * length
    mat6(4,6) = -(kmat2(4,1)*c00%vec(2) + kmat2(4,3)*c00%vec(4)) * length

    mat6(5,1) = -(c00%vec(2)*kmat2(2,1) + c00%vec(4)*kmat2(4,1)) * length
    mat6(5,2) = -(c00%vec(2)*kmat2(2,2) + c00%vec(4)*kmat2(4,2)) * length
    mat6(5,3) = -(c00%vec(2)*kmat2(2,3) + c00%vec(4)*kmat2(4,3)) * length
    mat6(5,4) = -(c00%vec(2)*kmat2(2,4) + c00%vec(4)*kmat2(4,4)) * length

    if (k3l /= 0) ele%coupled = .true.

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
    if (ele%value(harmon$) /= 0) ele%value(rf_wavelength$) =  &
                                   param%total_length / ele%value(harmon$)
    if (ele%value(volt$) /= 0) then
      if (ele%value(harmon$) == 0) then
        type *, 'ERROR IN MAKE_MAT6: "HARMON" ATTRIBUTE NOT SET FOR RF.'
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

    ele%coupled = .false.

!--------------------------------------------------------
! beam-beam interaction

  case (beambeam$)         

    n_slice = nint(ele%value(n_slice$))
    if (n_slice == 0) then
      ele%value(n_slice$) = 1.0  ! revert to default
      n_slice = 1
    elseif (n_slice < 1) then
      type *, 'ERROR IN MAKE_MAT6: N_SLICE FOR BEAMBEAM ELEMENT IS NEGATIVE'
      call type_ele (ele, .true., 0, .false., .false.)
      call exit
    endif

    if (ele%value(charge$) == 0 .or. param%n_part == 0) then
      ele%coupled = .false.
      ele%value(bbi_const$) = 0
      return
    endif

    if (ele%value(sig_x$) == 0 .or. ele%value(sig_y$) == 0) then
      type *, 'ERROR IN MAKE_MAT6: ZERO SIGMA IN BEAMBEAM ELEMENT!'
      call type_ele(ele, .true., 0, .false., .false.)
      call exit
    endif

    ele%value(bbi_const$) = -param%n_part * e_mass * ele%value(charge$) * &
      r_e /  (2 * pi * param%energy * (ele%value(sig_x$) + ele%value(sig_y$)))

! factor of 2 in orb.z.pos since relative motion of the two beams is 2*c_light

    if (n_slice == 1) then
      call bbi_kick_matrix (ele, c00, 0.0, mat6)
    else
      call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)

      s_pos = 0          ! start at IP
      orb%x%vel = c00%x%vel - ele%value(x_pitch$)
      orb%y%vel = c00%y%vel - ele%value(y_pitch$)
      call mat_unit (mat4, 4, 4)

      do i = 1, n_slice + 1
        s_pos_old = s_pos  ! current position
        s_pos = (z_slice(i) + c00%z%pos) / 2  ! position of slice relative to IP
        del_l = s_pos - s_pos_old
        mat4(1,1:4) = mat4(1,1:4) + del_l * mat4(2,1:4)
        mat4(3,1:4) = mat4(3,1:4) + del_l * mat4(4,1:4)
        if (i == n_slice + 1) exit
        orb%x%pos = c00%x%pos + s_pos * orb%x%vel
        orb%y%pos = c00%y%pos + s_pos * orb%y%vel
        call bbi_kick_matrix (ele, orb, s_pos, kmat1)
        mat4(2,1:4) = mat4(2,1:4) + kmat1(2,1) * mat4(1,1:4) + &
                                    kmat1(2,3) * mat4(3,1:4)
        mat4(4,1:4) = mat4(4,1:4) + kmat1(4,1) * mat4(1,1:4) + &
                                    kmat1(4,3) * mat4(3,1:4)
      enddo

      mat6(1:4,1:4) = mat4

    endif

!--------------------------------------------------------
! wiggler

  case (wiggler$) 

    call mat_unit (mat6, 6, 6)     ! make a unit matrix

    if (param%energy == 0) then
      ele%value(k1$) = 0
    else
      ele%value(k1$) = -0.5 * (0.2997 * ele%value(b_max$) / param%energy)**2
    endif

    if (ele%value(b_max$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = 3.3356 * param%energy / ele%value(b_max$)
    endif

    k1 = ele%value(k1$) / (1 + c00%z%vel)**2
  
! octuple correction to k1

    y_ave = (c00%y%pos + c11%y%pos) / 2
    k_z = 2 * pi * ele%value(n_pole$) / length
    k1 = k1 * (1 + 2 * (k_z * y_ave)**2)   

! correction for fact that wigglers with odd number of poles have end
! poles with modified bending radius

    n_pole = nint(ele%value(n_pole$))
    if (mod(n_pole, 2) == 1) then
      rho_bend = 4 * ele%value(rho$) / pi
      l_period = length / n_pole
      l_bend = 8 * l_period / pi**2
      l_drift = l_period - l_bend
      factor = sqrt(rho_bend**2 - (l_bend/2)**2)
      dx = 2 * (rho_bend - factor) + l_drift * l_bend / (2 * factor)
      factor = -2 * rho_bend * dx / (l_bend * (length-l_bend/2))
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

    if (ele%value(tilt$) == 0) then
      ele%coupled = .false.
    else
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
    call mat6_multipole (ele, param, c00, 1.0, ele%mat6, unit_multipole_matrix)
    if (unit_multipole_matrix) ele%coupled = .false.
    return

!--------------------------------------------------------
! loop

  case (loop$) 

! make a zero matrix because no transfer matrix may be calculated
! ahead of time

    mat6 = 0   
    r = ele%value(radius$)
    ele%value(diameter$) = 2 * r
    ele%value(r2$) = r**2
    ele%value(ri$) = r * ele%value(current$)
    ele%value(r2i$) = ele%value(r2$) * ele%value(current$)

!--------------------------------------------------------
! coil

  case (coil$)

! make a zero matrix because no transfer matrix may be calculated
! ahead of time

    mat6 = 0   

!--------------------------------------------------------
! accelerating solenoid with steerings
! WARNING: This 6x6 matrix may produce bad results at low energies!

  case (accel_sol$) 

    if ((ele%value(s_st1$) < 0.) .or.  &
        (ele%value(s_st1$) + ele%value(l_st1$) > ele%value(s_st2$)) .or.  &
        (ele%value(s_st2$) + ele%value(l_st2$) > ele%value(l$))) then
      type *, 'ERROR IN MAKE_MAT6: STEERINGS MUST NOT OVERLAP AND MUST BE',  &
        ' CONTAINED WITHIN'
      type *, 'THE ACCEL_SOL ELEMENT!'
      call type_ele(ele, .true., 0, .false., .false.)
      call exit
    endif

    call mat_unit (mat6, 6, 6)     ! make a unit matrix
    if (ele%value(volt$) /= 0) then
      if (ele%value(rf_wavelength$) == 0) then
        type *, 'ERROR IN MAKE_MAT6: RF IS ON BUT "RF_WAVELENGTH" NOT SET',  &
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
    call accel_sol_mat_calc (length, c_m, c_e, gamma_old, gamma_new, 0, 0,  &
                                            c00, mat4, vec_st)
    mat4 = mat6(1:4,1:4)

!--------------------------------------------------------
! redefinition of energy

  case (define_energy$) 

    call mat_unit (mat6, 6, 6)     ! make a unit matrix
    param%energy = ele%value(new_energy$)
    ele%coupled = .false.

!--------------------------------------------------------
! rbends are not allowed internally

  case (rbend$) 

    type *, 'ERROR IN MAKE_MAT6: RBEND ELEMENTS NOT ALLOWED INTERNALLY!'
    call err_exit

!--------------------------------------------------------
! unrecognized element

  case (custom$)

    call custom_make_mat6 (ele, param, c0, c1)
    return

!--------------------------------------------------------
! unrecognized element

  case default

    type *, 'ERROR IN MAKE_MAT6: UNKNOWN ELEMENT KEY',  &
                                    ele%key, '  ', key_name(ele%key)

  end select

!--------------------------------------------------------
! put in multipole components

8000 continue

  if (ele%nonzero_multipoles) then
    mat6_m = 0
    call mat6_multipole (ele, param, c00, 0.5, mat6_m, unit_multipole_matrix)
    mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)
    mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
    mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
  endif

end subroutine
                                                                   
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine mat6_multipole (ele, param, c00, factor, mat6, unit_matrix)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct) ele
  type (param_struct) param
  type (coord_struct) c00
                           
  real mat6(6,6), kmat1(4,4), factor
  real knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

  integer n

  logical unit_matrix

!                        

  call multipole_to_vecs (ele, param%particle, knl, tilt)

  if (c00%x%pos == 0 .and. c00%y%pos == 0 .and. knl(1) == 0) then
    unit_matrix = .true.
    return                  
  endif

  do n = 1, n_pole_maxx
    if (knl(n) /= 0) then
      unit_matrix = .false.
      call mat4_multipole (ele, knl(n), tilt(n), n, c00, kmat1)
      mat6(2:4:2, 1:3:2) = mat6(2:4:2, 1:3:2) + factor * kmat1(2:4:2, 1:3:2)
    endif
  enddo

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine QUAD_MAT_CALC (K1, LENGTH, MAT)
!
! Subroutine to initialize the transfer matrix for a quad
!-


subroutine quad_mat_calc (k1, length, mat)
    
  implicit none

  real length, mat(2,2), cx, sx
  real k1, sqrt_k, arg, arg2

!

  sqrt_k = sqrt(abs(k1))
  arg = sqrt_k * length

  if (arg < 1e-3) then
    arg2 = k1 * length**2
    cx = 1 - arg2 / 2
    sx = (1 - arg2 / 6) * length
  elseif (k1 < 0) then       ! focus
    cx = cos(arg)
    sx = sin(arg) / sqrt_k
  else                           ! defocus
    cx = cosh(arg)
    sx = sinh(arg) / sqrt_k
  endif

  mat(1,1) = cx
  mat(1,2) = sx
  mat(2,1) = k1 * sx
  mat(2,2) = cx

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine SOL_QUAD_MAT6_CALC (KS, K1, LENGTH, MAT6, ORB)
!
! Subroutine to calculate the transfer matrix for a combination 
! solenoid/quadrupole element.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ks      [Real]       Solenoid strength
!   k1      [Real]       Quadrupole strength
!   length  [Real]       Sol_quad length
!   orb(6)  [Real]       Orbit at beginning of the sol_quad.
!
! Output
!   mat6(6,6) [Real]  Transfer matrix across the sol_quad
!-

subroutine sol_quad_mat6_calc (ks, k1, s_len, m, orb)

  use bmad_struct

  implicit none

  real ks, k1, s_len
  real m(6,6)
  real orb(6)

  integer i, j
  integer order

  real ks2, s, c, snh, csh
  real darg1, alpha, alpha2, beta, beta2, f, q, r, a, b
  real df, dalpha2, dalpha, dbeta2, dbeta, darg
  real dC, dCsh, dS, dSnh, dq, dr, da, db
  real ks3, fp, fm, dfm, dfp, df_f, ug
  real s1, s2, snh1, snh2, dsnh1, dsnh2, ds1, ds2
  real coef1, coef2, dcoef1, dcoef2, ks4

  real t5(4,4), t6(4,4)

! Calc
          
  ks2 = ks**2
  ks3 = ks2 * ks 
  ks4 = ks2*ks2
  f = sqrt(ks4 + 4*k1**2)
  ug = 1 / (4*f)
  alpha2 = (f + ks2) / 2; alpha = sqrt(alpha2)
  beta2  = (f - ks2) / 2; beta  = sqrt(beta2)
  S = sin(alpha*s_len)                              
  C = cos(alpha*s_len)
  Snh = sinh(beta*s_len)
  Csh = cosh(beta*s_len)
  q = f + 2*k1 - ks2
  r = f - 2*k1 + ks2
  a = f + 2*k1 + ks2
  b = f - 2*k1 - ks2
  fp = f + 2*k1
  fm = f - 2*k1

  S1 = S * alpha
  S2 = S / alpha

  Snh1 = Snh * beta
  Snh2 = Snh / beta

  coef1 = ks2*r + 4*k1*a
  coef2 = ks2*q + 4*k1*b

  call mat_unit(m, 6, 6)
               
  m(1,1) = 2*ug * (fp*C + fm*Csh)
  m(1,2) = (2*ug/k1) * (q*S1 - r*Snh1)
  m(1,3) = (ks*ug/k1) * (-b*S1 + a*Snh1)
  m(1,4) = 4*ug*ks * (-C + Csh)

  m(2,1) = -(ug/2) * (coef1*S2 + coef2*Snh2)
  m(2,2) = m(1,1)             
  m(2,3) = ug*ks3 * (C - Csh)
  m(2,4) = ug*ks * (a*S2 + b*Snh2)

  m(3,1) = -m(2,4)
  m(3,2) = -m(1,4)
  m(3,3) = 2*ug * (fm*C + fp*Csh)  
  m(3,4) = 2*ug * (r*S2 + q*Snh2)

  m(4,1) = -m(2,3)     
  m(4,2) = -m(1,3)
  m(4,3) = (ug/(2*k1)) * (-coef2*S1 + coef1*Snh1)
  m(4,4) = m(3,3)

!

  if (all(orb(1:4) == 0)) return

  df      = -2 * (ks4 + 2*k1**2) / f
  dalpha2 = df/2 - ks2
  dalpha  = (df/2 - ks2)/(2*alpha)
  dbeta2  = ks2 + df/2
  dbeta   = (ks2 + df/2)/(2*beta)
  darg    = s_len*dalpha
  darg1   = s_len*dbeta         
  dC      = -darg*S
  dCsh    = darg1*Snh
  dS      = darg*C
  dSnh    = darg1*Csh
  dq      = -2*k1 + 2*ks2 + df
  dr      =  2*k1 - 2*ks2 + df
  da      = -2*k1 - 2*ks2 + df
  db      =  2*k1 + 2*ks2 + df
  dfp = df - 2*k1
  dfm = df + 2*k1
  df_f =  -df/f

  dS1 = dS * alpha + S * dalpha
  dS2 = dS / alpha - S * dalpha / alpha2

  dSnh1 = dSnh * beta + Snh * dbeta
  dSnh2 = dSnh / beta - Snh * dbeta / beta2

  dcoef1 = -2*ks2*r + ks2*dr - 4*k1*a + 4*k1*da
  dcoef2 = -2*ks2*q + ks2*dq - 4*k1*b + 4*k1*db                     

  t6(1,1) = m(1,1)*df_f + 2*ug*(fp*dC + C*dfp + fm*dCsh + Csh*dfm)
  t6(1,2) = m(1,2)*df_f + (2*ug/k1) * (dq*S1 + q*dS1 - dr*Snh1 - r*dSnh1)
  t6(1,3) = m(1,3)*df_f + (ks*ug/k1)*(-db*S1 - b*dS1 + da*Snh1 + a*dSnh1)
  t6(1,4) = m(1,4)*(df_f - 2) + 4*ks*ug*(-dC + dCsh) 

  t6(2,1) = m(2,1)*(df_f + 1) - &
              (ug/2)*(dcoef1*S2 + coef1*dS2 + dcoef2*Snh2 + coef2*dSnh2)
  t6(2,2) = t6(1,1)
  t6(2,3) = m(2,3)*(df_f - 2) + ks3*ug*(dC - dCsh) 
  t6(2,4) = m(2,4)*(df_f - 1) + ug*ks*(da*S2 + a*dS2 + db*Snh2 + b*dSnh2)

  t6(3,1) = -t6(2,4)
  t6(3,2) = -t6(1,4)
  t6(3,3) = m(3,3)*df_f + 2*ug*(fm*dC + C*dfm + fp*dCsh + Csh*dfp)
  t6(3,4) = m(3,4)*(df_f - 1) + 2*ug*(dr*S2 + r*dS2 + dq*Snh2 + q*dSnh2)

  t6(4,1) = -t6(2,3)        
  t6(4,2) = -t6(1,3)
  t6(4,3) = m(4,3)*(df_f + 2) + &
               (ug/(2*k1))*(-dcoef2*S1 - coef2*dS1 + dcoef1*Snh1 + coef1*dSnh1)
  t6(4,4) = t6(3,3)

!

  m(1:4,6) = matmul(t6(1:4,1:4), orb(1:4))
  m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(1,2) - m(4,6)*m(1,3) + m(3,6)*m(1,4)
  m(5,2) = -m(2,6)*m(2,1) + m(1,6)*m(2,2) - m(4,6)*m(2,3) + m(3,6)*m(2,4)
  m(5,3) = -m(2,6)*m(3,1) + m(1,6)*m(3,2) - m(4,6)*m(3,3) + m(3,6)*m(3,4)
  m(5,4) = -m(2,6)*m(4,1) + m(1,6)*m(4,2) - m(4,6)*m(4,3) + m(3,6)*m(4,4)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! MAT4_MULTIPOLE (ELE, KNL, TILT, N, C0, KICK_MAT)
!
! Subroutine to find the kick from a multipole
!
! Input:
!     C0   -- Coord_struct: coordinates of particle
!     KNL  -- Real: Strength of multipole
!     TILT -- Real: Tilt of multipole
!
! Output:
!     KICK_MAT(4,4) -- Real: Kick matrix
!-


subroutine mat4_multipole (ele, knl, tilt, n, c0, kick_mat)
                  
  use bmad_struct
  implicit none
  type (ele_struct)  ele
  type (coord_struct)  c0

  real x_pos, y_pos, x, y, knl, tilt, c(0:n_pole_maxx, 0:n_pole_maxx)
  real sin_ang, cos_ang, mexp, mat(2,2), rot(2,2)
  real kick_mat(4,4)

  integer m, n

  logical init_needed / .true. /

! init

  if (init_needed) then
    call multipole_c_init (c, n_pole_maxx)
    init_needed = .false.
  endif

  kick_mat = 0
  forall (m = 1:4) kick_mat(m,m) = 1

! simple case

  if (knl == 0 .or. (c0%x%pos == 0 .and. c0%y%pos == 0 .and. n /= 1)) then
    kick_mat(2:4:2, 1:3:2) = 0
    return
  endif

! get position of particle in frame of multipole

  x_pos = c0%x%pos
  y_pos = c0%y%pos
           
  if (tilt == 0) then
    x = x_pos
    y = y_pos
  else
    sin_ang = sin(tilt)
    cos_ang = cos(tilt)
    x =  x_pos * cos_ang + y_pos * sin_ang
    y = -x_pos * sin_ang + y_pos * cos_ang
  endif

! compute kick matrix

  mat = 0

  do m = 0, n, 2
    mat(1,1) = mat(1,1) +  &
                    knl * (n-m) * c(n, m) * mexp(x, n-m-1) * mexp(y, m)
    mat(1,2) = mat(1,2) +  &
                    knl * m * c(n, m) * mexp(x, n-m) * mexp (y, m-1)
  enddo

  do m = 1, n, 2
    mat(2,1) = mat(2,1) +  &
                    knl * (n-m) * c(n, m) * mexp(x, n-m-1) * mexp(y, m)
    mat(2,2) = mat(2,2) +  &
                    knl * m * c(n, m) * mexp(x, n-m) * mexp(y, m-1)
  enddo

! transform back to lab frame

  if (tilt /= 0) then
    rot(1,1) =  cos_ang
    rot(1,2) = -sin_ang
    rot(2,1) =  sin_ang
    rot(2,2) =  cos_ang
    mat = matmul(rot, mat)
    rot(1,2) =  sin_ang
    rot(2,1) = -sin_ang
    mat = matmul (mat, rot)
  endif

  kick_mat(2,1) = mat(1,1)
  kick_mat(2,3) = mat(1,2)
  kick_mat(4,1) = mat(2,1)
  kick_mat(4,3) = mat(2,2)

  return
  end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


  function mexp (x, m)

  real x, mexp
  integer m

!

  if (m < 0) then
    mexp = 0
  elseif (m == 0) then
    mexp = 1
  else
    mexp = x**m
  endif

  return
  end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


subroutine bbi_kick_matrix (ele, orb, s_pos, mat6)

  use bmad_struct
  implicit none

  type (ele_struct)  ele
  type (coord_struct)  orb

  real x_pos, y_pos, del, sig_x, sig_y, coef, garbage, s_pos
  real ratio, k0_x, k1_x, k0_y, k1_y, mat6(6,6), beta

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

  call mat_unit(mat6, 6, 6)
  mat6(2,1) = coef * (k1_x - k0_x) / (ele%value(n_slice$) * del * sig_x)
  mat6(4,3) = coef * (k1_y - k0_y) / (ele%value(n_slice$) * del * sig_y)

  ele%coupled = .false.

  return

  end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


subroutine bbi_slice_calc (n_slice, sig_z, z_slice)

  implicit none

  integer i, n_slice, n_slice_old / 0 /

  real sig_z, z_slice(*), y, inverse, z_norm(100)

  external probability_funct

  save z_norm

!

  if (n_slice == 1) then
    z_slice(1) = 0
  elseif (n_slice > 1) then
    do i = 1, n_slice
      if (n_slice /= n_slice_old) then
        y = (i - 0.5) / n_slice - 0.5
        z_norm(i) = inverse(probability_funct, y, -5.0, 5.0, 1.0e-5)
      endif
      z_slice(i) = sig_z * z_norm(i)
    enddo
    n_slice_old = n_slice
  else
    type *, 'ERROR IN BBI_SLICE_CALC: N_SLICE IS NEGATIVE:', n_slice
    call err_exit
  endif

  z_slice(n_slice+1) = 0

  end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+      
! Subroutine tilt_mat6 (mat6, tilt)
!
! Subroutine to tilt a 6x6 matrix.
!-

subroutine tilt_mat6 (mat6, tilt)

  use bmad_struct

  implicit none

  real tilt, mat6(6,6), mm(6,6)
  real c, s, c2, cs, s2

!

  if (tilt == 0) return

  c = cos(tilt)
  s = sin(tilt)

  mm(1,1:6) = c * mat6(1,1:6) - s * mat6(3,1:6)
  mm(2,1:6) = c * mat6(2,1:6) - s * mat6(4,1:6)
  mm(3,1:6) = c * mat6(3,1:6) + s * mat6(1,1:6)
  mm(4,1:6) = c * mat6(4,1:6) + s * mat6(2,1:6)
  mm(5,1:6) =     mat6(5,1:6)
  mm(6,1:6) =     mat6(6,1:6)

  mat6(1:6,1) = c * mm(1:6,1) - s * mm(1:6,3)
  mat6(1:6,2) = c * mm(1:6,2) - s * mm(1:6,4)
  mat6(1:6,3) = c * mm(1:6,3) + s * mm(1:6,1)
  mat6(1:6,4) = c * mm(1:6,4) + s * mm(1:6,2)
  mat6(1:6,5) =     mm(1:6,5)
  mat6(1:6,6) =     mm(1:6,6)
                     
end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine solenoid_mat_calc (ks, length, mat4)

  implicit none

  real ks, length, kss, c, s, c2, s2, cs
  real mat4(4,4)

  kss = ks / 2

  c = cos(kss*length)
  s = sin(kss*length)
  c2 = c*c
  s2 = s*s
  cs = c*s

  mat4(1,1) = c2
  mat4(1,2) = cs / kss
  mat4(1,3) = cs
  mat4(1,4) = s2 / kss
  mat4(2,1) = -kss * cs
  mat4(2,2) = c2
  mat4(2,3) = -kss * s2 
  mat4(2,4) = cs
  mat4(3,1) = -cs
  mat4(3,2) = -s2 / kss
  mat4(3,3) = c2
  mat4(3,4) = cs / kss
  mat4(4,1) = kss * s2
  mat4(4,2) = -cs
  mat4(4,3) = -kss * cs
  mat4(4,4) = c2

end subroutine
