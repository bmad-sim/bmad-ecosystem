!+
! Subroutine track1_accel_sol (start, ele, param, end)
!
! Subroutine to track through an accel_sol element
!
! Modules Needed:
!   use bmad
!
! Input:
!   START  -- Coord_struct: Starting position
!   ELE    -- Ele_struct: Element
!   PARAM  -- Param_struct:
!     %APERTURE_LIMIT_ON -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   END   -- Coord_struct: End position
!   PARAM
!     %LOST -- Set .true. If the particle is outside the aperture and
!                %APERTURE_LIMIT_ON is set. Also: %LOST is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %APERTURE_LIMIT_ON.
!-
                
subroutine track1_accel_sol (start, ele, param, end) 

  use bmad

  implicit none

  type (coord_struct)  start, end, c0
  type (ele_struct)  ele
  type (param_struct)  param

  real(rdef) gamma_b, gamma_new, gamma_old, l_over_gamma, ll, ls(5), s_cumul
  real(rdef) s_grand_cum, vec_st(4), x_beg_lim, y_beg_lim, x_lim, y_lim
  real(rdef) x_lim_chng_rate, y_lim_chng_rate, phase, mat4(4,4), length
  real(rdef) b_x(5), b_y(5), beta_b, beta_s, c_e, c_m, en_gain, gam_inv2_b

  real(rdef), parameter :: beta_crit$ = 0.999

  integer i, j

  logical check_x_ap, check_y_ap, inside_segment

! linac stuff
! beta_b is the total speed in units of c_light (before entering the element)
! beta_s is the longitudinal speed in units of c_light

  gamma_b = param%energy * (end%z%vel + 1) / e_mass
  beta_b = sqrt(1 - 1 / gamma_b**2)
  gam_inv2_b = 1.0 / gamma_b**2
  if (gam_inv2_b <= 0.001) then
    beta_s = 1/sqrt(1 + end%x%vel**2 + end%y%vel**2) *(1-gam_inv2_b/2)
  else
    beta_s = sqrt((1 - gam_inv2_b) / (1 + end%x%vel**2 + end%y%vel**2))
  endif

! Calculation of end.z.vel

  end = start
  length = ele%value(l$)

  call offset_particle (ele, param, end, set$)

  if (ele%value(volt$) /= 0) then
    phase = ele%value(lag$) + end%z%pos  &
                      / (ele%value(rf_wavelength$) * beta_s)
    en_gain = ele%value(volt$) * sin(twopi * phase) / 1.e9
    if ((en_gain + gamma_b * e_mass) <= e_mass) then
      param%lost = .true.
      return
    else
      end%z%vel = end%z%vel + en_gain / param%energy
      c_e = en_gain / (e_mass * length)
    endif
  else
    c_e = 0.0
  endif

! Beginning fringe

  c_m = param%particle * c_light * ele%value(b_z$) / (e_mass * 1.e9)
  call mat_unit(mat4, 4, 4)
  mat4(2,3) = c_m / 2 *  &
                  sqrt((1 + end%x%vel**2 + end%y%vel**2) / (gamma_b**2 - 1))
  mat4(4,1) = -mat4(2,3)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

! Segment before first steerings:

  ls(1) = ele%value(s_st1$)
  b_x(1) = 0
  b_y(1) = 0

! Segment with first steerings:

  ls(2) = ele%value(l_st1$)
  b_x(2) = ele%value(b_x1$)
  b_y(2) = ele%value(b_y1$)

! Segment between steerings:

  ls(3) = ele%value(s_st2$) - (ele%value(s_st1$) + ele%value(l_st1$))
  b_x(3) = 0
  b_y(3) = 0

! Segment with second steerings:

  ls(4)= ele%value(l_st2$)
  b_x(4) = ele%value(b_x2$)
  b_y(4) = ele%value(b_y2$)

! Segment after second steerings:

  ls(5) = length - (ele%value(s_st2$) + ele%value(l_st2$))
  b_x(5) = 0
  b_y(5) = 0

  gamma_old = gamma_b
  gamma_new = gamma_old

  s_grand_cum = 0
  x_beg_lim = ele%value(x_beg_limit$)
  y_beg_lim = ele%value(y_beg_limit$)
  x_lim_chng_rate = (ele%value(x_limit$) - x_beg_lim) / length
  y_lim_chng_rate = (ele%value(y_limit$) - y_beg_lim) / length
  check_x_ap = .false.
  check_y_ap = .false.
  if (abs(c_m) > 0.001) then
    l_over_gamma = pi / (10 * abs(c_m))
    if (ele%value(x_limit$) * x_beg_lim /= 0) check_x_ap = .true.
    if (ele%value(y_limit$) * y_beg_lim /= 0) check_y_ap = .true.
  else
    l_over_gamma = 100 * pi
  endif

  do j = 1, 5

    s_cumul = 0
    inside_segment = .true.
    do while (inside_segment)
      gamma_old = gamma_new
      ll = l_over_gamma * gamma_old
      if (s_cumul + ll >= ls(j)) then
        ll = ls(j) - s_cumul
        inside_segment = .false.
      endif
      s_cumul = s_cumul + ll
      s_grand_cum = s_grand_cum + ll

      if (ll > 1.e-5) then
        gamma_new = gamma_old + c_e * ll
        call accel_sol_mat_calc (ll, c_m, c_e, gamma_old, gamma_new,  &
        b_x(j), b_y(j), end, mat4, vec_st)
        end%vec(1:4) = matmul(mat4, end%vec(1:4))
        do i = 1, 4
          end%vec(i) = end%vec(i) + param%particle * vec_st(i)
        enddo

! Calculation of end.z.pos

        if (abs(c_e) > 0.001) then
          end%z%pos = end%z%pos - ll + (sqrt(1 + (c_e * ll  &
                          + sqrt(gamma_old**2 - 1))**2) - gamma_old) / c_e
        else
          end%z%pos = end%z%pos + ll * (sqrt(1 - 1/gamma_old**2) - 1  &
              - c_e * ll * (gamma_old**2 - 0.5)/gamma_old**3)
        endif

        if (param%aperture_limit_on) then
          if (check_x_ap) then
            x_lim = x_lim_chng_rate * s_grand_cum + x_beg_lim
            if (abs(end%x%pos) > x_lim) param%lost = .true.
          endif
          if (check_y_ap) then
            y_lim = y_lim_chng_rate * s_grand_cum + y_beg_lim
            if (abs(end%y%pos) > y_lim) param%lost = .true.
          endif
          if (param%lost) return
        endif

      endif
    enddo

  enddo

! Ending fringe

  call mat_unit(mat4, 4, 4)
  mat4(4,1) = c_m / 2 *  &
                sqrt((1 + end%x%vel**2 + end%y%vel**2) / (gamma_new**2 - 1))
  mat4(2,3) = -mat4(4,1)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

  call offset_particle (ele, param, end, unset$)

end subroutine
