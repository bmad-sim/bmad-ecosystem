!+
! Module accel_sol
!
! This module contains antequated accel_sol routines.
!-

#include "CESR_platform.inc"

module accel_sol

  use bmad_struct
  use bmad_interface

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_accel_sol (start, ele, param, end)
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
                
subroutine track_a_accel_sol (start, ele, param, end) 

  implicit none

  type (coord_struct)  start, end, c0
  type (ele_struct)  ele
  type (param_struct)  param

  real(rp) gamma_b, gamma_new, gamma_old, l_over_gamma, ll, ls(5), s_cumul
  real(rp) s_grand_cum, vec_st(4), x_beg_lim, y_beg_lim, x_lim, y_lim
  real(rp) x_lim_chng_rate, y_lim_chng_rate, phase, mat4(4,4), length
  real(rp) b_x(5), b_y(5), beta_b, beta_s, c_e, c_m, en_gain, gam_inv2_b

  real(rp), parameter :: beta_crit$ = 0.999

  integer i, j

  logical check_x_ap, check_y_ap, inside_segment

! linac stuff
! beta_b is the total speed in units of c_light (before entering the element)
! beta_s is the longitudinal speed in units of c_light

  gamma_b = ele%value(beam_energy$) * (end%vec(6) + 1) / m_electron
  beta_b = sqrt(1 - 1 / gamma_b**2)
  gam_inv2_b = 1.0 / gamma_b**2
  if (gam_inv2_b <= 0.001) then
    beta_s = 1/sqrt(1 + end%vec(2)**2 + end%vec(4)**2) *(1-gam_inv2_b/2)
  else
    beta_s = sqrt((1 - gam_inv2_b) / (1 + end%vec(2)**2 + end%vec(4)**2))
  endif

! Calculation of end.z.vel

  end = start
  length = ele%value(l$)

  call offset_particle (ele, param, end, set$)

  if (ele%value(voltage$) /= 0) then
    phase = twopi * ele%value(phi0$) + end%vec(5)  &
                      / (ele%value(rf_wavelength$) * beta_s)
    en_gain = ele%value(voltage$) * sin(twopi * phase)
    if ((en_gain + gamma_b * m_electron) <= m_electron) then
      param%lost = .true.
      return
    else
      end%vec(6) = end%vec(6) + en_gain / ele%value(beam_energy$)
      c_e = en_gain / (m_electron * length)
    endif
  else
    c_e = 0.0
  endif

! Beginning fringe

  c_m = param%particle * c_light * ele%value(b_z$) / m_electron
  call mat_make_unit(mat4)
  mat4(2,3) = c_m / 2 *  &
                  sqrt((1 + end%vec(2)**2 + end%vec(4)**2) / (gamma_b**2 - 1))
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
                                    b_x(j), b_y(j), end%vec, mat4, vec_st)
        end%vec(1:4) = matmul(mat4, end%vec(1:4))
        do i = 1, 4
          end%vec(i) = end%vec(i) + param%particle * vec_st(i)
        enddo

! Calculation of end.z.pos

        if (abs(c_e) > 0.001) then
          end%vec(5) = end%vec(5) - ll + (sqrt(1 + (c_e * ll  &
                          + sqrt(gamma_old**2 - 1))**2) - gamma_old) / c_e
        else
          end%vec(5) = end%vec(5) + ll * (sqrt(1 - 1/gamma_old**2) - 1  &
              - c_e * ll * (gamma_old**2 - 0.5)/gamma_old**3)
        endif

        if (param%aperture_limit_on) then
          if (check_x_ap) then
            x_lim = x_lim_chng_rate * s_grand_cum + x_beg_lim
            if (abs(end%vec(1)) > x_lim) param%lost = .true.
          endif
          if (check_y_ap) then
            y_lim = y_lim_chng_rate * s_grand_cum + y_beg_lim
            if (abs(end%vec(3)) > y_lim) param%lost = .true.
          endif
          if (param%lost) return
        endif

      endif
    enddo

  enddo

! Ending fringe

  call mat_make_unit(mat4)
  mat4(4,1) = c_m / 2 *  &
                sqrt((1 + end%vec(2)**2 + end%vec(4)**2) / (gamma_new**2 - 1))
  mat4(2,3) = -mat4(4,1)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

  call offset_particle (ele, param, end, unset$)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine change_basis (coord, ref_energy, ref_z, to_cart, time_disp)
!
!   Subroutine to convert accelerator coordinates (x, x', y, y', z, z') to
! cartesian coordinates and time derivatives (x, x_dot, y, y_dot, z, z_dot) or
! vice-versa.  The conversions assume no small-angle nor other approximations
! whatsoever EXCEPT that there is no (local) bend--the design orbit is locally
! straight.
!   Note that if the Lorentz factor (gamma) is around 100 or higher, converting
! _from_ accelerator coordinates and then back _to_ accelerator coordinates will
! produce a slight change in z'.
!   Also, due to rounding errors in iterations, particles may acquire a speed
! higher than c.  Such speeds will be reset (if TO_CART is false) so that the
! Lorentz factor is about 1000.
!  These two problems are corrected if coil_track is the calling routine.
! -- Created by Daniel Fromowitz, November 1998.
!
! Modules Needed:
!   use bmad
!
! Input:
!   coord      -- Coord_struct: Coordinates of particle
!   ref_energy -- Real(rp): Reference energy of beam
!   ref_z      -- Real(rp): Reference longitudinal position of beam
!   to_cart    -- Logical: True if converting to cartesian coordinates
!                          False if converting to accelerator coordinates
!   If to_cart == .false.:
!     time_disp -- Real(rp): Time displacement of particle
!
! Output:
!   coord  -- Coord_struct: Converted coordinates
!   If to_cart == .true.:
!     time_disp -- Real(rp): Time displacement of particle
!-

subroutine change_basis (coord, ref_energy, ref_z, to_cart, time_disp)

  implicit none

  type (coord_struct)  coord, saved
  real(rp) beta2, ref_energy, ref_z, time_disp, rvar
  logical to_cart

!

  saved%vec(1) = coord%vec(1)
  saved%vec(3) = coord%vec(3)

  if (to_cart) then
    rvar = (m_electron / (ref_energy * (coord%vec(6)+1)))**2
    if (rvar <= 0.001) then
      saved%vec(6) = c_light / sqrt(1 + coord%vec(2)**2 + coord%vec(4)**2)  &
        * (1 - rvar/2)
    else
      saved%vec(6) = c_light * sqrt((1 - rvar)  &
        / (1 + coord%vec(2)**2 + coord%vec(4)**2))
    endif
    saved%vec(2) = coord%vec(2) * saved%vec(6)
    saved%vec(4) = coord%vec(4) * saved%vec(6)
!         Convert z from a relative position to an absolute position:
    saved%vec(5) = ref_z
    time_disp =  - coord%vec(5) / saved%vec(6)
  else
    saved%vec(2) = coord%vec(2) / coord%vec(6)
    saved%vec(4) = coord%vec(4) / coord%vec(6)
    beta2 = (coord%vec(2)**2 + coord%vec(4)**2 + coord%vec(6)**2) / c_light**2
!         Prevent sqrt(negative number) due to rounding errors.  Set such
!         particles to have gamma = 1000 (approximately).
    if (beta2 >= 0.999999) beta2 = 0.999999
    saved%vec(6) = m_electron / (ref_energy * sqrt(1 - beta2)) - 1
!         Convert z from an absolute position to a relative position:
    saved%vec(5) = - coord%vec(6) * time_disp
  endif
  coord = saved

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine ACCEL_SOL_MAT_CALC (LS, C_M, C_E, GAMMA_OLD, GAMMA_NEW, B_X, B_Y,
!   COORD, MAT4, VEC_ST)
!
!   Subroutine to calculate the 4x4 transfer matrix (excluding steerings) for a
! segment of an accelerating solenoid.  A vector is also calculated for the
! steerings.
! -- Created by Daniel Fromowitz, September 1999.
!
! Input:
!     LS        -- Real(rp): length of the segment
!     C_M       -- Real(rp): constant proportional to the longitudinal magnetic
!                         field
!     C_E       -- Real(rp): constant proportional to the electric field
!     GAMMA_OLD -- Real(rp): Lorentz factor at beginning of segment
!     GAMMA_NEW -- Real(rp): Lorentz factor at end of segment
!     B_X       -- Real(rp): Horizontal field of transverse steering
!     B_Y       -- Real(rp): Vertical field of transverse steering
!     COORD(6)  -- Real(rp): Starting position
!
! Output:
!     MAT4(4,4) -- Real(rp): 4x4 transfer matrix excluding steerings
!     VEC_ST(4) -- Real(rp): Vector due to steerings (assuming positrons)
!-

subroutine accel_sol_mat_calc (ls, c_m, c_e, gamma_old, gamma_new, b_x,  &
                                                     b_y, coord, mat4, vec_st)

  implicit none

  
  real(rp) ls, c_m, c_e, gamma_old, gamma_new, b_x, b_y, mat4(4,4), vec_st(4)
  real(rp) coef, cosr, sinr, denom, ratio, ratio_c_m, sinr_c_m, onecosr_c_m
  real(rp) mat_st(4,2), coord(6)
  integer i

  if (abs(c_e) > 0.001) then
    ratio_c_m = log(gamma_new / gamma_old) / c_e
    ratio = c_m * ratio_c_m
  else
    ratio_c_m = ls / gamma_old * (1 - c_e * ls / (2 * gamma_old))
    ratio = c_m * ratio_c_m
  endif
  if (abs(c_m) > 0.001) then
    sinr_c_m = sin(ratio) / c_m
    onecosr_c_m = (1 - cos(ratio)) / c_m
  else
    sinr_c_m = ratio_c_m
    onecosr_c_m = c_m * ratio_c_m**2 / 2
  endif
  sinr = sin(ratio)
  cosr = cos(ratio)

  mat4(1,1) = 1
  mat4(1,2) = gamma_old * sinr_c_m
  mat4(1,3) = 0
  mat4(1,4) = gamma_old * onecosr_c_m
  mat4(2,1) = 0
  mat4(2,2) = cos(ratio) * gamma_old / gamma_new
  mat4(2,3) = 0
  mat4(2,4) = sin(ratio) * gamma_old / gamma_new
  mat4(3,1) = 0
  mat4(3,2) = -mat4(1,4)
  mat4(3,3) = 1
  mat4(3,4) = mat4(1,2)
  mat4(4,1) = 0
  mat4(4,2) = -mat4(2,4)
  mat4(4,3) = 0
  mat4(4,4) = mat4(2,2)

!     Steerings:

  if ((b_x /= 0.0) .or. (b_y /= 0.0)) then
    denom = c_e**2 + c_m**2
    if (denom > 2.e-6) then
      coef = c_light / m_electron / denom
      mat_st(1,1) = coef *  &
                    (c_m * ls - gamma_old * (c_e * onecosr_c_m + sinr))
      mat_st(1,2) = coef * (gamma_old * (cosr + c_e * sinr_c_m) - gamma_new)
      mat_st(2,1) = coef *  &
                   (c_m - gamma_old / gamma_new * (c_e * sinr + c_m * cosr))
      mat_st(2,2) = coef *  &
                   (gamma_old / gamma_new * (c_e * cosr - c_m * sinr) - c_e)
    else
      coef = c_light / m_electron  &
        * sqrt((1 + coord(2)**2 + coord(4)**2) / (gamma_old**2 - 1))
      mat_st(1,1) = 0
      mat_st(1,2) = -coef * ls**2 / 2
      mat_st(2,1) = 0
      mat_st(2,2) = -coef * ls
    endif
    mat_st(3,1) = -mat_st(1,2)
    mat_st(3,2) =  mat_st(1,1)
    mat_st(4,1) = -mat_st(2,2)
    mat_st(4,2) =  mat_st(2,1)
    do i = 1, 4
      vec_st(i) = mat_st(i,1) * b_x + mat_st(i,2) * b_y
    enddo
  else
    do i = 1, 4
      vec_st(i) = 0
    enddo
  endif

end subroutine

end module
