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
!     LS        -- Real: length of the segment
!     C_M       -- Real: constant proportional to the longitudinal magnetic
!                         field
!     C_E       -- Real: constant proportional to the electric field
!     GAMMA_OLD -- Real: Lorentz factor at beginning of segment
!     GAMMA_NEW -- Real: Lorentz factor at end of segment
!     B_X       -- Real: Horizontal field of transverse steering
!     B_Y       -- Real: Vertical field of transverse steering
!     COORD     -- Coord_struct: Starting position
!
! Output:
!     MAT4(4,4) -- Real: 4x4 transfer matrix excluding steerings
!     VEC_ST(4) -- Real: Vector due to steerings (assuming positrons)
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:47  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine accel_sol_mat_calc (ls, c_m, c_e, gamma_old, gamma_new, b_x,  &
    b_y, coord, mat4, vec_st)
  use bmad_struct
  implicit none

  type (coord_struct)  coord
  real ls, c_m, c_e, gamma_old, gamma_new, b_x, b_y, mat4(4,4), vec_st(4)
  real coef, cosr, sinr, denom, ratio, ratio_c_m, sinr_c_m, onecosr_c_m
  real mat_st(4,2)
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
      coef = c_light / (e_mass * 1.e9) / denom
      mat_st(1,1) = coef *  &
                    (c_m * ls - gamma_old * (c_e * onecosr_c_m + sinr))
      mat_st(1,2) = coef * (gamma_old * (cosr + c_e * sinr_c_m) - gamma_new)
      mat_st(2,1) = coef *  &
                   (c_m - gamma_old / gamma_new * (c_e * sinr + c_m * cosr))
      mat_st(2,2) = coef *  &
                   (gamma_old / gamma_new * (c_e * cosr - c_m * sinr) - c_e)
    else
      coef = c_light / (e_mass * 1.e9)  &
        * sqrt((1 + coord%x%vel**2 + coord%y%vel**2) / (gamma_old**2 - 1))
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

  return
  end
