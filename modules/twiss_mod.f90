#include "CESR_platform.inc" 

module twiss_mod

  use matrix_mod

  type twiss_struct
    real(rp) beta, alpha, gamma, phi, eta, etap
    real(rp) sigma, sigma_p, emit, norm_emit
  end type

contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_symp_decouple (t0, tol, stat, u, v, ubar, vbar, g,
!                                                     twiss1, twiss2, type_out)
!
! Subroutine to find the symplectic eigen modes of the
! one turn 4x4 coupled transfer matrix T0.
!
! Modules needed:
!   use bmad
!
! Input:
!   t0(4,4)  -- Real(rp): Input matrix
!   type_out -- Logical: If .true. then an error message is typed out
!               for a non ok$ STAT
!   tol      -- Real(rp): tollerence for nonsymplectiy
!
! Output:
!   stat    -- Integer: status of results:
!                       OK$, IN_STOP_BAND$, UNSTABLE$, NON_SYMPLECTIC$
!   twiss1  -- Twiss_struct: Twiss params for the "upper left" mode.
!      %phi   -- Rotation angle in radians, 0 < %PHI < twopi
!   twiss2  -- Twiss_struct: Twiss params for the "lower right" mode.
!   u(4,4), v(4,4), ubar(4,4), vbar(4,4), g(4,4) -- Real(rp): See
!                 MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.
!
!-

subroutine mat_symp_decouple(t0, tol, stat, U, V, Ubar, Vbar, G,  &
                                             twiss1, twiss2, type_out)

  implicit none

  type (twiss_struct)  twiss1, twiss2

  integer i, j, stat

  real(rp) t0(4,4), unit4(4,4), U(4,4), V(4,4), V_inv(4,4)
  real(rp) Ubar(4,4), Vbar(4,4), G(4,4), G_inv(4,4)
  real(rp) t0_11(2,2), t0_12(2,2), t0_21(2,2), t0_22(2,2)
  real(rp) c(2,2), c_conj(2,2), H(2,2), temp2(2,2)
  real(rp) g1(2,2), g2(2,2), g1_inv(2,2), g2_inv(2,2)
  real(rp) gamma, det_H, det,  trace_t0_diff, denom
  real(rp) scaler, tol

  logical type_out

! define some matraces
! remember that array storage is column first!

  data unit4 /  1,  0,  0,  0,  &
             0,  1,  0,  0,  &
             0,  0,  1,  0,  &
             0,  0,  0,  1 /

! check input matrix

  if (mat_symp_error (t0) > tol) then
    stat = non_symplectic$
    if (type_out) then
      print *, 'ERROR IN MAT_SYMP_DECOUPLE: NON-SYMPLECTIC INPUT MATRIX'
    endif
    return
  endif

! load submatrices

  do i = 1, 2
    do j = 1, 2
      t0_11(i, j) = t0(i, j)       ! = M matrix (MGB eq 6)
      t0_12(i, j) = t0(i, j+2)     ! = m matrix
      t0_21(i, j) = t0(i+2, j)     ! = n matrix
      t0_22(i, j) = t0(i+2, j+2)   ! = N matrix
    enddo
  enddo

! Construct H matrix (MGB eq 12)

  call mat_symp_conj (t0_21, temp2)
  H = t0_12 + temp2


! Compute traces and
! compute DET_H and determine if we are in a stop band (MGB Eq. 14)

  trace_t0_diff = (t0_11(1,1) + t0_11(2,2)) - (t0_22(1,1) + t0_22(2,2))
  call mat_det (H, det_H)
  denom = trace_t0_diff**2 + 4.0 * det_H

  if (denom <= 0) then
    stat = in_stop_band$
    u(1,1) = 1 - denom   ! fake so matrix looks unstable
    u(2,2) = 1 - denom
    return
  endif

! Compute GAMMA (MGB Eq. 14)

  gamma = sqrt(0.5 + 0.5 * sqrt(trace_t0_diff**2 / denom))

! Construct C matrix (MGB Eq. 13 with Eq. 14) and symplectic conjugate.

  scaler = -sign(1.0_rp, trace_t0_diff) / (gamma * sqrt(denom))
  c = scaler * H
  call mat_symp_conj (c, c_conj)

! Compute matrix V and inverse V_INV (MGB Eq. 10)

  V = gamma * unit4
  V(1:2,3:4) = c
  V(3:4,1:2) = -c_conj

  call mat_symp_conj (V, V_inv)

! Compute uncoupled matrix U (MGB Eq. 10)

  U = matmul (matmul (V_inv, t0), V)

  if (mat_symp_error(U) > tol) then
    stat = non_symplectic$
    if (type_out) then
      print *, 'ERROR IN MAT_SYMP_DECOUPLE: NON-SYMPLECTIC U MATRIX'
    endif
    return
  endif

! check that the eigen modes are stable

  if (abs(U(1,1) + U(2,2)) > 2 .or.  &
                               abs(U(3,3) + U(4,4)) > 2) then
    stat = unstable$
    return
  endif

! calculate twiss parameters for U sub matrices

  call twiss_from_mat2 (U(1:2,1:2), det, twiss1, stat, tol, .false.)
  if (stat /= ok$) then
    if (type_out) print *,  &
      'ERROR IN MAT_SYMP_DECOUPLE: UNABLE TO COMPUTE "A" mode TWISS'
    return
  endif

  call twiss_from_mat2 (U(3:4,3:4), det, twiss2, stat, tol, .false.)
  if (stat /= ok$) then
    if (type_out) print *,  &
      'ERROR IN MAT_SYMP_DECOUPLE: UNABLE TO COMPUTE "B" mode TWISS'
    return
  endif

! Compute normalized uncoupled matrix Ubar

! First compute G matrix and inverse G_INV (PPB/DLR. Eq. 3+)

  g1(1, 1) = 1 / sqrt(twiss1%beta)
  g1(1, 2) = 0.0
  g1(2, 1) = twiss1%alpha / sqrt(twiss1%beta)
  g1(2, 2) = sqrt(twiss1%beta)

  g2(1, 1) = 1 / sqrt(twiss2%beta)
  g2(1, 2) = 0.0
  g2(2, 1) = twiss2%alpha / sqrt(twiss2%beta)
  g2(2, 2) = sqrt(twiss2%beta)

  G = 0
  G(1:2,1:2) = g1
  G(3:4,3:4) = g2

  call mat_symp_conj (g1, g1_inv)
  call mat_symp_conj (g2, g2_inv)

  G_inv = 0
  G_inv(1:2,1:2) = g1_inv
  G_inv(3:4,3:4) = g2_inv

! Compute Vbar (PPB/DLR Eq. 3++)

  Vbar = matmul (matmul (G, V), G_inv)

! compute Ubar

  Ubar = matmul (matmul (G, U), G_inv)

  if (mat_symp_error (Ubar) > tol) then
    stat = non_symplectic$
    if (type_out) then
      print *, 'ERROR IN MAT_SYMP_DECOUPLE: NON-SYMPLECTIC UBAR MATRIX'
    endif
    return
  endif

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine twiss_from_mat2 (mat, det, twiss, stat, tol, type_out)
!
! Subroutine to extract the twiss parameters from one-turn
! 2x2 matrices
!
! Modules needed:
!   use bmad
!
! Input:
!   mat(:, :)   -- Real(rp): Input matrix
!   type_out    -- Logical: If .true. then an error message is typed out
!                            for a non ok$ STAT
!   tol         -- Real(rp): tollerence for nonsymplectiy
!
! Output:
!   stat -- Integer: status of results:
!                         OK$, UNSTABLE$, NON_SYMPLECTIC$
!   det    -- Real(rp): Determinant of matrix. Should be = 1
!   twiss  -- Twiss_struct: Twiss parameters
!                            TWISS.PHI is in radians, 0 < TWISS.PHI < twopi
!-

subroutine twiss_from_mat2 (mat, det, twiss, stat, tol, type_out)

  implicit none

  type (twiss_struct)  twiss

  integer stat
  real(rp) mat(:, :), t_cos, t_sin, det, tol, radical
  logical type_out

!

  stat = ok$

  t_cos = (mat(1,1) + mat(2,2)) / 2.0
  det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)

  if (abs(t_cos) >= 1.0) then
    if (type_out) print *, 'ERROR IN TWISS_FROM_MAT: UNSTABLE MATRIX'
    stat = unstable$
    return
  endif

  if (abs(det - 1) > tol) then
    if (type_out) print *,  &
                    'WARNING IN TWISS_FROM_MAT: MATRIX DETERMINANT >< 1:', det
    stat = non_symplectic$
  endif

  t_sin = sign(1.0_rp, mat(1,2)) * sqrt(1.0 - t_cos**2)
  twiss%phi = atan2(t_sin, t_cos)
  if (twiss%phi < 0) twiss%phi = twiss%phi + twopi

  if (abs(t_sin) < 1.0e-7) then
    twiss%alpha = 0.0
    twiss%beta = 1.0
    twiss%gamma = 1.0
  else
    twiss%alpha = (mat(1,1) - mat(2,2)) / (2.0 * t_sin)
    radical = -(1 + twiss%alpha**2) * mat(1,2) / mat(2,1)
    if (radical <= 0) then
      if (type_out) print *, 'ERROR IN TWISS_FROM_MAT: NON-SYMPLECTIC MATRIX'
      stat = non_symplectic$
      return
    endif
    twiss%beta = sqrt(radical)
    twiss%gamma = (1 + twiss%alpha**2) / twiss%beta
  endif

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine twiss_to_1_turn_mat (twiss, phi, mat2)
!
! Subroutine to form the 2x2 1-turn transfer matrix from the twiss parameters.
!
! Modules needed:
!   use bmad
!
! Input:
!   twiss -- Twiss_struct: Structure holding the Twiss parameters.
!     %beta
!     %alpha
!   phi   -- Real(rp): Tune in radians.
!
! Output:
!   mat2(2,2) -- Real(rp): 1-turn matrix.
!-

subroutine twiss_to_1_turn_mat (twiss, phi, mat2)

  implicit none

  type (twiss_struct) twiss

  real(rp) phi, mat2(2,2), c, s

!

  c = cos(phi)
  s = sin(phi)

  mat2(1,1) =  c + s * twiss%alpha
  mat2(1,2) =  s * twiss%beta
  mat2(2,1) = -s * (1 + twiss%alpha**2) / twiss%beta
  mat2(2,2) =  c - s * twiss%alpha

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine make_g2_mats (twiss, g_mat, g_inv_mat)
!
! Subroutine make the matrices needed to go from normal mode coords
! to coordinates with the beta function removed.
!
! Modules Needed:
!   use bmad
!
! Input:
!   twiss        -- Twiss_struct: Twiss parameters.
!
! Output:
!   g_mat(2,2)     -- Real(rp): Normal mode to betaless coords.
!   g_inv_mat(2,2) -- Real(rp): The inverse of g_mat.
!-

subroutine make_g2_mats (twiss, g2_mat, g2_inv_mat)

  implicit none

  type (twiss_struct) twiss

  real(rp) g2_mat(2,2), g2_inv_mat(2,2)
  real(rp) sqrt_beta, alpha
!

  sqrt_beta = sqrt(twiss%beta)
  alpha     = twiss%alpha

  g2_mat(1,1) = 1 / sqrt_beta
  g2_mat(1,2) = 0
  g2_mat(2,1) = alpha / sqrt_beta
  g2_mat(2,2) = sqrt_beta

  g2_inv_mat = g2_mat
  g2_inv_mat(2,1) = -g2_mat(2,1)

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine transfer_mat2_from_twiss (twiss1, twiss2, mat)
!
! Subroutine to make a 2 x 2 transfer matrix from the twiss parameters
! at the end points.
!
! Modules Needed:
!   use bmad
!
! Input:
!   twiss1  -- Twiss_struct: Twiss parameters at the initial point.
!     %beta   -- Beta parameter.
!     %alpha  -- Alpha parameter.
!     %phi    -- Phase at initial point.
!   twiss2  -- Twiss_struct: Twiss parameters at the end point.
!     %beta   -- Beta parameter.
!     %alpha  -- Alpha parameter.
!     %phi    -- Phase at final point.
!
! Output:
!   mat(2,2) -- Real(rp): Transfer matrix between the two points.
!-

subroutine transfer_mat2_from_twiss (twiss1, twiss2, mat)

  implicit none

  type (twiss_struct) twiss1, twiss2

  real(rp) mat(2,2), a1, a2, b1, b2, sin21, cos21

!

  sin21 = sin(twiss2%phi - twiss1%phi)
  cos21 = cos(twiss2%phi - twiss1%phi)
  b1 = twiss1%beta;  b2 = twiss2%beta
  a1 = twiss1%alpha; a2 = twiss2%alpha

  mat(1,1) = sqrt(b2/b1) * (cos21 + a1 * sin21)
  mat(1,2) = sqrt(b1 * b2) * sin21
  mat(2,1) = -((a2-a1) * cos21 + (1 + a1 * a2) * sin21) / (sqrt(b1 * b2))
  mat(2,2) = sqrt(b1/b2) * (cos21 - a2 * sin21)

end subroutine

end module
