!+
! Subroutine MAT_SYMP_DECOUPLE (T0, TOL, STAT, U, V, UBAR, VBAR, G,
!                                                     TWISS1, TWISS2, TYPE_OUT)
!
! Subroutine to find the symplectic eigen modes of the
! one turn 4x4 coupled transfer matrix T0.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   T0(4,4)  -- Real: Input matrix
!   TYPE_OUT -- Logical: If .true. then an error message is typed out
!               for a non ok$ STAT
!   TOL      -- Real: tollerence for nonsymplectiy
!
! Output:
!   STAT    -- Integer: status of results:
!                       OK$, IN_STOP_BAND$, UNSTABLE$, NON_SYMPLECTIC$
!   TWISS1  -- Twiss_struct: Twiss params for the "upper left" mode.
!      %PHI   -- Rotation angle in radians, 0 < %PHI < twopi
!   TWISS2  -- Twiss_struct: Twiss params for the "lower right" mode.
!   U(4,4), V(4,4), UBAR(4,4), VBAR(4,4), G(4,4) -- Real: See
!                 MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.
!
!-


subroutine mat_symp_decouple(t0, tol, stat, U, V, Ubar, Vbar, G,  &
                                             twiss1, twiss2, type_out)

  use bmad_struct
  use bmad_interface

  implicit none

  type (twiss_struct)  twiss1, twiss2

  integer i, j, stat

  real t0(4,4), unit4(4,4), U(4,4), V(4,4), V_inv(4,4)
  real temp4(4,4), tol, error
  real Ubar(4,4), Vbar(4,4), G(4,4), G_inv(4,4)
  real t0_11(2,2), t0_12(2,2), t0_21(2,2), t0_22(2,2)
  real c(2,2), c_conj(2,2), H(2,2), temp2(2,2)
  real g1(2,2), g2(2,2), g1_inv(2,2), g2_inv(2,2)
  real gamma, det_H, det,  trace_t0_diff, denom
  real scaler

  logical type_out

! define some matraces
! remember that array storage is column first!

  data unit4 /  1,  0,  0,  0,  &
             0,  1,  0,  0,  &
             0,  0,  1,  0,  &
             0,  0,  0,  1 /

! check input matrix

  call mat_symp_check (t0, error)
  if (error > tol) then
    stat = non_symplectic$
    if (type_out) then
      type *, 'ERROR IN MAT_SYMP_DECOUPLE: NON-SYMPLECTIC INPUT MATRIX'
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

  call mat_symp_conj (t0_21, temp2, 2, 2)
  H = t0_12 + temp2


! Compute traces and
! compute DET_H and determine if we are in a stop band (MGB Eq. 14)

  trace_t0_diff = (t0_11(1,1) + t0_11(2,2)) - (t0_22(1,1) + t0_22(2,2))
  call mat_det (H, det_H, 2, 2)
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

  scaler = -sign(1.0, trace_t0_diff) / (gamma * sqrt(denom))
  c = scaler * H
  call mat_symp_conj (c, c_conj, 2, 2)

! Compute matrix V and inverse V_INV (MGB Eq. 10)

  V = gamma * unit4
  V(1:2,3:4) = c
  V(3:4,1:2) = -c_conj

  call mat_symp_conj (V, V_inv, 4, 4)

! Compute uncoupled matrix U (MGB Eq. 10)

  U = matmul (matmul (V_inv, t0), V)

  call mat_symp_check (U, error)
  if (error > tol) then
    stat = non_symplectic$
    if (type_out) then
      type *, 'ERROR IN MAT_SYMP_DECOUPLE: NON-SYMPLECTIC U MATRIX'
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
    if (type_out) type *,  &
      'ERROR IN MAT_SYMP_DECOUPLE: UNABLE TO COMPUTE "A" mode TWISS'
    return
  endif

  call twiss_from_mat2 (U(3:4,3:4), det, twiss2, stat, tol, .false.)
  if (stat /= ok$) then
    if (type_out) type *,  &
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

  call mat_symp_conj (g1, g1_inv, 2, 2)
  call mat_symp_conj (g2, g2_inv, 2, 2)

  G_inv = 0
  G_inv(1:2,1:2) = g1_inv
  G_inv(3:4,3:4) = g2_inv

! Compute Vbar (PPB/DLR Eq. 3++)

  Vbar = matmul (matmul (G, V), G_inv)

! compute Ubar

  Ubar = matmul (matmul (G, U), G_inv)

  call mat_symp_check (Ubar, error)
  if (error > tol) then
    stat = non_symplectic$
    if (type_out) then
      type *, 'ERROR IN MAT_SYMP_DECOUPLE: NON-SYMPLECTIC UBAR MATRIX'
    endif
    return
  endif

!

  return
  end


