!+
! Subroutine mat_symp_decouple (t0, stat, u, v, ubar, vbar, g, twiss1, twiss2, gamma, type_out)
!
! Subroutine to find the symplectic eigen modes of the
! one turn 4x4 coupled transfer matrix T0.
!
! Input:
!   t0(4,4)  -- Real(rp): Input matrix
!   type_out -- Logical: If .true. then an error message is typed out
!               for a non ok$ STAT
!
! Output:
!   stat      -- Integer: status of results: ok$, in_stop_band$, or unstable$
!   twiss1    -- Twiss_struct: Twiss params for the "upper left" mode.
!      %phi     -- Rotation angle in radians, 0 < %PHI < twopi
!   twiss2    -- Twiss_struct: Twiss params for the "lower right" mode.
!   u(4,4), v(4,4), ubar(4,4), vbar(4,4), g(4,4) 
!             -- Real(rp): See MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.
!   gamma     -- Real(rp): gamma_c factor.
!-

subroutine mat_symp_decouple(t0, stat, U, V, Ubar, Vbar, G,  twiss1, twiss2, gamma, type_out)

use equal_mod, dummy => mat_symp_decouple

implicit none

type (twiss_struct)  twiss1, twiss2

integer i, j, stat

real(rp) t0(4,4), U(4,4), V(4,4)
real(rp) Ubar(4,4), Vbar(4,4), G(4,4), G_inv(4,4)
real(rp) t0_11(2,2), t0_12(2,2), t0_21(2,2), t0_22(2,2)
real(rp) c(2,2), H(2,2)
real(rp) g1(2,2), g2(2,2)
real(rp) gamma, det_H, det,  trace_t0_diff, denom, scalar
logical type_out
character(20), parameter :: r_name = 'mat_symp_decouple'

! Load submatrices

do i = 1, 2
  do j = 1, 2
    t0_11(i, j) = t0(i, j)       ! = M matrix (MGB eq 6)
    t0_12(i, j) = t0(i, j+2)     ! = m matrix
    t0_21(i, j) = t0(i+2, j)     ! = n matrix
    t0_22(i, j) = t0(i+2, j+2)   ! = N matrix
  enddo
enddo

!--------------

if (all(t0(1:2,3:4) == 0)) then  ! Uncoupled
  gamma = 1
  c = 0

else
  ! Construct H matrix (MGB eq 12)

  H = t0_12 + mat_symp_conj (t0_21)

  ! Compute traces and
  ! Compute DET_H and determine if we are in a stop band (MGB Eq. 14)

  trace_t0_diff = (t0_11(1,1) + t0_11(2,2)) - (t0_22(1,1) + t0_22(2,2))
  det_h = determinant (H)
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

  scalar = -sign(1.0_rp, trace_t0_diff) / (gamma * sqrt(denom))
  c = scalar * H
endif

! Compute matrix V (MGB Eq. 10)

V = 0
forall (i = 1:4) v(i,i) = gamma 
V(1:2,3:4) = c
V(3:4,1:2) = -mat_symp_conj(c)

! Compute uncoupled matrix U (MGB Eq. 10)

U = matmul (matmul (mat_symp_conj(V), t0), V)

! check that the eigen modes are stable

if (abs(U(1,1) + U(2,2)) > 2 .or. abs(U(3,3) + U(4,4)) > 2) then
  if (abs(U(1,1) + U(2,2)) > 2 .and. abs(U(3,3) + U(4,4)) > 2) then
    stat = unstable$
  elseif (abs(U(1,1) + U(2,2)) > 2) then
    stat = unstable_a$
  else
    stat = unstable_b$
  endif
  return
endif

! calculate twiss parameters for U sub matrices

call twiss_from_mat2 (U(1:2,1:2), twiss1, stat, .false.)
if (stat /= ok$) then
  stat = unstable_a$
  if (type_out) call out_io (s_warn$, r_name, 'UNABLE TO COMPUTE "A" mode TWISS')
  return
endif

call twiss_from_mat2 (U(3:4,3:4), twiss2, stat, .false.)
if (stat /= ok$) then
  stat = unstable_b$
  if (type_out) call out_io (s_warn$, r_name, 'UNABLE TO COMPUTE "B" mode TWISS')
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

G_inv = 0
G_inv(1:2,1:2) = mat_symp_conj(g1)
G_inv(3:4,3:4) = mat_symp_conj(g2)

! Compute Vbar (PPB/DLR Eq. 3++)

Vbar = matmul (matmul (G, V), G_inv)

! compute Ubar

Ubar = matmul (matmul (G, U), G_inv)

end subroutine mat_symp_decouple

