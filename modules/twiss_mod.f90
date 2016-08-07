module twiss_mod

use bmad_struct

contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_symp_decouple (t0, stat, u, v, ubar, vbar, g, twiss1, twiss2, gamma, type_out)
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

! check input matrix

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

H = t0_12 + mat_symp_conj (t0_21)

! Compute traces and
! compute DET_H and determine if we are in a stop band (MGB Eq. 14)

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

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine twiss_from_mat2 (mat, twiss, stat, type_out)
!
! Subroutine to extract the twiss parameters from one-turn
! 2x2 matrices
!
! Modules needed:
!   use bmad
!
! Input:
!   mat(:, :)   -- Real(rp): Input matrix.
!   type_out    -- Logical: If .true. then an error message is typed out
!                            for a non ok$ STAT
!
! Output:
!   stat -- Integer: status of results:
!                         OK$, UNSTABLE$, NON_SYMPLECTIC$
!   twiss  -- Twiss_struct: Twiss parameters
!                            TWISS.PHI is in radians, 0 < TWISS.PHI < twopi
!-

subroutine twiss_from_mat2 (mat_in, twiss, stat, type_out)

implicit none

type (twiss_struct)  twiss

integer stat
real(rp) mat_in(:, :), t_cos, t_sin, det, radical, mat(2,2)
logical type_out
character(16), parameter :: r_name = 'twiss_from_mat2'

!

stat = ok$

det = mat_in(1,1) * mat_in(2,2) - mat_in(1,2) * mat_in(2,1)
mat = mat_in / det

t_cos = (mat(1,1) + mat(2,2)) / 2.0
if (abs(t_cos) >= 1.0) then
  if (type_out) call out_io (s_warn$, r_name, 'UNSTABLE MATRIX')
  stat = unstable$
  return
endif

t_sin = sign(1.0_rp, mat(1,2)) * sqrt(1.0 - t_cos**2)
twiss%phi = atan2(t_sin, t_cos)
if (twiss%phi < 0) twiss%phi = twiss%phi + twopi

twiss%alpha = (mat(1,1) - mat(2,2)) / (2.0 * t_sin)
radical = -(1 + twiss%alpha**2) * mat(1,2) / mat(2,1)
if (radical <= 0) then
  if (type_out) call out_io (s_warn$, r_name, 'NON-SYMPLECTIC MATRIX')
  stat = non_symplectic$
  return
endif
twiss%beta = sqrt(radical)
twiss%gamma = (1 + twiss%alpha**2) / twiss%beta

end subroutine twiss_from_mat2

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

end subroutine twiss_to_1_turn_mat

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

end subroutine make_g2_mats

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

end subroutine transfer_mat2_from_twiss

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function average_twiss (frac1, twiss1, twiss2) result (ave_twiss)
!
! Routine to average twiss parameters.
!
! Input:
!   frac1  -- real(rp): Fraction of twiss1 to use in the average.
!   twiss1 -- twiss_struct: Twiss parameters to average.
!   twiss1 -- twiss_struct: Twiss parameters to average.
!
! Output:
!   ave_twiss -- twiss_struct: Average twiss.
!-

function average_twiss (frac1, twiss1, twiss2) result (ave_twiss)

implicit none

type (twiss_struct) twiss1, twiss2, ave_twiss
real(rp) frac1

!

ave_twiss%beta      = frac1 * twiss1%beta      + (1-frac1) * twiss2%beta
ave_twiss%alpha     = frac1 * twiss1%alpha     + (1-frac1) * twiss2%alpha
ave_twiss%gamma     = frac1 * twiss1%gamma     + (1-frac1) * twiss2%gamma
ave_twiss%phi       = frac1 * twiss1%phi       + (1-frac1) * twiss2%phi
ave_twiss%eta       = frac1 * twiss1%eta       + (1-frac1) * twiss2%eta
ave_twiss%etap      = frac1 * twiss1%etap      + (1-frac1) * twiss2%etap
ave_twiss%sigma     = frac1 * twiss1%sigma     + (1-frac1) * twiss2%sigma
ave_twiss%sigma_p   = frac1 * twiss1%sigma_p   + (1-frac1) * twiss2%sigma_p
ave_twiss%emit      = frac1 * twiss1%emit      + (1-frac1) * twiss2%emit
ave_twiss%norm_emit = frac1 * twiss1%norm_emit + (1-frac1) * twiss2%norm_emit

end function average_twiss

end module
