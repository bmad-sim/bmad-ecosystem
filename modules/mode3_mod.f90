module mode3_mod

use twiss_and_track_mod

implicit none

real(rp), parameter :: m = 0.707106781d0  ! 1/sqrt(2)
real(rp), parameter :: o = 0.0d0  ! for compact code
real(rp), parameter :: l = 1.0d0  ! for compact code

real(rp), parameter :: Qr(6,6) = reshape( [m, m, o, o, o, o, o, o, o, o, o, o, o, o, m, m, o, o, &
                                           o, o, o, o, o, o, o, o, o, o, m, m, o, o, o, o, o, o], [6,6] )
real(rp), parameter :: Qi(6,6) = reshape( [o, o, o, o, o, o, m, -m, o, o, o, o, o, o, o, o, o, o, &
                                           o, o, m, -m, o, o, o, o, o, o, o, o, o, o, o, o, m, -m], [6,6] )
real(rp), parameter :: S(6,6) = reshape( [o, -l, o, o, o, o, l, o, o, o, o, o, &
                                          o, o, o, -l, o, o, o, o, l, o, o, o, &
                                          o, o, o, o, o, -l, o, o, o, o, l, o], [6,6] )
real(rp), parameter :: I2(2, 2) = reshape( [1, 0, 0, 1], [2, 2] )

private m, o, l
private Qr, Qi
private S

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine t6_to_B123(N, abz_tunes, B1, B2, B3)
!
! This decomposes the one-turn matrix according to Equation 56 from 
! "Alternative approach to general coupled linear optics" by A. Wolski. PRSTAB.
!
! Note that a sigma matrix can be assembeled from:  sigma = B1*emit_a + B2*emit_b + B3*emit_c
!
! Input:
!   t6(6,6)     -- real(rp): 1-turn transfer matrix.  RF assumed to be on.
!   abz_tunes(3) -- real(rp): a-mode and b-mode tunes.  Used to order eigensystem.
!
! Output:
!   B1(6,6)     -- real(rp): Beta matrix associated with a-mode.
!   B2(6,6)     -- real(rp): Beta matrix associated with b-mode.
!   B3(6,6)     -- real(rp): Beta matrix associated with c-mode.
!-

subroutine t6_to_B123(t6, abz_tunes, B1, B2, B3)

real(rp) t6(6,6)
real(rp) abz_tunes(3)
real(rp) B1(6,6), B2(6,6), B3(6,6)
logical err_flag

real(rp) mat_tunes(3)
real(rp) T1(6,6), T2(6,6), T3(6,6)
real(rp) eval_r(6), eval_i(6)
real(rp) evec_r(6,6), evec_i(6,6)
integer i

character(*), parameter :: r_name = 't6_to_b123'

!

T1 = 0.0d0
T2 = 0.0d0
T3 = 0.0d0

T1(1,2) = 1.0d0
T1(2,1) = 1.0d0
T2(3,4) = 1.0d0
T2(4,3) = 1.0d0
T3(5,6) = 1.0d0
T3(6,5) = 1.0d0

call eigen_decomp_6mat(t6, eval_r, eval_i, evec_r, evec_i, err_flag, mat_tunes)
if (err_flag) then
  call out_io (s_error$, r_name, "Error received from eigen_decomp_6mat.")
  B1 = 0.0d0
  B2 = 0.0d0
  B3 = 0.0d0
  return
endif

call order_evecs_by_tune(evec_r, evec_i, eval_r, eval_i, mat_tunes, abz_tunes, err_flag)
if (err_flag) then
  call out_io (s_error$, r_name, "Error received from order_evecs_by_tune.")
  B1 = 0.0d0
  B2 = 0.0d0
  B3 = 0.0d0
  return
endif

call normalize_evecs(evec_r, evec_i)

B1 = matmul(evec_r, matmul(T1, transpose(evec_r))) - matmul(evec_i, matmul(T1, transpose(evec_i))) 
B2 = matmul(evec_r, matmul(T2, transpose(evec_r))) - matmul(evec_i, matmul(T2, transpose(evec_i))) 
B3 = matmul(evec_r, matmul(T3, transpose(evec_r))) - matmul(evec_i, matmul(T3, transpose(evec_i))) 

end subroutine t6_to_B123

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine normal_mode3_calc (mat, tune, B, HV, above_transition)
!
! Does an Eigen decomposition of the 1-turn transfer matrix (mat) and generates
! B, V, H.
!
! If the above_transition argument is present and false, then the 3rd (z) mode is assumed
! to have a positive slip factor (z-mode rotates counter clockwise in phase space).
! Default is True ==> z-mode has a negative slip factor so the mode rotates clock-wise in phase space.
!
! Input:
!  mat(6,6)            -- real(rp): 1-turn transfer matrix
!  above_transition    -- logical, optional:  If present and false, then z-mode assumes positive slip factor.
!                                             Else negative slip factor assumed.
!
! Output:
!  tune(3)             -- real(rp): Tunes of the 3 normal modes (radians)
!  B(6,6)              -- real(rp): B is block diagonal and related to the normal mode Twiss parameters.
!  HV(6,6)             -- real(rp): Transforms from normal mode coordinates to canonical coordinates: x = H.V.a
!
!-

subroutine normal_mode3_calc (t6, tunes, B, HV, above_transition)

real(rp) t6(6,6)
real(rp) tunes(3)
real(rp) B(6,6)
real(rp) HV(6,6)
logical, optional :: above_transition

real(rp) N(6,6)
real(rp) V(6,6)
real(rp) H(6,6)
logical err_flag

character(*), parameter :: r_name = 'normal_mode3_calc'

!

call make_N(t6, N, err_flag, tunes_out=tunes)
if ( err_flag ) then
  call out_io (s_error$, r_name, "Error received from make_N.")
  tunes = 0.0d0
  B = 0.0d0
  HV = 0.0d0
  return
endif

if (.not. logic_option(.false., above_transition)) tunes(3) = tunes(3) + twopi

call make_HVBP (N, B, V, H)
HV = matmul(H, V)

!  HV=mat_symp_conj(HV) 
B = mat_symp_conj(B)  !for legacy compatability

end subroutine normal_mode3_calc

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine make_HVBP(N, B, V, H, Vbar, Hbar)
!
! Parameterizes the eigen-decomposition of the 6x6 transfer matrix into HVBP as defined in:
! "From the beam-envelop matrix to synchrotron-radiation integrals" by Ohmi, Hirata, and Oide.
!
! M = N.U.Inverse[N] where U is block diagonal and the blocks are 2x2 rotation matrices.
! N = H.V.B.P
! P has the same free parameters as B
! B "Twiss matrix" has 6 free parameters (Twiss alphas and betas)
! B blocks have the form /     sqrt(beta)         0       \
!                        \ -alpha/sqrt(beta) 1/sqrt(beta) /
! V "Teng matrix" has 4 free parameters (xy, xpy, ypx, and pxpy coupling)
! H "Dispersion matrix" has 8 free parameters (xz, xpz, pxz, pxpz, yz, ypz, pyz, pypz coupling)
! 
!
! Input:
!  N(6,6)              -- real(rp): Matrix of eigenvectors prepared by make_N
!
! Output:
!  B(6,6)              -- real(rp): Block diagonal matrix of Twiss parameters
!  V(6,6)              -- real(rp): horizontal-vertical coupling information
!  H(6,6)              -- real(rp): horizontal-longitudinal and vertical-longitudinal coupling information
!  Vbar(6,6)           -- real(rp), optional: mat_symp_conj(B).V.B
!  Hbar(6,6)           -- real(rp), optional: mat_symp_conj(B).H.B
!
!-

subroutine make_HVBP (N, B, V, H, Vbar, Hbar)

real(rp) N(6,6)
real(rp) B(6,6)
real(rp) R(6,6)
real(rp) H(6,6)
real(rp), optional :: Vbar(6,6)
real(rp), optional :: Hbar(6,6)

integer i

! Note: the variables are named here to according to the convention in the above mentioned paper.
real(rp) V(6,6)
real(rp) a, ax, ay
real(rp) BcPc(2, 2)
real(rp) Hx(2, 2)
real(rp) Hy(2, 2)
real(rp) VBP(6,6)
real(rp) mu
real(rp) BbPb(2, 2)
real(rp) BaPa(2, 2)
real(rp) V2(2, 2)
real(rp) BP(6,6)
real(rp) cospa, sinpa
real(rp) cospb, sinpb
real(rp) cospc, sinpc
real(rp) Pinv(6,6)

a = sqrt(abs(determinant(N(5:6, 5:6))))
BcPc = N(5:6, 5:6) / a
Hx = matmul(N(1:2, 5:6), mat_symp_conj(BcPc))
Hy = matmul(N(3:4, 5:6), mat_symp_conj(BcPc))
ax = determinant(Hx)/(1.0d0+a)  !shorthand
ay = determinant(Hy)/(1.0d0+a)  !shorthand

H(1:2, 1:2) = (1.0d0-ax)*I2
H(3:4, 3:4) = (1.0d0-ay)*I2
H(5:6, 5:6) = a*I2
H(1:2, 5:6) = Hx
H(3:4, 5:6) = Hy
H(5:6, 1:2) = -1.0d0*mat_symp_conj(Hx)
H(5:6, 3:4) = -1.0d0*mat_symp_conj(Hy)
H(1:2, 3:4) = -1.0d0 * matmul(Hx, mat_symp_conj(Hy)) / (1.0d0 + a)
H(3:4, 1:2) = -1.0d0 * matmul(Hy, mat_symp_conj(Hx)) / (1.0d0 + a)

VBP = matmul(mat_symp_conj(H), N)

mu = sqrt(abs(determinant(VBP(1:2, 1:2))))
BaPa = VBP(1:2, 1:2)/mu
BbPb = VBP(3:4, 3:4)/mu
V2 = matmul(VBP(1:2, 3:4), mat_symp_conj(BbPb))

V = 0.0d0
V(1:2, 1:2) = mu*I2
V(3:4, 3:4) = mu*I2
V(5:6, 5:6) = I2
V(1:2, 3:4) = V2
V(3:4, 1:2) = -1.0d0*mat_symp_conj(V2)

BP = 0.0d0
BP(1:2, 1:2) = BaPa
BP(3:4, 3:4) = BbPb
BP(5:6, 5:6) = BcPc

!- The following convention for P, puts B (the Twiss matrix) into the form where the upper right element is zero.
cospa = 1.0d0 / sqrt(1.0d0 + (BP(1,2)/BP(1,1))**2)
sinpa = 1.0d0 * BP(1,2) / BP(1,1) * cospa
cospb = 1.0d0 / sqrt(1.0d0 + (BP(3,4)/BP(3,3))**2)
sinpb = 1.0d0 * BP(3,4) / BP(3,3) * cospb
cospc = 1.0d0 / sqrt(1.0d0 + (BP(5,6)/BP(5,5))**2)
sinpc = 1.0d0 * BP(5,6) / BP(5,5) * cospc
Pinv = 0.0d0
Pinv(1,1) = cospa
Pinv(2,2) = cospa
Pinv(1,2) = -1.0d0 * sinpa
Pinv(2,1) = sinpa
Pinv(3,3) = cospb
Pinv(4,4) = cospb
Pinv(3,4) = -1.0d0 * sinpb
Pinv(4,3) = sinpb
Pinv(5,5) = cospc
Pinv(6,6) = cospc
Pinv(5,6) = -1.0d0 * sinpc
Pinv(6,5) = sinpc

B = matmul(BP, Pinv)

if ( present(Vbar) ) Vbar = matmul(mat_symp_conj(B), matmul(V, B))
if ( present(Hbar) ) Hbar = matmul(mat_symp_conj(B), matmul(H, B))

end subroutine make_HVBP

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine xyz_to_action(ring, loc, X, J, err_flag)
!
! Given the canonical phase space coordinates X of a particle, this returns
! a vector from which Ja, Jb, Jc can be easily extracted.
!
! The J vector looks like:
! J = (sqrt(2Ja)cos(phia), -sqrt(2Ja)sin(phia), sqrt(2Jb)cos(phib), -sqrt(2Jb)sin(phib), sqrt(2Jc)cos(phic), -sqrt(2Jc)sin(phic))
!
! NOTE: Mind the negative sign on the sin terms when computing phase advances or calculating frequency spectra.
!
! J is obtained from:
! J = N_inv . X
! Where N_inv is from the Eigen decomposition of the 1-turn transfer matrix.
!
! The normal mode invariant actions can be obtained from J as, 
! Ja = (J(1)**2 + J(2)**2)/2.0d0
! Jb = (J(3)**2 + J(4)**2)/2.0d0
! Jc = (J(5)**2 + J(6)**2)/2.0d0
!
! Input:
!  ring     -- lat_struct: lattice
!  loc      -- class(integer or real(rp)): location in ring: either element index or s coordinate.
!  X(1:6)   -- real(rp): canonical phase space coordinates of the particle
!
! Output:
!  J(1:6)   -- real(rp): Vector containing normal mode invariants and phases
!  err_flag -- logical: Set to true on error.  Often means Eigen decomposition failed.
!-

subroutine xyz_to_action(ring, loc, X, J, err_flag)

type(lat_struct) ring
class(*) :: loc
real(rp) X(6)
real(rp) J(6)
logical err_flag

real(rp) t6(6,6)
real(rp) N(6,6)
real(rp) abz_tunes(3)

integer ix_use

type(ele_struct) ele_at_s

character(*), parameter :: r_name = 'xyz_to_action'

!

abz_tunes(1) = ring%a%tune
abz_tunes(2) = ring%b%tune
abz_tunes(3) = ring%z%tune

select type(loc)
type is (integer)
  call transfer_matrix_calc (ring, t6, ix1=loc, one_turn=.true.)
type is (real(rp))
  ix_use = element_at_s(ring, loc, .false.)
  call transfer_matrix_calc (ring, t6, ix1=ix_use, one_turn=.true.)
  call twiss_and_track_at_s (ring, loc, ele_at_s)
  t6 = matmul(ele_at_s%mat6, matmul(mat_symp_conj(ring%ele(ix_use)%mat6), matmul(t6, matmul(ring%ele(ix_use)%mat6, mat_symp_conj(ele_at_s%mat6)))))
end select

if (all(abs(abz_tunes).gt. 0.0001)) then
  call make_N(t6, N, err_flag, abz_tunes)
else
  call make_N(t6, N, err_flag)
endif
if ( err_flag ) then
  call out_io (s_error$, r_name, "Error received from make_N.")
  J = 0.0d0
  return
endif

J = matmul(mat_symp_conj(N), X)
end subroutine xyz_to_action

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine action_to_xyz(ring, ix, J, X, err_flag)
!
! Given the normal mode invariants and phases J of a particle, returns the canonical coordinates.
!
! The J vector looks like:
! J = (sqrt(2Ja)cos(phia), -sqrt(2Ja)sin(phia), sqrt(2Jb)cos(phib), -sqrt(2Jb)sin(phib), sqrt(2Jc)cos(phic), -sqrt(2Jc)sin(phic))
!
! X is obtained from:
! X = N . J
! Where N is from the Eigen decomposition of the 1-turn transfer matrix.
!
! Input:
!  ring     -- lat_struct: lattice
!  ix       -- integer: element index at which to calculate J
!  J(1:6)   -- real(rp): Vector containing normal mode invariants and phases
!
! Output:
!  X(1:6)   -- real(rp): canonical phase space coordinates of the particle
!  err_flag    -- logical: Set to true on error.  Often means Eigen decomposition failed.
!-
subroutine action_to_xyz(ring, ix, J, X, err_flag)

type(lat_struct) ring
integer ix
real(rp) J(1:6)
real(rp) X(1:6)
logical err_flag 

real(rp) t6(1:6, 1:6)
real(rp) N(1:6, 1:6)
real(rp) Ninv(1:6, 1:6)
real(rp) abz_tunes(3)
real(rp) gamma(3)
integer i

character(*), parameter :: r_name = 'action_to_xyz'

!

abz_tunes(1) = ring%a%tune
abz_tunes(2) = ring%b%tune
abz_tunes(3) = ring%z%tune

call transfer_matrix_calc (ring, t6, ix1=ix, one_turn=.true.)
call make_N(t6, N, err_flag, abz_tunes)
if ( err_flag ) then
  call out_io (s_error$, r_name, "Error received from make_N.")
  X = 0.0d0
  return
endif

X = matmul(N, J)

end subroutine action_to_xyz

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine eigen_decomp_6mat(mat, eval_r, eval_i, evec_r, evec_i, tunes, err_flag)
!
! Compute eigenvalues and eigenvectors of a real 6x6 matrix.
! The evals and evecs are in general complex.
!
! Input:
!   mat(6,6)     - real(rp):  6x6 real matrix.  Usually a transfer matrix or sigma matrix.
!
! Output:
!   eval_r(6)    - real(rp):  real part of eigenvalues.
!   eval_i(6)    - real(rp):  complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  real part of eigenvectors arranged down columns.
!   evec_i(6,6)  - real(rp):  complex part of eigenvectors arranged down columns.
!   err_flag     - logical, optional: set to true if an error has occured.
!   tunes(3)     - real(rp):  Mode tunes, in radians.
!-

subroutine eigen_decomp_6mat(mat, eval_r, eval_i, evec_r, evec_i, err_flag, tunes)

use la_precision, only: wp => dp
use f95_lapack

real(rp) mat(6,6)
real(rp) A(6,6)
real(rp) VR(6,6)
real(rp) eval_r(6), eval_i(6)
real(rp) evec_r(6,6), evec_i(6,6)
logical err_flag
real(rp), optional :: tunes(3)

integer i_error
integer pair1, pair2, pair3, pairIndexes(6)
integer i
complex(rp) evec(6,6)
complex(rp) check_mat(6,6)

character(*), parameter :: r_name = 'eigen_decomp_6mat'

!call mat_eigen (mat, eval_r, eval_i, evec_r, evec_i, err_flag)
!evec_r = transpose(evec_r)
!evec_i = transpose(evec_i)

err_flag = .true.

A = mat  !LA_GEEV destroys the contents of its first argument.
CALL la_geev(A, eval_r, eval_i, VR=VR, INFO=i_error)
if ( i_error /= 0 ) THEN
  call out_io (s_fatal$, r_name, "la_geev returned error: \i0\ ", i_error)
  if (global_com%exit_on_error) call err_exit
  eval_r = 0.0d0
  eval_i = 0.0d0
  evec_r = 0.0d0
  evec_i = 0.0d0
  return
endif

evec_r(:, 1) = VR(:, 1)
evec_r(:, 2) = VR(:, 1)
evec_r(:, 3) = VR(:, 3)
evec_r(:, 4) = VR(:, 3)
evec_r(:, 5) = VR(:, 5)
evec_r(:, 6) = VR(:, 5)
evec_i(:, 1) = VR(:, 2)
evec_i(:, 2) = -VR(:, 2)
evec_i(:, 3) = VR(:, 4)
evec_i(:, 4) = -VR(:, 4)
evec_i(:, 5) = VR(:, 6)
evec_i(:, 6) = -VR(:, 6)

evec = cmplx(evec_r, evec_i, rp)
check_mat = matmul(transpose(conjg(evec)), matmul(S, evec))

if ( aimag(check_mat(1,1)) > 0.0 ) then
  evec_i(:, 1) = -evec_i(:, 1)
  evec_i(:, 2) = -evec_i(:, 2)
  eval_i(1) = -eval_i(1)
  eval_i(2) = -eval_i(2)
endif
if ( aimag(check_mat(3,3)) > 0.0 ) then
  evec_i(:, 3) = -evec_i(:, 3)
  evec_i(:, 4) = -evec_i(:, 4)
  eval_i(3) = -eval_i(3)
  eval_i(4) = -eval_i(4)
endif
if ( aimag(check_mat(5,5)) > 0.0 ) then
  evec_i(:, 5) = -evec_i(:, 5)
  evec_i(:, 6) = -evec_i(:, 6)
  eval_i(5) = -eval_i(5)
  eval_i(6) = -eval_i(6)
endif

if (present(tunes)) then
  tunes(1) = MyTan(-eval_i(1), eval_r(1))
  tunes(2) = MyTan(-eval_i(3), eval_r(3))
  tunes(3) = MyTan(-eval_i(5), eval_r(5))
endif

err_flag = .false.

!---------------------------------------------------------
contains

function MyTan(y, x) result(arg)
  !For a complex number x+iy graphed on an xhat, yhat plane, this routine returns the angle
  !between (x, y) and the +x axis, measured counter-clockwise.  There is a branch cut along +x.
  !This routine returns a number between 0 and 2pi.
  real(rp) x, y, arg

  if (y .ge. 0) then
    arg = atan2(y, x)
  else
    arg = atan2(y, x) + 2.0d0*pi
  endif
end function MyTan

end subroutine eigen_decomp_6mat

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine order_evecs_by_N_similarity(evec_r, evec_i, eval_r, eval_i, mat_tunes, Nmat)
! 
! This subroutine orderes the eigensystem such that Nmat.mat_symp_conj(N) is closest
! to the identity.  Nmat is supplied externally.
!
! Input:
!   eval_r(6)    - real(rp):  real part of eigenvalues.
!   eval_i(6)    - real(rp):  complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  real part of eigenvectors arranged down columns.
!   evec_i(6,6)  - real(rp):  complex part of eigenvectors arranged down columns.
!   mat_tunes(3) - real(rp):  Three normal mode tunes, in radians.
!   Nmat(6,6)    - real(rp):  Normalized, real eigen matrix from make_N.
!
! Output:
!   eval_r(6)    - real(rp):  Ordered real part of eigenvalues.
!   eval_i(6)    - real(rp):  Ordered complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  Ordered real part of eigenvectors arranged down columns.
!   evec_i(6,6)  - real(rp):  Ordered complex part of eigenvectors arranged down columns.
!   mat_tunes(3) - real(rp):  Ordered normal mode tunes, in radians.
!-

subroutine order_evecs_by_N_similarity(evec_r, evec_i, eval_r, eval_i, mat_tunes, Nmat)

real(rp) evec_r(6,6), evec_i(6,6)
real(rp) eval_r(6), eval_i(6)
real(rp) mat_tunes(3)
real(rp) Nmat(6,6)

real(rp) evec_r_local(6,6), evec_i_local(6,6)
real(rp) Nlocal(6,6)
integer pair1, pair2, pair3
integer pairIndexes(6)
integer iterations(6, 3)
real(rp) check_mat(6,6)
real(rp) checks(6)
integer i, best

iterations(1, :) = [1, 2, 3]
iterations(2, :) = [1, 3, 2]
iterations(3, :) = [2, 1, 3]
iterations(4, :) = [2, 3, 1]
iterations(5, :) = [3, 1, 2]
iterations(6, :) = [3, 2, 1]

do i=1, 6
  pair1 = iterations(i, 1)
  pair2 = iterations(i, 2)
  pair3 = iterations(i, 3)

  pairIndexes = [ 2*pair1-1, 2*pair1, 2*pair2-1, 2*pair2, 2*pair3-1, 2*pair3 ]
  evec_r_local = evec_r(:, pairIndexes)
  evec_i_local = evec_i(:, pairIndexes)
  call normalize_evecs(evec_r_local, evec_i_local)

  Nlocal = matmul(evec_r_local, Qr) - matmul(evec_i_local, Qi)

  check_mat = matmul(mat_symp_conj(Nlocal), Nmat)

  checks(i) = abs(trace(abs(check_mat)) - 6.0d0)
enddo
best= minloc(checks, 1)
pair1 = iterations(best, 1)
pair2 = iterations(best, 2)
pair3 = iterations(best, 3)
pairIndexes = [ 2*pair1-1, 2*pair1, 2*pair2-1, 2*pair2, 2*pair3-1, 2*pair3 ]

evec_r = evec_r(:, pairIndexes)
evec_i = evec_i(:, pairIndexes)
eval_r = eval_r(pairIndexes)
eval_i = eval_i(pairIndexes)
mat_tunes = mat_tunes([pair1, pair2, pair3])

!--------------------------------
contains

function trace(mat)
  real(rp) trace
  real(rp) mat(:, :)
  integer i
  trace = 0.0d0
  do i=1, size(mat(:, 1))
    trace = trace + mat(i, i) 
  enddo
end function trace

end subroutine order_evecs_by_N_similarity

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine order_evecs_by_plane_dominance(evec_r, evec_i, eval_r, eval_i, mat_tunes)
! 
! This subroutine orderes the eigensystem according to which modes dominate the horizontal, 
! vertical, and longitudinal planes.  This subroutine works well in machines
! that are not strongly coupled.  In machines with strong coupling, where the relation
! between the three eigenmodes a, b, c and the three lab coordinates x, y, z can change
! through the machine, this subroutine will not provide consistent ordering.
!
! Input:
!   eval_r(6)    - real(rp):  real part of eigenvalues.
!   eval_i(6)    - real(rp):  complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  real part of eigenvectors arranged down columns.
!   evec_i(6,6)  - real(rp):  complex part of eigenvectors arranged down columns.
!   mat_tunes(3) - real(rp), optional:  Three normal mode tunes, in radians.
!
! Output:
!   eval_r(6)    - real(rp):  Ordered real part of eigenvalues.
!   eval_i(6)    - real(rp):  Ordered complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  Ordered real part of eigenvectors.
!   evec_i(6,6)  - real(rp):  Ordered complex part of eigenvectors.
!   mat_tunes(3) - real(rp), optional:  Reordered same as evecs.
!-

subroutine order_evecs_by_plane_dominance(evec_r, evec_i, eval_r, eval_i, mat_tunes)

real(rp) evec_r(6,6), evec_i(6,6)
real(rp) eval_r(6), eval_i(6)
real(rp), optional :: mat_tunes(3)
real(rp) vec(3)
integer pair1, pair2, pair3
integer pairindexes(6)

pair1 = maxloc([abs(evec_r(1,1)), abs(evec_r(1,3)), abs(evec_r(1,5))], 1)

vec = [abs(evec_r(3,1)), abs(evec_r(3,3)), abs(evec_r(3,5))]
vec(pair1) = -1
pair2 = maxloc(vec, 1)

pair3 = 6 - pair1 - pair2

pairIndexes = [2*pair1-1, 2*pair1, 2*pair2-1, 2*pair2, 2*pair3-1, 2*pair3]

evec_r = evec_r(:, pairIndexes)
evec_i = evec_i(:, pairIndexes)
eval_r = eval_r(pairIndexes)
eval_i = eval_i(pairIndexes)

if(present(mat_tunes)) then
  mat_tunes = mat_tunes([pair1, pair2, pair3])
endif

end subroutine order_evecs_by_plane_dominance

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine order_evecs_by_tune(evec_r, evec_i, eval_r, eval_i, mat_tunes, abz_tunes, err_flag)
! 
! This subroutine orders the eigensystem by matching the tunes of the eigensystem to 
! externally supplied tunes abz_tunes.  abz_tunes is in radians.
!
! Input:
!   eval_r(6)    - real(rp):  real part of eigenvalues.
!   eval_i(6)    - real(rp):  complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  real part of eigenvectors arranged down columns.
!   evec_i(6,6)  - real(rp):  complex part of eigenvectors arranged down columns.
!   mat_tunes(3) - real(rp):  Three normal mode tunes, in radians.
!   abz_tunes(3)  - real(rp):  Tunes to order eigensystem by.
!
! Output:
!   eval_r(6)    - real(rp):  Ordered real part of eigenvalues.
!   eval_i(6)    - real(rp):  Ordered complex part of eigenvalues.
!   evec_r(6,6)  - real(rp):  Ordered real part of eigenvectors.
!   evec_i(6,6)  - real(rp):  Ordered complex part of eigenvectors.
!   err_flag     - logical, optional:  Set to true if an error occured.
!-

subroutine order_evecs_by_tune (evec_r, evec_i, eval_r, eval_i, mat_tunes, abz_tunes, err_flag)

real(rp) evec_r(6,6), evec_i(6,6)
real(rp) eval_r(6), eval_i(6)
real(rp) mat_tunes(3)
logical err_flag

real(rp) abz_tunes(3)
real(rp) val(6)
integer pairindexes(6)

character(*), parameter :: r_name = 'order_evecs_by_tune'

! Order eigenvector pairs

err_flag = .true.

if ( any(abs(mat_tunes(1:3)) .lt. 0.0001) ) then
  call out_io (s_fatal$, r_name, "mat_tunes is not fully populated.  Printing mat_tunes.")
  write(*, '(3f14.5)') mat_tunes(1:3)
  if (global_com%exit_on_error) call err_exit
  return
endif
if ( any(abs(abz_tunes(1:3)) .lt. 0.0001) ) then
  call out_io (s_fatal$, r_name, "abz_tunes is not fully populated.  Printing abz_tunes.")
  write(*, '(3f14.5)') abz_tunes(1:3)
  if (global_com%exit_on_error) call err_exit
  return
endif

val(1) = ( (mat_tunes(1)-abz_tunes(1)) / (mat_tunes(1)+abz_tunes(1)) )**2 + &
         ( (mat_tunes(2)-abz_tunes(2)) / (mat_tunes(2)+abz_tunes(2)) )**2 + &
         ( (mat_tunes(3)-abz_tunes(3)) / (mat_tunes(3)+abz_tunes(3)) )**2
val(2) = ( (mat_tunes(1)-abz_tunes(1)) / (mat_tunes(1)+abz_tunes(1)) )**2 + &
         ( (mat_tunes(2)-abz_tunes(3)) / (mat_tunes(2)+abz_tunes(3)) )**2 + &
         ( (mat_tunes(3)-abz_tunes(2)) / (mat_tunes(3)+abz_tunes(2)) )**2
val(3) = ( (mat_tunes(1)-abz_tunes(2)) / (mat_tunes(1)+abz_tunes(2)) )**2 + &
         ( (mat_tunes(2)-abz_tunes(1)) / (mat_tunes(2)+abz_tunes(1)) )**2 + &
         ( (mat_tunes(3)-abz_tunes(3)) / (mat_tunes(3)+abz_tunes(3)) )**2
val(4) = ( (mat_tunes(1)-abz_tunes(2)) / (mat_tunes(1)+abz_tunes(2)) )**2 + &
         ( (mat_tunes(2)-abz_tunes(3)) / (mat_tunes(2)+abz_tunes(3)) )**2 + &
         ( (mat_tunes(3)-abz_tunes(1)) / (mat_tunes(3)+abz_tunes(1)) )**2
val(5) = ( (mat_tunes(1)-abz_tunes(3)) / (mat_tunes(1)+abz_tunes(3)) )**2 + &
         ( (mat_tunes(2)-abz_tunes(1)) / (mat_tunes(2)+abz_tunes(1)) )**2 + &
         ( (mat_tunes(3)-abz_tunes(2)) / (mat_tunes(3)+abz_tunes(2)) )**2
val(6) = ( (mat_tunes(1)-abz_tunes(3)) / (mat_tunes(1)+abz_tunes(3)) )**2 + &
         ( (mat_tunes(2)-abz_tunes(2)) / (mat_tunes(2)+abz_tunes(2)) )**2 + &
         ( (mat_tunes(3)-abz_tunes(1)) / (mat_tunes(3)+abz_tunes(1)) )**2

if (minval(val, 1) .gt. 0.1) then
  call out_io (s_fatal$, r_name, "Unable to match mat_tunes.  Printing mat_tunes and abz_tunes.")
  write(*, '(3f14.5)') mat_tunes(1:3)
  write(*, '(3f14.5)') abz_tunes(1:3)
  return
endif

select case(minloc(val, 1))
  case(1)
    pairindexes = [1, 2, 3, 4, 5, 6]
    mat_tunes = mat_tunes([1, 2, 3])
  case(2)
    pairindexes = [1, 2, 5, 6, 3, 4]
    mat_tunes = mat_tunes([1, 3, 2])
  case(3)
    pairindexes = [3, 4, 1, 2, 5, 6]
    mat_tunes = mat_tunes([2, 1, 3])
  case(4)
    pairindexes = [3, 4, 5, 6, 1, 2]
    mat_tunes = mat_tunes([2, 3, 1])
  case(5)
    pairindexes = [5, 6, 1, 2, 3, 4]
    mat_tunes = mat_tunes([3, 1, 2])
  case(6)
    pairindexes = [5, 6, 3, 4, 1, 2]
    mat_tunes = mat_tunes([3, 2, 1])
end select

evec_r = evec_r(:, pairindexes)
evec_i = evec_i(:, pairIndexes)
eval_r = eval_r(pairIndexes)
eval_i = eval_i(pairIndexes)

err_flag = .false.

end subroutine order_evecs_by_tune

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine make_N(t6, N, err_flag, abz_tunes, tunes_out)
!
! Given a 1-turn transfer matrix, this returns N and its inverse Ninv.
! N converts between normal invarients and phases and canonical coordinates:
! X = N.J
!
! N is obtained from the Eigen decomposition of the 1-turn transfer matrix.
! It is obtained by applying certain normalizations to the matrix of Eigen vectors, then making
! the result real using Q.
!
! If abz_tunes is present, then the eigensystem is ordered by matching the tunes.
! If abz_tunes is not present, then the eigensystem is ordered by plane dominance.
!
! It is assumed that the synchrotron tune is less than pi.
!
! Input:
!  t6(6,6)             -- real(rp): 1-turn transfer matrix
!  abz_tunes(3)        -- real(rp), optional: a-mode is abz_tunes(1), b-mode is abz_tunes(2), synch tune is abz_tunes(3)
!
! Output:
!  N(6,6)              -- real(rp): X = N.J
!  err_flag            -- logical: Set to true on error.  Often means Eigen decomposition failed.
!  tunes_out(3)        -- real(rp), optional: Fractional tune (in radians) of the 3 normal modes of t6.
!-

subroutine make_N(t6, N, err_flag, abz_tunes, tunes_out)

real(rp) t6(6,6)
real(rp) N(6,6)
logical err_flag
real(rp), optional :: abz_tunes(3)
real(rp), optional :: tunes_out(3)

real(rp) eval_r(6), eval_i(6)
real(rp) evec_r(6,6), evec_i(6,6)
real(rp) mat_tunes(3)
real(rp) mat(6,6)
integer i

character(*), parameter :: r_name = 'make_N'

!

call eigen_decomp_6mat(t6, eval_r, eval_i, evec_r, evec_i, err_flag, mat_tunes)
if (err_flag) then
  call out_io (s_error$, r_name, "Error received from eigen_decomp_6mat.")
  N = 0.0d0
  if (present(abz_tunes)) abz_tunes = 0.0d0
  if (present(tunes_out)) tunes_out = 0.0d0
  return
endif

if ( present(abz_tunes) ) then
  call order_evecs_by_tune(evec_r, evec_i, eval_r, eval_i, mat_tunes, abz_tunes, err_flag)
  if (err_flag) then
    call out_io (s_error$, r_name, "Error received from order_evecs_by_tune.")
    N = 0.0d0
    if (present(abz_tunes)) abz_tunes = 0.0d0
    if (present(tunes_out)) tunes_out = 0.0d0
    return
  endif
else
  call order_evecs_by_plane_dominance(evec_r, evec_i, eval_r, eval_i, mat_tunes)
endif

! if (abs(mat_tunes(3)) .gt. pi) then  !assume synchrotron tune less than pi
!   mat_tunes(3) = mat_tunes(3) - twopi
! endif

call normalize_evecs(evec_r, evec_i)

N = matmul(evec_r, Qr) - matmul(evec_i, Qi)

if (present(tunes_out)) then
  tunes_out = mat_tunes
endif
end subroutine make_N

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine get_emit_from_sigma_mat(sigma_mat, normal, Nmat, err_flag)
!
! Given a beam envelop sigma matrix sigma_mat, this returns the 3 normal mode
! emittances.
!
! The normal mode emittance of the sigma matrix are the eigenvalues of
! sigma_mat . S
!
! If Nmat is present, then the modes are ordered such that the eigensystem most
! closely resembles Nmat.  If Nmat is not present, then the modes are ordered
! according to which plane they dominate.
!
!     / 0  1  0  0  0  0 \
!     |-1  0  0  0  0  0 |
! S = | 0  0  0  1  0  0 |
!     | 0  0 -1  0  0  0 |
!     | 0  0  0  0  0  1 |
!     \ 0  0  0  0 -1  0 /
!
! Input:
!  sigma_mat(6,6)   -- real(rp): beam envelop sigma matrix
!  Nmat(6,6)        -- real(rp): If present, then the emittanced will be ordered such that
!                                the eigensystem most closely resembles Nmat.
! Output:
!  normal(3)        -- real(rp): normal mode emittances
!  err_flag         -- logical: Set to true if something went wrong.  Otherwise set to false.
!-

subroutine get_emit_from_sigma_mat(sigma_mat, normal, Nmat, err_flag)

real(rp) sigma_mat(6,6)
real(rp) normal(3)
logical err_flag
real(rp), optional :: Nmat(6,6)

real(rp) tunes(3)
real(rp) eval_r(6), eval_i(6)
real(rp) evec_r(6,6), evec_i(6,6)
real(rp) sigmas(6,6)

integer i

character(*), parameter :: r_name = 'get_emit_from_sigma_mat'

!
sigmaS = matmul(sigma_mat, S)

call eigen_decomp_6mat(sigmas, eval_r, eval_i, evec_r, evec_i, err_flag, tunes)
if (err_flag) then
  call out_io (s_error$, r_name, "Error received from eigen_decomp_6mat.")
  normal = 0.0d0
  return
endif

if (present(Nmat)) then
  call order_evecs_by_N_similarity(evec_r, evec_i, eval_r, eval_i, tunes, Nmat)
else
  call order_evecs_by_plane_dominance(evec_r, evec_i, eval_r, eval_i, tunes)
endif

normal(1) = abs(eval_i(1))
normal(2) = abs(eval_i(3))
normal(3) = abs(eval_i(5))

end subroutine get_emit_from_sigma_mat

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine beam_tilts(S, angle_xy, angle_xz, angle_yz, angle_xpz, angle_ypz)
!
! Given a 6x6 matrix of second-order moments, this routine returns
! the beam tilts.
!
! angle_xy is obtained from the projection of the beam envelop into the
! xy plane.  The angle is that between the major axis of the projected
! beam envelope and the +x axis.  Positive angles are measured towards the
! +y axis.
!
! angle_xz is obtained from the projection of the beam envelop into the
! xy plane.  The angle is that between the major axis of the projected beam envelope
! and the +z axis.  Positive angles are measured towards the +x axis.
!
! angle_yz is obtained from the projection of the beam envelop into the
! yz plane.  The angle is that between the major axis of the projected beam envelope
! and the +z axis.  Positive angles are measured towards the +y axis.
!
! Input:
!   S(6,6)              -- real(rp): matrix of second order moments of beam envelope
!
! Output:
!   angle_xy            -- real(rp): transverse tilt of beam envelope
!   angle_xz            -- real(rp): horizontal crabbing of beam envelope
!   angle_yz            -- real(rp): vertical crabbing of beam envelope
!   angle_xpz           -- real(rp): x-pz coupling
!   angle_ypz           -- real(rp): y-pz coupling
!-

subroutine beam_tilts(S, angle_xy, angle_xz, angle_yz, angle_xpz, angle_ypz)

real(rp) S(6,6)
real(rp) angle_xy, angle_xz, angle_yz
real(rp) angle_xpz, angle_ypz

angle_xy  = 0.5_rp * atan2( 2.0d0*S(1, 3), S(1,1)-S(3,3) )
angle_xz  = 0.5_rp * atan2( 2.0d0*S(1, 5), S(5,5)-S(1,1) )
angle_yz  = 0.5_rp * atan2( 2.0d0*S(3, 5), S(5,5)-S(3,3) )
angle_xpz = 0.5_rp * atan2( 2.0d0*S(1, 6), S(6,6)-S(1,1) )
angle_ypz = 0.5_rp * atan2( 2.0d0*S(3, 6), S(6,6)-S(3,3) )

end subroutine beam_tilts

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine make_smat_from_abc(t6, mode, sigma_mat, err_flag, Nout)
!
! Given the 1-turn transfer matrix and a normal_modes_struct containing the normal mode
! emittances, this routine returns the beam envelop sigma matrix.
!
! sigma_mat = N.D.transpose(N)
! equivalent to: sigma_mat.S = N.D.mat_symp_conj(N)
!
! One way to populate mode%a%tune and mode%b%tune:
!   mode%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi, twopi)
!   mode%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi, twopi)
!
! Input:
!  t6(6,6)          -- real(rp): 1-turn transfer matrix
!  mode             -- normal_modes_struct: normal mode emittances
!      %a%emittance -- real(rp): a-mode emittance
!      %b%emittance -- real(rp): b-mode emittance
!      %z%emittance -- real(rp): z-mode emittance
!      %a%tune      -- real(rp): a-mode tune.  Used to associate emittances with the proper mode.
!      %b%tune      -- real(rp): b-mode tune.  Used to associate emittances with the proper mode.
!      %z%tune      -- real(rp): z-mode tune.  Used to associate emittances with the proper mode.
!
! Output:
!  sigma_mat(6,6)   -- real(rp): beam envelop sigma matrix
!  err_flag         -- logical:  set to true if something goes wrong.  Usually means Eigen decomposition of the 1-turn matrix failed.
!  Nout(6,6)        -- real(rp), optional: Contains the normalized eigenvectors that were used to make the sigma matrix.
!-

subroutine make_smat_from_abc(t6, mode, sigma_mat, err_flag, Nout)

real(rp) t6(6,6)
type(normal_modes_struct) mode
real(rp) sigma_mat(6,6)
logical err_flag
real(rp), optional :: Nout(6,6)

real(rp) N(6,6)
real(rp) D(6,6)
real(rp) abz_tunes(3)

integer i

character(*), parameter :: r_name = 'make_smat_from_abc'

!

abz_tunes(1) = mode%a%tune
abz_tunes(2) = mode%b%tune
abz_tunes(3) = mode%z%tune

call make_N(t6, N, err_flag, abz_tunes)
if (err_flag) then
  call out_io (s_error$, r_name, "Error received from make_N.")
  sigma_mat = 0.0d0
  if (present(Nout)) Nout = 0.0d0
  return
endif

if (present(Nout)) Nout = N

! make_N takes the normal mode tunes and sorts N such that the first two columns are associated with the a-mode, 
! second two with b-mode, and last two with z-mode.

D=0.0d0
D(1,1) = mode%a%emittance
D(2,2) = mode%a%emittance
D(3,3) = mode%b%emittance
D(4,4) = mode%b%emittance
D(5,5) = mode%z%emittance
D(6,6) = mode%z%emittance

sigma_mat = matmul( matmul(N, D), transpose(N) )

end subroutine make_smat_from_abc

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine normalize_evecs(evec_r, evec_i)
!
! Normalizes eigenvectors such that transpose(E).S.E = iS, where E = evec_r + i evec_i
!
! Input:
!   evec_r(6,6)  - real(rp):  real part of eigenvectors arranged down columns.
!   evec_i(6,6)  - real(rp):  complex part of eigenvectors arranged down columns.
! Output:
!   evec_r(6,6)  - real(rp):  Eigensystem normalized to be symplectic.
!   evec_i(6,6)  - real(rp):  Eigensystem normalized to be symplectic.
!-

subroutine normalize_evecs(evec_r, evec_i)

real(rp) evec_r(6,6), evec_i(6,6)
real(rp) norm1, norm2, norm3
complex(rp) evec(6,6)
complex(rp) mat(6,6)

evec = cmplx(evec_r, evec_i, rp)
mat = matmul(transpose(conjg(evec)), matmul(S, evec))

norm1 = sqrt(abs(mat(1,1)))
norm2 = sqrt(abs(mat(3,3)))
norm3 = sqrt(abs(mat(5,5)))

evec_r(:, 1:2) = evec_r(:, 1:2) / norm1
evec_i(:, 1:2) = evec_i(:, 1:2) / norm1
evec_r(:, 3:4) = evec_r(:, 3:4) / norm2
evec_i(:, 3:4) = evec_i(:, 3:4) / norm2
evec_r(:, 5:6) = evec_r(:, 5:6) / norm3
evec_i(:, 5:6) = evec_i(:, 5:6) / norm3

end subroutine normalize_evecs

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine project_emit_to_xyz(ring, ix, mode, sigma_x, sigma_y, sigma_z)
!
! Obtains the projected x, y, and z beamsizes by building the sigma matrix
! from the normal mode emittances and 1-turn transfer matrix.
! These projectes beamsize are what would be seen by instrumentation.
!
! This method of projecting takes into account transverse and longitudinal coupling.
!
! This method of obtaining the projected beam sizes is from "Alternitive approach to general
! coupled linear optics" by Andrzej Wolski.
!
! The normal mode emittances used to generate a beam envelop sigma matrix from the 
! 1-turn transfer matrix.  The projected sizes are from the 1, 1 3, 3 and 5, 5 elements of
! the sigma matrix.
!
! Input:
!  ring             -- lat_struct: the storage ring
!  ix               -- integer: element at which to make the projection
!  mode             -- normal_modes_struct: normal mode emittances
!      %a%emittance -- real(rp): a-mode emittance
!      %b%emittance -- real(rp): b-mode emittance
!      %z%emittance -- real(rp): z-mode emittance
!      %a%tune      -- real(rp): a-mode tune.  Used to associate emittances with the proper mode.
!      %b%tune      -- real(rp): b-mode tune.  Used to associate emittances with the proper mode.
!
! Output:
!  sigma_x          -- real(rp): projected horizontal beamsize
!  sigma_y          -- real(rp): projected vertical beamsize
!  sigma_z          -- real(rp): projected longitudinal beamsize
!-

subroutine project_emit_to_xyz(ring, ix, mode, sigma_x, sigma_y, sigma_z)

type(lat_struct) ring
integer ix
type(normal_modes_struct) mode
real(rp) sigma_x, sigma_y, sigma_z
logical err_flag

real(rp) t6(6,6)
real(rp) sigma_mat(6,6)

character(*), parameter :: r_name = 'project_emit_to_xyz'

!

call transfer_matrix_calc (ring, t6, ix1=ix, one_turn=.true.)
call make_smat_from_abc(t6, mode, sigma_mat, err_flag)
if (err_flag) then
  call out_io (s_error$, r_name, "Error received from make_smat_from_abc.")
  sigma_x = 0.0d0
  sigma_y = 0.0d0
  sigma_z = 0.0d0
  return
endif

sigma_x = sqrt(sigma_mat(1, 1))
sigma_y = sqrt(sigma_mat(3, 3))
sigma_z = sqrt(sigma_mat(5, 5))

end subroutine project_emit_to_xyz

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine twiss3_propagate_all (lat, ix_branch)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
! Subroutine from original mode3_mod.
!
! Input:
!   lat       -- lat_struct: Lattice
!   ix_branch -- integer, optional :: Branch index. 0 = default.
!-

subroutine twiss3_propagate_all (lat, ix_branch)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer, optional :: ix_branch
integer i
logical err_flag

!

branch => lat%branch(integer_option(0, ix_branch))

do i = 1, branch%n_ele_track
  call twiss3_propagate1 (branch%ele(i-1), branch%ele(i), err_flag)
enddo

end subroutine twiss3_propagate_all

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine twiss3_propagate1 (ele1, ele2, err_flag)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
! Subroutine from original mode3_mod.
!-

subroutine twiss3_propagate1 (ele1, ele2, err_flag)

type mat2_struct
  real(rp) m(2, 2)
end type

type (ele_struct) ele1, ele2
type (mat2_struct) w(3)

real(rp) gamma(3), tv(6,6), w_inv(2, 2), radx

integer i, ik
logical err, err_flag

!

err_flag = .true.

if (.not. associated(ele2%mode3)) allocate(ele2%mode3)

tv = matmul (ele2%mat6, ele1%mode3%v)

do i = 1, 3
  ik = 2 * i - 1
  w(i)%m = tv(ik:ik+1, ik:ik+1)
  radx = determinant (w(i)%m)
  if (radx < 0) return
  gamma(i) = SQRT(radx)
  w(i)%m = w(i)%m / gamma(i)
  w_inv = mat_symp_conj(w(i)%m)
  ele2%mode3%v(1:6, ik:ik+1) = matmul(tv(1:6, ik:ik+1), w_inv)
enddo

ele2%mode3%x%eta = ele2%mode3%v(1, 6)
ele2%mode3%y%eta = ele2%mode3%v(3, 6)

ele2%mode3%x%etap = ele2%mode3%v(2, 6)
ele2%mode3%y%etap = ele2%mode3%v(4, 6)

call twiss1_propagate (ele1%mode3%a, w(1)%m, ele2%key, ele2%value(l$), ele2%mode3%a, err)
if (err) return
call twiss1_propagate (ele1%mode3%b, w(2)%m, ele2%key, ele2%value(l$), ele2%mode3%b, err)
if (err) return
call twiss1_propagate (ele1%mode3%c, w(3)%m, ele2%key, 0.0_rp,         ele2%mode3%c, err)
if (err) return

err_flag = .false.

end subroutine twiss3_propagate1 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine twiss3_from_twiss2 (ele)
!
! Routine to calculate the 3D Twiss parameters given the 2D transverse Twiss parameters and some 
! longitudinal parameters.
! Also see: twiss3_at_start
!
! Input:
!   ele     -- ele_struct: Lattice element at which the calculation is made.
!
! Output:
!   ele     -- ele_struct: Element
!-

subroutine twiss3_from_twiss2 (ele)

type (ele_struct) ele
real(rp) tune3(3), d_mat(6,6)

!

if (.not. associated(ele%mode3)) allocate(ele%mode3)

ele%mode3%a = ele%a
ele%mode3%b = ele%b
ele%mode3%c = ele%z

call mat_make_unit(ele%mode3%v)
call make_v_mats(ele, ele%mode3%v(1:4,1:4))

ele%mode3%v(1,1:6) = ele%mode3%v(1,1:6) + ele%mode3%x%eta

call mat_make_unit(d_mat)
d_mat(1:4,6) = [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap]
d_mat(5,1:4) = [-ele%x%etap, ele%x%eta, -ele%y%etap, ele%y%eta]
ele%mode3%v = matmul(d_mat, ele%mode3%v)

end subroutine twiss3_from_twiss2

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine twiss3_at_start (lat, error, ix_branch, tune3)
!
! Subroutine to calculate the 3D twiss parameters of the three modes of the full 6D 1-turn transfer matrix.
! This routine is for lattices with closed geometries. For open lattices see: twiss3_from_twiss2.
!
! Note: The rf must be on for this calculation.
!
! Input:
!   lat         -- lat_struct: Lattice with
!   ix_branch   -- integer, optional: Branch index. 0 = default.
!
! Output:
!   lat       -- lat-struct:
!     %branch(ix_branch)%ele(0)  -- Ele_struct: Starting element
!       %mode3    -- Mode3_struct: Structure holding the normal modes.
!         %v(6,6)    -- Real(rp): V coupling matrix.
!         %a            -- Twiss_struct: "a" normal mode Twiss parameters.
!         %b            -- Twiss_struct: "b" normal mode Twiss parameters.
!         %c            -- Twiss_struct: "c" normal mode Twiss parameters.
!   error     -- Logical: Set True if there is no RF. False otherwise.
!   tune3(3)  -- real(rp), optional: Normal mode tunes
!-

subroutine twiss3_at_start (lat, err_flag, ix_branch, tune3)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

real(rp), optional :: tune3(3)
real(rp) g(6,6), t3(3)

integer, optional :: ix_branch
integer n

logical err_flag
character(20) :: r_name = 'twiss3_at_start'

!

branch => lat%branch(integer_option(0, ix_branch))
ele => branch%ele(0)
err_flag = .true.

if (.not. associated(ele%mode3)) allocate(ele%mode3)

call transfer_matrix_calc (lat, branch%param%t1_with_RF, ix_branch = ix_branch, one_turn=.true.)
if (all(branch%param%t1_with_RF(6, 1:5) == 0)) then
  call out_io (s_error$, r_name, 'RF IS OFF FOR THE MODE3 CALCULATION!')
  return
endif
call normal_mode3_calc (branch%param%t1_with_RF, t3, g, ele%mode3%v)

if (present(tune3)) tune3 = t3

ele%mode3%x%eta = ele%mode3%v(1,6)
ele%mode3%y%eta = ele%mode3%v(3,6)

ele%mode3%x%etap = ele%mode3%v(2,6)
ele%mode3%y%etap = ele%mode3%v(4,6)

call mode1_calc (g(1:2, 1:2), ele%mode3%a)
call mode1_calc (g(3:4, 3:4), ele%mode3%b)
call mode1_calc (g(5:6, 5:6), ele%mode3%c)

err_flag = .false.

!-------------------------------------------------------------------------------------
contains

subroutine mode1_calc (gg, twiss)

type(twiss_struct) twiss
real(rp) gg(:, :)

!

twiss%beta = gg(2, 2)**2
twiss%alpha = gg(2, 1) * gg(2, 2)
twiss%gamma = (1 + twiss%alpha**2) / twiss%beta
twiss%phi = 0

end subroutine mode1_calc

end subroutine twiss3_at_start 

end module mode3_mod





