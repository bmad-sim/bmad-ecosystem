MODULE mode3_mod
USE bmad

REAL(rp), PARAMETER :: m = 0.707106781d0  ! 1/sqrt(2)
REAL(rp), PARAMETER :: o = 0.0d0  ! for compact code
REAL(rp), PARAMETER :: l = 1.0d0  ! for compact code

REAL(rp), PARAMETER :: Qr(6,6) = RESHAPE( [m,m,o,o,o,o, o,o,o,o,o,o, o,o,m,m,o,o,   &
                                           o,o,o,o,o,o, o,o,o,o,m,m, o,o,o,o,o,o],[6,6] )
REAL(rp), PARAMETER :: Qi(6,6) = RESHAPE( [o,o,o,o,o,o, m,-m,o,o,o,o, o,o,o,o,o,o,  &
                                           o,o,m,-m,o,o, o,o,o,o,o,o, o,o,o,o,m,-m],[6,6] )
REAL(rp), PARAMETER :: Qinv_r(6,6) = RESHAPE( [m,o,o,o,o,o, m,o,o,o,o,o, o,o,m,o,o,o,  &
                                               o,o,m,o,o,o, o,o,o,o,m,o, o,o,o,o,m,o],[6,6] )
REAL(rp), PARAMETER :: Qinv_i(6,6) = RESHAPE( [o,-m,o,o,o,o, o,m,o,o,o,o, o,o,o,-m,o,o,  &
                                               o,o,o,m,o,o, o,o,o,o,o,-m, o,o,o,o,o,m],[6,6] )
REAL(rp), PARAMETER :: S(6,6) = RESHAPE( [o,-l,o,o,o,o, l,o,o,o,o,o,  &
                                          o,o,o,-l,o,o, o,o,l,o,o,o,  &
                                          o,o,o,o,o,-l, o,o,o,o,l,o],[6,6] )

PRIVATE m, o, l
PRIVATE Qr, Qi
PRIVATE Qinv_r, Qinv_i
PRIVATE S

CONTAINS

!+
! Subroutine normal_mode3_calc (mat, tune, G, V, synchrotron_motion)
!
! Does an Eigen decomposition of the 1-turn transfer matrix (mat) and generates
! G and V.  Betatron and synchrotron tunes are places in tune.
!
! Input:
!  mat(6,6)            -- real(rp): 1-turn transfer matrix
!  synchrotron_motion  -- real(rp), optional: Default is to treat tune(3) as if it were a synchrotron tune:
!                                             tune(3) = 2pi - tune(3).  IF present and true, then this
!                                             correction is overridden.
! Output:
!  tune(3)             -- real(rp): Tunes of the 3 normal modes (radians)
!  G(6,6)              -- real(rp): G is block diagonal and related to the normal mode Twiss parameters.
!  V(6,6)              -- real(rp): Converts from normal mode coordinates to canonical coordinates: x = V.a
!
! Note: Vbar can be obtained from a direct call to make_V.
!-
SUBROUTINE normal_mode3_calc (mat, tune, G, V, synchrotron_motion)
  USE bmad

  IMPLICIT NONE

  REAL(rp) mat(6,6)
  REAL(rp) tune(3)
  REAL(rp) G(6,6)
  REAL(rp) V(6,6)
  LOGICAL, OPTIONAL :: synchrotron_motion

  INTEGER i

  REAL(rp) N(6,6)
  REAL(rp) Ninv(6,6)
  REAL(rp) gamma(3)

  LOGICAL error

  CALL make_N(mat, N, Ninv, gamma, error, tune, synchrotron_motion)
  CALL make_G(N,gamma,G)
  CALL make_V(N,gamma,V)

END SUBROUTINE normal_mode3_calc

!+
! Subroutine xyz_to_action(ring,ix,X,J,error)
!
! Given the canonical phase space coordinates X of a particle, this returns
! a vector from which Ja, Jb, Jc can be easily extracted.
!
! The J vector looks like:
! J = (sqrt(2Ja)cos(phia), -sqrt(2Ja)sin(phia), sqrt(2Jb)cos(phib), -sqrt(2Jb)sin(phib), sqrt(2Jc)cos(phic), -sqrt(2Jc)sin(phic))
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
!  ix       -- integer: element index at which to calculate J
!  X(1:6)   -- real(rp): canonical phase space coordinates of the particle
!
! Output:
!  J(1:6)   -- real(rp): Vector containing normal mode invariants and phases
!  error    -- logical: Set to true on error.  Often means Eigen decomposition failed.
!
!-
SUBROUTINE xyz_to_action(ring,ix,X,J,error)
  USE bmad

  IMPLICIT none

  TYPE(lat_struct) ring
  INTEGER ix
  REAL(rp) J(1:6)
  REAL(rp) X(1:6)
  REAL(rp) t6(1:6,1:6)
  REAL(rp) N(1:6,1:6)
  REAL(rp) Ninv(1:6,1:6)
  REAL(rp) eval_r(1:6), eval_i(1:6)
  REAL(rp) evec_r(1:6,1:6), evec_i(1:6,1:6)
  REAL(rp) nrml_evec_r(1:6,1:6), nrml_evec_i(1:6,1:6)
  REAL(rp) evec_inv_r(1:6,1:6), evec_inv_i(1:6,1:6)
  REAL(rp) gamma(3)

  INTEGER i

  LOGICAL error

  CALL transfer_matrix_calc (ring, .true., t6, ix1=ix)
  CALL make_N(t6, N, Ninv, gamma, error)
  IF( error ) THEN
    RETURN
  ENDIF

  J = MATMUL(Ninv,X)
END SUBROUTINE xyz_to_action

!+
! Subroutine action_to_xyz(ring,ix,J,X,error)
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
!  error    -- logical: Set to true on error.  Often means Eigen decomposition failed.
!
!-
SUBROUTINE action_to_xyz(ring,ix,J,X,error)
  USE bmad
  USE eigen_mod

  IMPLICIT none

  TYPE(lat_struct) ring
  INTEGER ix
  REAL(rp) J(1:6)
  REAL(rp) X(1:6)
  REAL(rp) t6(1:6,1:6)
  REAL(rp) N(1:6,1:6)
  REAL(rp) Ninv(1:6,1:6)
  REAL(rp) gamma(3)

  INTEGER i

  LOGICAL error

  CALL transfer_matrix_calc (ring, .true., t6, ix1=ix)
  CALL make_N(t6, N, Ninv, gamma, error)
  IF( error ) THEN
    RETURN
  ENDIF

  X = MATMUL(N,J)
END SUBROUTINE action_to_xyz

!+
! Subroutine make_N(t6,N,Ninv,error,tunes,synchrotron_motion)
!
! Given a 1-turn transfer matrix, this returns N and its inverse Ninv.
! N converts between normal invarients and phases and canonical coordinates:
! X = N.J
!
! N is obtained from the Eigen decomposition of the 1-turn transfer matrix.
! It is obtained by applying certain normalizations to the matrix of Eigen vectors, then making
! the result real using Q.
!
! Input:
!  t6(6,6)             -- real(rp): 1-turn transfer matrix
!  synchrotron_motion  -- real(rp), optional: Default is to treat tune(3) as if it were a synchrotron tune:
!                                             tune(3) = 2pi - tune(3).  IF present and true, then this
!                                             correction is overridden.
! Output:
!  N(6,6)              -- real(rp): X = N.J
!  Ninv(6,6)           -- real(rp): J = Ninv.X
!  error               -- logical: Set to true on error.  Often means Eigen decomposition failed.
!  gamma(3)            -- real(rp): gamma1, gamma2, and gamma3, related to calculation of G
!  tunes(3)            -- real(rp): Tune of the 3 normal modes of the beam (radians)
!-
SUBROUTINE make_N(t6,N,Ninv,gamma,error,tunes,synchrotron_motion)
  USE bmad
  USE eigen_mod

  IMPLICIT NONE

  REAL(rp) t6(6,6)
  REAL(rp) eval_r(6), eval_i(6)
  REAL(rp) evec_r(6,6), evec_i(6,6)
  REAL(rp) adj_evec_r(6,6), adj_evec_i(6,6)
  REAL(rp) nrml_evec_r(6,6), nrml_evec_i(6,6)
  REAL(rp) evec_inv_r(6,6), evec_inv_i(6,6)
  REAL(rp) N(6,6), Ninv(6,6)
  REAL(rp) gamma(3)
  REAL(rp), OPTIONAL :: tunes(3)
  LOGICAL, OPTIONAL :: synchrotron_motion

  LOGICAL error

  INTEGER i

  IF( error ) THEN
    RETURN
  ENDIF

  CALL mat_eigen (t6, eval_r, eval_i, evec_r, evec_i, error)
  evec_r = TRANSPOSE(evec_r)
  evec_i = TRANSPOSE(evec_i)

  IF( PRESENT(tunes) ) THEN
    !Calculate fractional tunes from eigenvalues.  eigenvectors are needed to determine if tune is
    !above or below 0.5
    CALL tunes_from_evals(eval_r, eval_i, evec_r, evec_i, tunes)
  ENDIF
  IF( PRESENT(tunes) .and. PRESENT(synchrotron_motion) ) THEN
    IF (logic_option(.false., synchrotron_motion)) tunes(3) = twopi-tunes(3)
  ENDIF

  !Adjust evec complex phase and normalize to make symplectic
  CALL normalize_and_phase_emat(evec_r, evec_i, adj_evec_r, adj_evec_i)
  CALL symp_emat_columns(adj_evec_r, adj_evec_i, nrml_evec_r, nrml_evec_i)

  !Get inverse so we can easily obtain Ninv
  !Because matrix is symplectic, its cplx symp conj is its inverse.
  CALL cplx_symp_conj(nrml_evec_r, nrml_evec_i, evec_inv_r, evec_inv_i)

  !Transform to real basis
  N = MATMUL(nrml_evec_r,Qr) - MATMUL(nrml_evec_i,Qi)
  Ninv = MATMUL(Qinv_r,evec_inv_r) - MATMUL(Qinv_i,evec_inv_i)

  gamma(1) = SQRT(ABS(N(1,1)*N(2,2) - N(1,2)*N(2,1)))
  gamma(2) = SQRT(ABS(N(3,3)*N(4,4) - N(3,4)*N(4,3)))
  gamma(3) = SQRT(ABS(N(5,5)*N(6,6) - N(5,6)*N(6,5)))

  !Ncheck = MATMUL(nrml_evec_r,Qi) + MATMUL(nrml_evec_i,Qr)  !this should be zero
END SUBROUTINE make_N

!+
! Subroutine tunes_from_evals(eval_r,eval_i,evec_r,evec_i,tunes)
!
! Calculates tunes from complex eigenvalues.  Uses determinant of 2x2 blocks down diagonal
! to determine if tune is above or below 0.5.
! Tunes are in units of radians per turn.
!
! Input:
!   eval_r(6)       -- REAL(rp) real part of eigenvalues
!   eval_i(6)       -- REAL(rp) imaginary part of eigenvalues
!   evec_r(6,6)     -- REAL(rp) real part of eigenvector matrix
!   evec_i(6,6)     -- REAL(rp) imaginary part of eigenvector matrix
! Output:
!   tunes(3)        -- REAL(rp) normal mode fractional tunes.  fractional_tune(i) = tunes(i)/(2 pi)
!-
SUBROUTINE tunes_from_evals(eval_r,eval_i,evec_r,evec_i,tunes)
  USE bmad

  IMPLICIT NONE

  REAL(rp) eval_r(6)
  REAL(rp) eval_i(6)
  REAL(rp) evec_r(6,6)
  REAL(rp) evec_i(6,6)
  REAL(rp) tunes(3)

  REAL(rp) det1, det2, det3


  det1 = AIMAG( CMPLX(evec_r(1,1),evec_i(1,1))*CMPLX(evec_r(2,2),evec_i(2,2)) - &
                CMPLX(evec_r(1,2),evec_i(1,2))*CMPLX(evec_r(2,1),evec_i(2,1)) )
  det2 = AIMAG( CMPLX(evec_r(3,3),evec_i(3,3))*CMPLX(evec_r(4,4),evec_i(4,4)) - &
                CMPLX(evec_r(3,4),evec_i(3,4))*CMPLX(evec_r(4,3),evec_i(4,3)) )
  det3 = AIMAG( CMPLX(evec_r(5,5),evec_i(5,5))*CMPLX(evec_r(6,6),evec_i(6,6)) - &
                CMPLX(evec_r(5,6),evec_i(5,6))*CMPLX(evec_r(6,5),evec_i(6,5)) )

  IF( det1 .GT. 0 ) THEN
    tunes(1) = twopi - atan2(eval_i(1),eval_r(1))
  ELSE
    tunes(1) = atan2(eval_i(1),eval_r(1))
  ENDIF
  IF( det2 .GT. 0 ) THEN
    tunes(2) = twopi - atan2(eval_i(3),eval_r(3))
  ELSE
    tunes(2) = atan2(eval_i(3),eval_r(3))
  ENDIF
  IF( det3 .GT. 0 ) THEN
    tunes(3) = twopi - atan2(eval_i(5),eval_r(5))
  ELSE
    tunes(3) = atan2(eval_i(5),eval_r(5))
  ENDIF
END SUBROUTINE tunes_from_evals

!+
! Subroutine make_G(N,gamma,G,Ginv)
!
! Construct G from N and and the normal mode gammas.
!
! G is block diagonal.  It is populated with the daggers of the N matrix, each normalized by the gamma numbers.
! See code for details.
!
! G contains the normal-mode Twiss functions.  Each 2x2 block looks like:
! 
! /  1/SQRT(beta)           0      \
! \  alpha/SQRT(beta)   SQRT(beta) /
!
! Input:
!  N(6,6)     -- real(rp): N matrix, which is a normalizes, phased, and symplectic matrix of the eigenvectors of the 1-turn matrix.
!  gamma(3)   -- real(rp): gamma numbers related to construction of N.
! Output:
!  G(6,6)     -- real(rp): Block diagonal 
!  Ginv(6,6)  -- real(rp): Inverse of G
!
!-
SUBROUTINE make_G(N,gamma,G,Ginv)
  USE bmad

  IMPLICIT NONE

  REAL(rp) N(6,6)
  REAL(rp) gamma(3)
  REAL(rp) G(6,6)
  REAL(rp), OPTIONAL :: Ginv(6,6)

  G = 0.0_rp

  G(1:2,1:2) = dagger2(N(1:2,1:2))/gamma(1)
  G(3:4,3:4) = dagger2(N(3:4,3:4))/gamma(2)
  G(5:6,5:6) = dagger2(N(5:6,5:6))/gamma(3)

  IF( PRESENT(Ginv) ) THEN
    Ginv = 0.0_rp
    Ginv(1:2,1:2) = dagger2(G(1:2,1:2))
    Ginv(3:4,3:4) = dagger2(G(3:4,3:4))
    Ginv(5:6,5:6) = dagger2(G(5:6,5:6))
  ENDIF
END SUBROUTINE make_G

!+
! Subroutine make_V(N,gamma,V,Vbar)
!
! Makes V matrix from N via G.
!
! V converts from normal mode coordinates to canonical coordinates: x = V.a
!
! Vbar is useful for analyzing the coupling properties of a storage ring.
!
! Input:
!  N(6,6)    -- real(rp): A form of the eigen matrix of the 1-turn transfer matrix
!  gamma(3)  -- real(rp): gamma numbers related to construction of N.
! Output:
!  V(6,6)    -- real(rp): Matrix that converts from normal mode coordintes to canonical coordinates: x=Va
!  Vbar(6,6) -- real(rp): Useful in analyzing the coupling properties of a storage ring.
!-
SUBROUTINE make_V(N,gamma,V,Vbar)
  USE bmad

  IMPLICIT NONE

  REAL(rp) N(6,6)
  REAL(rp) gamma(3)
  REAL(rp) G(6,6)
  REAL(rp) Ginv(6,6)
  REAL(rp) V(6,6)
  REAL(rp), OPTIONAL :: Vbar(6,6)

  CALL make_G(N,gamma,G,Ginv)
  V = MATMUL(N,G) 

  IF( PRESENT(Vbar) ) THEN
    Vbar = MATMUL(G,MATMUL(V,MATMUL(TRANSPOSE(S),MATMUL(TRANSPOSE(G),S))))
  ENDIF
END SUBROUTINE make_V

!+
! Subroutine normal_sigma_mat(sigma_mat,normal)
!
! Given a beam envelop sigma matrix sigma_mat, this returns the 3 normal mode
! emittances.
!
! The normal mode emittance of the sigma matrix are the eigenvalues of
! sigma_mat . S
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
! Output:
!  normal(3)        -- real(rp): normal mode emittances
!-
SUBROUTINE normal_sigma_mat(sigma_mat,normal)
  USE bmad
  USE eigen_mod

  IMPLICIT none

  REAL(rp) sigma_mat(1:6,1:6)
  REAL(rp) normal(1:3)
  REAL(rp) eval_r(1:6), eval_i(1:6)
  REAL(rp) evec_r(1:6,1:6), evec_i(1:6,1:6)
  REAL(rp) sigmaS(1:6,1:6)

  LOGICAL error

  sigmaS = MATMUL(sigma_mat,S)

  CALL mat_eigen(sigmaS,eval_r,eval_i,evec_r,evec_i,error)

  normal(1) = ABS(eval_i(1))
  normal(2) = ABS(eval_i(3))
  normal(3) = ABS(eval_i(5))
END SUBROUTINE normal_sigma_mat

!+
! Subrouting get_abc_from_updated_smat(ring, ix, sigma_mat, normal, error)
!
! This subroutine is experimental.  It obtains the normal mode emittances from the sigma matrix using the 
! eigenvectors of the 1-turn transfer matrix.  It is used to obtain the emittances from a sigma matrix
! that has been perturbed.  Under ordinary circumstances, the eigenvectors of the sigma matrix and 1-turn matrix
! should be the same.  However, if the sigma matrix has been perturbed, then they may not be the same.  The question
! is then:  Are the normal mode emittances the eigenvalues of the sigma matrix? Or should the eigenvectors
! of the 1-turn matrix be applied to the sigmatrix, and the values of the resulting almost-diagonal matrix
! taken as the emittance?  Or maybe the emittance is not well-defined if the sigma-matrix is perturbed.
!
! Input:
!  ring           -- lat_struct: the ring
!  ix             -- integer: element at which to do the transformation
!  sigma_mat(6,6) -- real(rp): beam envelop sigma matrix (possibly perturbed)
! Output:
!  normal(3)      -- real(rp): something like the normal mode emittances of the sigma matrix
!  error          -- logical:  set to true if something goes wrong.  Usually means Eigen decomposition of the 1-turn matrix failed.
!-
SUBROUTINE get_abc_from_updated_smat(ring, ix, sigma_mat, normal, error)

  USE bmad
  USE eigen_mod

  TYPE(lat_struct) ring
  INTEGER ix
  REAL(rp) normal(3)
  REAL(rp) sigma_mat(6,6)

  REAL(rp) eval_r(6)
  REAL(rp) eval_i(6)
  REAL(rp) evec_r(6,6)
  REAL(rp) evec_i(6,6)
  REAL(rp) evec_inv_r(6,6)
  REAL(rp) evec_inv_i(6,6)
  REAL(rp) t6(6,6)
  REAL(rp) smatS(6,6)

  REAL(rp) term3(6,6)
  REAL(rp) term4(6,6)
  REAL(rp) termSum(6,6)

  LOGICAL ok, error
  INTEGER i

  error = .false.

  CALL transfer_matrix_calc (ring, .true., t6, ix1=ix)

  CALL mat_eigen (t6, eval_r, eval_i, evec_r, evec_i, error)
  IF( error ) THEN
    WRITE(*,'(A,I6,A)') "BAD: Eigenvectors of transfer matrix not found for element ", ix, ring%ele(ix)%name
    RETURN
  ENDIF

  evec_r = TRANSPOSE(evec_r)
  evec_i = TRANSPOSE(evec_i)
  CALL symp_emat_columns(evec_r, evec_i, evec_r, evec_i)
  CALL cplx_symp_conj(evec_r, evec_i, evec_inv_r, evec_inv_i)

  smatS = MATMUL(sigma_mat,S)

  term3 = MATMUL(evec_inv_r,MATMUL(smatS,evec_i))
  term4 = MATMUL(evec_inv_i,MATMUL(smatS,evec_r))
  termSum = term3+term4

  normal(1) = termSum(2,2)
  normal(2) = termSum(4,4)
  normal(3) = termSum(6,6)

END SUBROUTINE get_abc_from_updated_smat

!+
! Subroutine make_smat_from_abc(t6, mode, sigma_mat, error)
!
! Given the 1-turn transfer matrix and a normal_modes_struct containing the normal mode
! emittances, this routine returns the beam envelop sigma matrix.
!
! Input:
!  t6(6,6)          -- real(rp): 1-turn transfer matrix
!  mode             -- normal_modes_struct: normal mode emittances
!      %a%emittance -- real(rp): a-mode emittance
!      %b%emittance -- real(rp): b-mode emittance
!      %z%emittance -- real(rp): z-mode emittance
! Output:
!  sigma_mat(6,6)   -- real(rp): beam envelop sigma matrix
!  error            -- logical:  set to true if something goes wrong.  Usually means Eigen decomposition of the 1-turn matrix failed.
!-
SUBROUTINE make_smat_from_abc(t6, mode, sigma_mat, error)

  USE bmad
  USE eigen_mod

  TYPE(normal_modes_struct) mode
  REAL(rp) sigma_mat(6,6)

  REAL(rp) eval_r(6)
  REAL(rp) eval_i(6)
  REAL(rp) evec_r(6,6)
  REAL(rp) evec_i(6,6)
  REAL(rp) evec_inv_r(6,6)
  REAL(rp) evec_inv_i(6,6)
  REAL(rp) t6(6,6)
  REAL(rp) Drl(6,6)

  LOGICAL ok, error
  INTEGER i

  error = .false.

  CALL mat_eigen (t6, eval_r, eval_i, evec_r, evec_i, error)
  IF( error ) THEN
    WRITE(*,'(A,I6,A)') "BAD: Eigenvectors of transfer matrix not found for element."
    RETURN
  ENDIF

  evec_r = TRANSPOSE(evec_r)
  evec_i = TRANSPOSE(evec_i)
  CALL symp_emat_columns(evec_r, evec_i, evec_r, evec_i)
  CALL cplx_symp_conj(evec_r, evec_i, evec_inv_r, evec_inv_i)

  Drl = 0.0d0
  Drl(1,1) = -mode%a%emittance
  Drl(2,2) =  mode%a%emittance
  Drl(3,3) = -mode%b%emittance
  Drl(4,4) =  mode%b%emittance
  Drl(5,5) = -mode%z%emittance
  Drl(6,6) =  mode%z%emittance
  sigma_mat = MATMUL((MATMUL(MATMUL(evec_i,Drl),evec_inv_r) + MATMUL(MATMUL(evec_r,Drl),evec_inv_i)),S)

END SUBROUTINE make_smat_from_abc

!+
!
! Subroutine cplx_symp_conj(evec_r, evec_i, symp_evec_r, symp_evec_i)
!
! Return the complex symplectic conjugate of a 6x6 matrix.
! This is useful because if the input matrix is symplectic with respect to S, then the complex symplectic conjugate
! is the matrix inverse.
!
! This is much more robust and accurate than calculating the inverse by back substitution.
!
! Input:
!  evec_r(6,6)      -- real(rp): real part of input matrix
!  evec_i(6,6)      -- real(rp): imaginary part of input matrix
! Output:
!  symp_evec_r(6,6) -- real(rp): real part of the symplectic complex conjugate of the input matrix
!  symp_evec_i(6,6) -- real(rp): imaginary part of the symplectic complex conjugate of the input matrix
!
!-
SUBROUTINE cplx_symp_conj(evec_r, evec_i, symp_evec_r, symp_evec_i)

  USE bmad

  IMPLICIT NONE

  REAL(rp) evec_r(6,6)
  REAL(rp) evec_i(6,6)
  REAL(rp) symp_evec_r(6,6)
  REAL(rp) symp_evec_i(6,6)

  symp_evec_r(1:2,1:2) = dagger2(evec_i(1:2,1:2))
  symp_evec_r(1:2,3:4) = dagger2(evec_i(3:4,1:2))
  symp_evec_r(1:2,5:6) = dagger2(evec_i(5:6,1:2))
  symp_evec_r(3:4,1:2) = dagger2(evec_i(1:2,3:4))
  symp_evec_r(3:4,3:4) = dagger2(evec_i(3:4,3:4))
  symp_evec_r(3:4,5:6) = dagger2(evec_i(5:6,3:4))
  symp_evec_r(5:6,1:2) = dagger2(evec_i(1:2,5:6))
  symp_evec_r(5:6,3:4) = dagger2(evec_i(3:4,5:6))
  symp_evec_r(5:6,5:6) = dagger2(evec_i(5:6,5:6))

  symp_evec_i(1:2,1:2) = -dagger2(evec_r(1:2,1:2))
  symp_evec_i(1:2,3:4) = -dagger2(evec_r(3:4,1:2))
  symp_evec_i(1:2,5:6) = -dagger2(evec_r(5:6,1:2))
  symp_evec_i(3:4,1:2) = -dagger2(evec_r(1:2,3:4))
  symp_evec_i(3:4,3:4) = -dagger2(evec_r(3:4,3:4))
  symp_evec_i(3:4,5:6) = -dagger2(evec_r(5:6,3:4))
  symp_evec_i(5:6,1:2) = -dagger2(evec_r(1:2,5:6))
  symp_evec_i(5:6,3:4) = -dagger2(evec_r(3:4,5:6))
  symp_evec_i(5:6,5:6) = -dagger2(evec_r(5:6,5:6))
END SUBROUTINE cplx_symp_conj

!+
! Function dagger2(A) RESULT(Ad)
!
! Return the complex symplectic conjugate of a 2x2 matrix.
!
! A_dagger = / A22  -A12 \
!            \ -A21  A11 /
!
! Input:
!  A(2,2)   -- real(rp): 2x2 matrix
! Output:
!  Ad(2,2)  -- real(rp): A_dagger
!-
FUNCTION dagger2(A) RESULT(Ad)
  USE bmad

  REAL(rp) A(2,2)
  REAL(rp) Ad(2,2)

  Ad(1,1) =  A(2,2)
  Ad(1,2) = -A(1,2)
  Ad(2,1) = -A(2,1)
  Ad(2,2) =  A(1,1)
END FUNCTION dagger2

!+
! Subroutine normalize_and_phase_emat(evec_r, evec_i, adj_evec_r, adj_evec_i)
!
! This subroutine assumes that the eigenvectors are arranged in complex-conjugate pairs: (e1 e1* e2 e2* e3 e3*)
!
! The following normalizations are applied to the matrix of eigenvectors.  The result is still a matrix of eigenvectors,
! but it is normalized and the phase adjusted such that G and V can be easily extracted.
! It is probably also unique, though this remains to be proven.
!
! 1) Normalize the eigenvectors such that their length is 1.
! 2) Adjust the phase of the eigenvectors such that the (1,1 and 1,2) and (3,3 and 3,4) and (5,5 and 5,6) elements are real.
! 3) Fixes the sign of the pairs so that the 1,1 3,3 and 5,5 elements are positive.
!
! Input:
!  evec_r(6,6)      -- real(rp): real part of Eigen matrix
!  evec_i(6,6)      -- real(rp): imaginary part of Eigen matrix
!  tunes(3)         -- real(rp): Tunes of the 3 normal modes.
! Output:
!  adj_evec_r(6,6)  -- real(rp): real part of Eigen matrix: normalized, phased, and signs fixed
!  adj_evec_i(6,6)  -- real(rp): imaginary part of Eigen matrix: normalized, phased, and signs fixed
!  tunes(3)         -- real(rp): Tunes of the 3 normal modes.  The ambiguity in tune is resolved so that det(G) = 1.
!
!-
SUBROUTINE normalize_and_phase_emat(evec_r, evec_i, adj_evec_r, adj_evec_i)
  REAL(rp) evec_r(6,6)
  REAL(rp) evec_i(6,6)
  REAL(rp) adj_evec_r(6,6)
  REAL(rp) adj_evec_i(6,6)

  REAL(rp) theta  
  REAL(rp) costh, sinth
  REAL(rp) nrml
  REAL(rp) fix

  INTEGER i, j, ix

  DO i=1,3
    ix = i*2-1

    ! Calculate norm of the eigenvector.  Result will be used to normalize length to 1.
    nrml = SQRT(SUM(evec_r(:,ix)*evec_r(:,ix) + evec_i(:,ix)*evec_i(:,ix)))

    ! For each element of the eigenvector, rotate the eigenvector in the complex plane
    ! by an angle that makes the (1,1 and 1,2) or (3,3 and 3,4) or (5,5 and 5,6) elements of the eigen matrix real.
    theta = -ATAN2(evec_i(ix,ix),evec_r(ix,ix))

    ! Apply the normalization and rotation
    costh = COS(theta+pi)
    sinth = SIN(theta+pi)
    adj_evec_r(:,ix) = (evec_r(:,ix)*costh - evec_i(:,ix)*sinth) / nrml
    adj_evec_i(:,ix) = (evec_r(:,ix)*sinth + evec_i(:,ix)*costh) / nrml
    costh = COS(-theta+pi)
    sinth = SIN(-theta+pi)
    adj_evec_r(:,ix+1) = (evec_r(:,ix+1)*costh - evec_i(:,ix+1)*sinth) / nrml
    adj_evec_i(:,ix+1) = (evec_r(:,ix+1)*sinth + evec_i(:,ix+1)*costh) / nrml
  ENDDO

  !Fix signs to make the 1,1 3,3 and 5,5 elements positive
  DO i=1,3
    ix = i*2-1
    fix = SIGN(1.0d0,adj_evec_r(ix,ix))
    adj_evec_r(:,ix) = fix*adj_evec_r(:,ix)
    adj_evec_r(:,ix+1) = fix*adj_evec_r(:,ix+1)
    adj_evec_i(:,ix) = fix*adj_evec_i(:,ix)
    adj_evec_i(:,ix+1) = fix*adj_evec_i(:,ix+1)
  ENDDO
END SUBROUTINE normalize_and_phase_emat

!+
! Subroutine symp_emat_columns(evec_r, evec_i, nrm_evec_r, nrm_evec_i)
!
! Applies a normalization to the columns of a 6x6 matrix of eigenvectors, where the columns are complex
! conjugate pairs.  The normalization is such that the sum of the determinants of the 2x2 blocks down each pair of rows is 1.  
! The resulting matrix is symplectic with respect to S,
!  tr(E).S.E = iS
!
! Input:
!  evec_r(6,6)      -- real(rp): real part of Eigen matrix
!  evec_i(6,6)      -- real(rp): imaginary part of Eigen matrix
! Output:
!  nrm_evec_r(6,6)  -- real(rp): real part of Eigen matrix made symplectic w.r.t. S
!  nrm_evec_i(6,6)  -- real(rp): imaginary part of Eigen matrix made symplectic w.r.t. S
!-
SUBROUTINE symp_emat_columns(evec_r, evec_i, nrm_evec_r, nrm_evec_i)

  USE bmad

  IMPLICIT NONE

  REAL(rp) evec_r(6,6)
  REAL(rp) evec_i(6,6)
  REAL(rp) part_r(1,6)
  REAL(rp) part_i(1,6)
  REAL(rp) nrm_evec_r(6,6)
  REAL(rp) nrm_evec_i(6,6)
  REAL(rp) norm
  INTEGER i

  DO i=1,5,2
    part_r = MATMUL(TRANSPOSE(evec_r(:,i:i)),S)  
    part_i = MATMUL(TRANSPOSE(evec_i(:,i:i)),S)  
    norm = ABS(DOT_PRODUCT(part_r(1,:),evec_i(:,i+1)) + DOT_PRODUCT(part_i(1,:),evec_r(:,i+1)))
    nrm_evec_r(:,i)   = evec_r(:,i)   / SQRT(norm)
    nrm_evec_i(:,i)   = evec_i(:,i)   / SQRT(norm)
    nrm_evec_r(:,i+1) = evec_r(:,i+1) / SQRT(norm)
    nrm_evec_i(:,i+1) = evec_i(:,i+1) / SQRT(norm)
  ENDDO

END SUBROUTINE symp_emat_columns

!+
! Subroutine project_via_EDES(ring, ix, mode, sigma_x, sigma_y, sigma_z)
!
! Obtains the projected x,y, and z beamsizes by building the sigma matrix
! from the normal mode emittances and 1-turn transfer matrix.
! These projectes beamsize are what would be seen by instrumentation.
!
! This method of projecting takes into account transverse and longitudinal coupling.
!
! This method of obtaining the projected beam sizes is from "Alternitive approach to general
! coupled linear optics" by Andrzej Wolski.
!
! The normal mode emittances used to generate a beam envelop sigma matrix from the 
! 1-turn transfer matrix.  The projected sizes are from the 1,1 3,3 and 5,5 elements of
! the sigma matrix.
!
! This gives the same result as the project_via_Vbar subroutine.
!
! Input:
!  ring             -- lat_struct: the storage ring
!  ix               -- integer: element at which to make the projection
!  mode             -- normal_modes_struct: normal mode emittances
!      %a%emittance -- real(rp): a-mode emittance
!      %b%emittance -- real(rp): b-mode emittance
!      %z%emittance -- real(rp): z-mode emittance
! Output:
!  sigma_x          -- real(rp): projected horizontal beamsize
!  sigma_y          -- real(rp): projected vertical beamsize
!  sigma_z          -- real(rp): projected longitudinal beamsize
!-
SUBROUTINE project_via_EDES(ring, ix, mode, sigma_x, sigma_y, sigma_z)

  USE bmad
  USE eigen_mod

  IMPLICIT NONE

  TYPE(lat_struct) ring
  INTEGER ix
  TYPE(normal_modes_struct) mode
  REAL(rp) sigma_x, sigma_y, sigma_z
  REAL(rp) t6(6,6)

  REAL(rp) sigma_mat(6,6)
  LOGICAL error

  CALL transfer_matrix_calc (ring, .true., t6, ix1=ix)
  CALL make_smat_from_abc(t6, mode, sigma_mat, error)

  sigma_x = SQRT(sigma_mat(1,1))
  sigma_y = SQRT(sigma_mat(3,3))
  sigma_z = SQRT(sigma_mat(5,5))
END SUBROUTINE project_via_EDES

!+
! Subroutine project_via_Vbar(ring, ix, mode, sigma_x, sigma_y, sigma_z)
!
! Obtains the projected x,y, and z beamsizes using Vbar and G6 and 
! the normal mode emittances.  This is a 6x6 extension of the 4x4 method described
! in "Linear analysis of coupled lattices" by D. Sagan and D. Rubin.
!
! This routine gives the same result as subroutine project_via_EDES.
!-
SUBROUTINE project_via_Vbar(ring, ix, mode, sigma_x, sigma_y, sigma_z)

  USE bmad

  TYPE(lat_struct) ring
  INTEGER ix
  TYPE(normal_modes_struct) mode
  REAL(rp) sigma_x, sigma_y, sigma_z

  REAL(rp) t6(1:6,1:6)
  REAL(rp) G6mat(1:6,1:6)
  REAL(rp) G6inv(1:6,1:6)
  REAL(rp) V6mat(1:6,1:6)
  REAL(rp) V6bar(1:6,1:6)
  REAL(rp) GiVb(1:6,1:6)
  
  REAL(rp) N(6,6)
  REAL(rp) Ninv(6,6)
  REAL(rp) gamma(3)

  REAL(rp) a_to_x, b_to_x, c_to_x
  REAL(rp) a_to_y, b_to_y, c_to_y
  REAL(rp) a_to_z, b_to_z, c_to_z
  REAL(rp) TermA, TermB, TermC

  LOGICAL ok, error

  !Calculate terms for horizontal projection at vBSM source point
  CALL transfer_matrix_calc (ring, .true., t6, ix1=ix)
  CALL make_N(t6, N, Ninv, gamma, error)
  IF(.not. ok) THEN
    WRITE(*,*) "BAD: make_N failed.  It is likely that the Eigen decomposition of the 1-turn matrix failed."
  ENDIF
  CALL make_G(N,gamma,G6mat,G6inv)
  CALL make_V(N,gamma,V6mat,V6bar)

  GiVb = MATMUL(G6inv,V6bar)

  !Terms for x-projection
  a_to_x = GiVb(1,1)**2+GiVb(1,2)**2
  b_to_x = GiVb(1,3)**2+GiVb(1,4)**2
  c_to_x = GiVb(1,5)**2+GiVb(1,6)**2

  !Terms for y-projection
  a_to_y = GiVb(3,1)**2+GiVb(3,2)**2
  b_to_y = GiVb(3,3)**2+GiVb(3,4)**2
  c_to_y = GiVb(3,5)**2+GiVb(3,6)**2

  !Terms for z-projection
  a_to_z = GiVb(5,1)**2+GiVb(5,2)**2
  b_to_z = GiVb(5,3)**2+GiVb(5,4)**2
  c_to_z = GiVb(5,5)**2+GiVb(5,6)**2

  !Calculate projected horizontal size
  TermA = mode%a%emittance * a_to_x
  TermB = mode%b%emittance * b_to_x
  TermC = mode%z%emittance * c_to_x
  sigma_x = SQRT(TermA + TermB + TermC)

  !Calculate projected vertical size
  TermA = mode%a%emittance * a_to_y
  TermB = mode%b%emittance * b_to_y
  TermC = mode%z%emittance * c_to_y
  sigma_y = SQRT(TermA + TermB + TermC)

  !Calculate projected bunch length
  TermA = mode%a%emittance * a_to_z
  TermB = mode%b%emittance * b_to_z
  TermC = mode%z%emittance * c_to_z
  sigma_z = SQRT(TermA + TermB + TermC)
END SUBROUTINE project_via_Vbar


!----------------------------------------------
! Subroutines below are from original mode3_mod
!----------------------------------------------

!+
! Subroutine twiss3_propagate_all (lat)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
!-

SUBROUTINE twiss3_propagate_all (lat)

IMPLICIT NONE

TYPE (lat_struct) lat

INTEGER i

DO i = 1, lat%n_ele_track
  call twiss3_propagate1 (lat%ele(i-1), lat%ele(i))
ENDDO

END SUBROUTINE twiss3_propagate_all

!+
! Subroutine twiss3_propagate1 (ele1, ele2)
!
! Subroutine to propagate the twiss parameters using all three normal modes.
!-

SUBROUTINE twiss3_propagate1 (ele1, ele2)

IMPLICIT NONE

TYPE mat2_struct
  REAL(rp) m(2,2)
END TYPE

TYPE (ele_struct) ele1, ele2
TYPE (mat2_struct) w(3)

REAL(rp) gamma(3), tv(6,6), w_inv(2,2)

INTEGER i, ik
LOGICAL err

!

IF (.NOT. ASSOCIATED(ele2%mode3)) ALLOCATE(ele2%mode3)

tv = MATMUL (ele2%mat6, ele1%mode3%v)

DO i = 1, 3
  ik = 2 * i - 1
  w(i)%m = tv(ik:ik+1,ik:ik+1)
  gamma(i) = SQRT(determinant (w(i)%m))
  w(i)%m = w(i)%m / gamma(i)
  call mat_symp_conj (w(i)%m, w_inv)
  ele2%mode3%v(1:6, ik:ik+1) = matmul(tv(1:6, ik:ik+1), w_inv)
ENDDO

ele2%mode3%x%eta = ele2%mode3%v(1,6)
ele2%mode3%y%eta = ele2%mode3%v(3,6)

ele2%mode3%x%etap = ele2%mode3%v(1,5)
ele2%mode3%y%etap = ele2%mode3%v(3,5)

call twiss1_propagate (ele1%mode3%a, w(1)%m,  ele2%value(l$), ele2%mode3%a, err)
call twiss1_propagate (ele1%mode3%b, w(2)%m,  ele2%value(l$), ele2%mode3%b, err)
call twiss1_propagate (ele1%mode3%c, w(3)%m,  0.0_rp,         ele2%mode3%c, err)

END SUBROUTINE

!+
! Subroutine twiss3_at_start (lat, error)
!
! Subroutine to calculate the twiss parameters of the three modes of the full 6D transfer
! matrix.
! Note: The rf must be on for this calculation.
!
! Modules needed:
!   use mode3_mod
!
! Input:
!   lat -- lat_struct: Lattice with
!
! Output:
!   lat   -- lat-struct:
!     %ele(0)  -- Ele_struct: Starting element
!       %mode3    -- Mode3_struct: Structure holding the normal modes.
!         %v(6,6)    -- Real(rp): V coupling matrix.
!         %a            -- Twiss_struct: "a" normal mode Twiss parameters.
!         %b            -- Twiss_struct: "b" normal mode Twiss parameters.
!         %c            -- Twiss_struct: "c" normal mode Twiss parameters.
!   error -- Logical: Set True if there is no RF. False otherwise.
!-

SUBROUTINE twiss3_at_start (lat, error)

IMPLICIT NONE

TYPE (lat_struct) lat
REAL(rp) g(6,6), tune3(3)
INTEGER n
LOGICAL error
CHARACTER(20) :: r_name = 'twiss3_at_start'

!

error = .true.

IF (.NOT. ASSOCIATED(lat%ele(0)%mode3)) ALLOCATE(lat%ele(0)%mode3)

CALL transfer_matrix_calc (lat, .true., lat%param%t1_with_RF)
if (ALL(lat%param%t1_with_RF(6,1:5) == 0)) then
  call out_io (s_error$, r_name, 'RF IS OFF FOR THE MODE3 CALCULATION!')
  RETURN
ENDIF
CALL normal_mode3_calc (lat%param%t1_with_RF, tune3, g, lat%ele(0)%mode3%v)

lat%ele(0)%mode3%x%eta = lat%ele(0)%mode3%v(1,6)
lat%ele(0)%mode3%y%eta = lat%ele(0)%mode3%v(3,6)

lat%ele(0)%mode3%x%etap = lat%ele(0)%mode3%v(1,5)
lat%ele(0)%mode3%y%etap = lat%ele(0)%mode3%v(3,5)

CALL mode1_calc (g(1:2, 1:2), tune3(1), lat%ele(0)%mode3%a)
CALL mode1_calc (g(3:4, 3:4), tune3(2), lat%ele(0)%mode3%b)
CALL mode1_calc (g(5:6, 5:6), tune3(3), lat%ele(0)%mode3%c)

error = .false.

!-------------------------------------------------------------------------------------
CONTAINS

  SUBROUTINE mode1_calc (gg, tune, twiss)

  TYPE (twiss_struct) twiss
  REAL(rp) gg(:,:), tune

  !

  twiss%beta = gg(2,2)**2
  twiss%alpha = gg(2,1) * gg(2,2)
  twiss%gamma = (1 + twiss%alpha**2) / twiss%beta
  twiss%phi = 0

  END SUBROUTINE

END SUBROUTINE

END MODULE mode3_mod
