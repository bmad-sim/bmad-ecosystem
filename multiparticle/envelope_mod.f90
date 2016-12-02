module envelope_mod

use twiss_and_track_mod
use, intrinsic :: iso_c_binding
use fgsl

implicit none

real(fgsl_double), parameter :: eps7 = 1.0d-7
integer(fgsl_size_t), parameter :: limit = 1000_fgsl_size_t

real(rp), parameter :: o = 0.0d0  ! for compact code
real(rp), parameter :: l = 1.0d0  ! for compact code
real(rp), parameter :: S6(6,6) = reshape( [o, -l, o, o, o, o, l, o, o, o, o, o, &
                                          o, o, o, -l, o, o, o, o, l, o, o, o, &
                                          o, o, o, o, o, -l, o, o, o, o, l, o], [6,6] )
real(rp), parameter :: S2(2,2)= reshape( [o, -l, l, o], [2,2] )
real(rp), parameter :: I2(2,2)= reshape( [l,  o, o, l], [2,2] )
real(rp), parameter :: I3(3,3) = reshape( [l, o, o, o, l, o, o, o, l], [3,3] )
real(rp), parameter :: I6(6,6) = reshape( [l, o, o, o, o, o, o, l, o, o, o, o, &
                                           o, o, l, o, o, o, o, o, o, l, o, o, &
                                           o, o, o, o, l, o, o, o, o, o, o, l], [6,6] )

contains

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine make_SR_mats(ele,co,etay,M,Bbar)
!
! Equation references are to "From the beam-envelope matrix to synchrotron-radiation
! integrals" by K. Ohmi, K. Hirata, and K. Oide.  Phys.Rev.E v49n1, January 1994.
!
! This subroutine makes the synchrotron radiation damping and diffusion matrices.
! The damping is contained in M, and the diffusion is contained in Bbar.
!
! See transport_with_sr_ele in this module for how to track through an element 
! using these matrices.
!
! IBS simulation requires a realistic vertical emittance, but ideal lattices
! do not simulate accurately the vertical emittance growth from synchrotron radiation.
! If etay is greater than zero, then the transfer matrix is modified: Mt = M.W
! W is a symplectic matrix that couples y and py to pz.
! Here is how to use etay:  With IBS turned off, adjust etay such that the equilibrium
! emittance is the desired value.  For example, for the SLS upgrade prototype dc12c,
! setting etay to 0.0033 results in a vertical emittance of 10 pm.
!
! Modules needed:
!   use envelope_mod
!
! Input:
!   ele                    -- ele_struct: 
!     %value(l$)           -- real(rp): element or slice length
!     %mat6(6,6)           -- real(rp): 6x6 transfer matrix
!     %value(b_field$)     -- real(rp): bend field
!     %value(b1_gradient$) -- real(rp): quadrupole field
!   co                     -- coord_struct:
!     %vec(6)              -- real(rp): p_z
!   etay                   -- real(rp): A vertical-dispersion like quantity that adds
!                                       vertical vertical excitation to the diffusion matrix Bbar.
! Output:
!   M(6,6)                 -- real(rp): transfer matrix with damping.
!   Bbar(6,6)              -- real(rp): diffusion matrix.
! 
!-
subroutine make_SR_mats(ele,co,etay,M,Bbar,Dbar)
  use bmad

  type(ele_struct) ele
  type(coord_struct) co
  real(rp) etay
  real(rp) M(:,:), Bbar(:,:)
  real(rp), optional :: Dbar(:,:)

  real(rp) Sigma_exit(6,6), Sigma_ent(6,6)
  real(rp) M0(6,6), W(6,6), Mt(6,6)
  real(rp) B0, B1, delta
  real(rp) l, gamma, rho

  integer i

  M0 = ele%mat6
  if(any(ele%key==(/sbend$,rbend$/))) then
    l = ele%value(l$)
    call convert_total_energy_to(ele%value(E_TOT$), -1, gamma)
    rho = ele%value(rho$)
    B0 = ele%value(b_field$)
    B1 = ele%value(b1_gradient$)
    delta = co%vec(6)

    M = M0 - l*matmul(M0,damping_matrix_D(gamma,rho,B0,B1,delta))  !Eqn. 80 & 81

    if(etay .gt. 0) then
      W = I6
      W(3,6) = -etay
      W(4,6) = -etay/5.0d0 !etay prime
      W(5,3) =  etay/5.0d0 !etay prime
      W(5,4) = -etay
      Mt = matmul(M,W)
      Bbar = l*matmul(Mt,matmul(diffusion_matrix_B(gamma,rho),transpose(Mt)))
    else
      Bbar = l*matmul(M,matmul(diffusion_matrix_B(gamma,rho),transpose(M)))  !Eqn. 88
    endif
    if(present(Dbar)) Dbar = l*matmul(M0,damping_matrix_D(gamma,rho,B0,B1,delta))
  else
    M = M0
    Bbar = 0.0d0
    if(present(Dbar)) Dbar = 0.0d0
  endif
end subroutine

!+
! Calculate diffusion matrix B from Ohmi Eqn 114 & 115.
!-
function diffusion_matrix_B(gamma,rho) result(mat)
  real(rp) mat(6,6)
  real(rp) gamma, rho
  real(rp), parameter :: D_const = 55.0d0/24.0d0/sqrt(3.0d0)*r_e*(1.055e-34)/(9.11e-31)/c_light

  mat = 0.0d0

  mat(6,6) = D_const * gamma**5 / abs(rho)**3
end function

!+
! Calculate damping matrix D from Ohmi Eqn 121.
!
! This routine assumes zero scalar potential.
!-
function damping_matrix_D(gamma, rho, B0, B1, delta) result(mat)
  real(rp) mat(6,6)
  real(rp) gamma, rho, B0, B1, delta
  real(rp) delta_f

  mat = 0.0d0

  delta_f = 1.0d0/(1.0d0+delta)

  mat(2,2) = delta_f
  mat(4,4) = delta_f
  mat(6,6) = 2.0d0*delta_f

  mat(6,1) = 1.0d0/rho + 2.0d0/B0*B1
  mat(6,3) = 2.0d0/B0*B1

  mat = energy_loss_rate(gamma,rho)*mat
end function

!+
! Calculate energy-loss rate from Ohmi Eqn 43.
!
! This implementation assumes dl/ds = 1
!-
function energy_loss_rate(gamma,rho) result(P)
  real(rp) P
  real(rp) gamma, rho
  real(rp), parameter :: P_const = 2.0d0/3.0d0*r_e

  P = P_const * gamma**3 / rho**2
end function

!+
! subroutine transport_with_sr(ele,M,Bbar,Sigma_ent,Sigma_exit)
!
! Transport a 6x6 beam envelope matrix through an element including
! damping and diffusion.
!
! Only bends include damping & diffusion.
!
! Input:
!   ele            -- integer: element key
!   M(6,6)         -- real(rp): transfer matrix with damping, populated by make_SR_mats
!   Bbar(6,6)      -- real(rp): diffusion matrix, populated by make_SR_mats
!   Sigma_ent(6,6) -- real(rp): incoming beam envelope matrix
! Output:
!   Sigma_exit(6,6) -- real(rp): outgoing beam envelope matrix
!-
subroutine transport_with_sr(ele,M,Bbar,Sigma_ent,Sigma_exit)
  use bmad

  type(ele_struct) ele
  real(rp) Sigma_exit(6,6), Sigma_ent(6,6)
  real(rp) M(:,:), Bbar(:,:)

  if(any(ele%key==(/sbend$,rbend$/))) then
    Sigma_exit = matmul(M,matmul(Sigma_ent,transpose(M))) + Bbar 
  else
    Sigma_exit = matmul(ele%mat6,matmul(Sigma_ent,transpose(ele%mat6)))
  endif
end subroutine

!+
! subroutine transport_with_sr_and_ibs(ele,M,Bbar,Sigma_ent,Sigma_exit)
!
! Transport a 6x6 beam envelope matrix through an element including
! damping and diffusion.
!
! Only bends include damping & diffusion.
!
! Input:
!   ele            -- integer: element key
!   M(6,6)         -- real(rp): transfer matrix with damping, populated by make_SR_mats
!   Bbar(6,6)      -- real(rp): diffusion matrix, populated by make_SR_mats
!   Sigma_ent(6,6) -- real(rp): incoming beam envelope matrix
! Output:
!   Sigma_exit(6,6) -- real(rp): outgoing beam envelope matrix
!-
subroutine transport_with_sr_and_ibs(ele,M,Bbar,Sigma_ent,Sigma_exit,tail_cut,tau_a,n_part)
  use bmad

  type(ele_struct) ele
  real(rp) Sigma_exit(6,6), Sigma_ent(6,6)
  real(rp) ibs_mat(6,6)
  real(rp) M(:,:), Bbar(:,:)
  real(rp) n_part, tau_a
  logical tail_cut

  call beam_envelope_ibs(Sigma_ent, ibs_mat, tail_cut, tau_a, ele%value(E_TOT$), n_part)

  if(any(ele%key==(/sbend$,rbend$/))) then
    Sigma_exit = matmul(M,matmul(Sigma_ent,transpose(M))) + Bbar + ibs_mat*ele%value(l$) 
  else
    Sigma_exit = matmul(ele%mat6,matmul(Sigma_ent,transpose(ele%mat6))) + ibs_mat*ele%value(l$)
  endif
end subroutine

!+
! subroutine make_V(M,V)
!
! For a one-turn transfer matrix M, this routine find the eigen matrix V.
! V is ordered such that its column pairs dominate the H,V, and then Z planes.
! It is normalized to be symplectic.
!
! Input:
!   M(6,6)     - real(rp): One turn transfer matrix.
! Output:
!   V(6,6)     - complex(rp): Matrix that diagonalizes M. 
!-
subroutine make_V(M,V)

  use mode3_mod, only: normalize_evecs, order_evecs_by_plane_dominance
  use f95_lapack

  implicit none

  real(rp) M(6,6)
  complex(rp) V(6,6)

  complex(rp) Vinv(6,6)

  real(rp) temp_mat(6,6), VR(6,6)
  real(rp) eval_r(6), eval_i(6)
  real(rp) evec_r(6,6), evec_i(6,6)
  real(rp) v1_r(6), v1_i(6)
  real(rp) v2_r(6), v2_i(6)
  real(rp) v3_r(6), v3_i(6)

  integer i
  integer i_error

  character(*), parameter :: r_name = 'make_V'

  temp_mat = M  !LA_GEEV destroys the contents of its first argument.
  CALL la_geev(temp_mat, eval_r, eval_i, VR=VR, INFO=i_error)
  if ( i_error /= 0 ) THEN
    call out_io (s_fatal$, r_name, "la_geev returned error: \i0\ ", i_error)
    if (global_com%exit_on_error) call err_exit
    eval_r = 0.0d0
    eval_i = 0.0d0
    evec_r = 0.0d0
    evec_i = 0.0d0
    return
  endif

  v1_r = VR(:, 1)  !v1_r
  v1_i = VR(:, 2)  !v1_i
  v2_r = VR(:, 3)  !v2_r
  v2_i = VR(:, 4)  !v2_i
  v3_r = VR(:, 5)  !v3_r
  v3_i = VR(:, 6)  !v3_i

  Vinv(:,1) = cmplx(v1_r,v1_i)
  Vinv(:,2) = (0.0d0,1.0d0)*conjg(Vinv(:,1))
  Vinv(:,3) = cmplx(v2_r,v2_i)
  Vinv(:,4) = (0.0d0,1.0d0)*conjg(Vinv(:,3))
  Vinv(:,5) = cmplx(v3_r,v3_i)
  Vinv(:,6) = (0.0d0,1.0d0)*conjg(Vinv(:,5))
  evec_r = real(Vinv)
  evec_i = aimag(Vinv)
  call order_evecs_by_plane_dominance(evec_r, evec_i, eval_r, eval_i)
  call normalize_evecs(evec_r, evec_i)
  Vinv = cmplx(evec_r,evec_i)
  V = mat_symp_conj_i(Vinv)   
end subroutine

!+
! subroutine envelope_radints(eles,co,alpha,emit)
!
! Calculates damping decrement and emittance of the three
! normal modes by integrating the diffusion and damping matrices.
! The damping times can obtained from alpha using:
!    tau = lattice_length/c_light/alpha
!
! alpha and emit are quantities for the three normal modes.
! alpha and emit are ordered by plane dominance.
!
! Only radiation from sbends and rbends is taken into account.
! The one-turn transfer matrix at each element (slice) is obtained
! by concatenating the individual element transfer matrices.
!
! Input:
!   eles(:)                      - ele_struct: array of element structures representing ring.
!          %mat6(6,6)            - real(rp): element transfer matrix.
!          %value(l$)            - real(rp): element (slice) length.
!          %value(rho$)          - real(rp): bending radius.
!          %value(b_field$)      - real(rp): bend field.
!          %value(b1_gradient$)  - real(rp): quadrupole field.
!   co(:)                        - coord_struct:  closed orbit at each element.
!        %vec(6)                 - real(rp): energy offset.
! Output:
!   alpha(3)                     - real(rp): Normal mode damping decrements.
!   emit(3)                      - real(rp): Normal mode emittances.
!-
subroutine envelope_radints(eles,co,alpha,emit)

  implicit none

  type(ele_struct) eles(:)
  type(coord_struct) co(:)
  real(rp) alpha(3), emit(3)

  real(rp) Dbar(6,6), Bbar(6,6)
  real(rp) one_turn_mat(6,6), temp_mat(6,6)
  real(rp) delta
  real(rp) gamma, l, rho, B0, B1

  complex(rp) Lambda(6,6), Theta(6,6)
  complex(rp) V(6,6)

  integer i, j

  one_turn_mat = I6
  do i=1,size(eles)
    one_turn_mat = matmul(eles(i)%mat6,one_turn_mat)
  enddo

  Lambda = (0.0d0,0.0d0)
  Theta = (0.0d0,0.0d0)
  do i=1,size(eles)
    if(any(eles(i)%key==(/sbend$,rbend$/))) then
      call convert_total_energy_to(eles(i)%value(E_TOT$), -1, gamma)
      delta = co(i)%vec(6)
      l = eles(i)%value(l$)
      rho = eles(i)%value(rho$)
      B0 = eles(i)%value(b_field$)
      B1 = eles(i)%value(b1_gradient$)
      call make_V(one_turn_mat,V)
      Lambda = Lambda + l*matmul(V,matmul(damping_matrix_D(gamma,rho,B0,B1,delta),mat_symp_conj_i(V)))
      Theta = Theta + l*matmul(V,matmul(diffusion_matrix_B(gamma,rho),conjg(transpose((V)))))
    endif
    one_turn_mat = matmul(eles(i)%mat6,matmul(one_turn_mat,mat_symp_conj(eles(i)%mat6)))
  enddo

  do i=1,3
    alpha(i) = real(Lambda(2*i-1,2*i-1))  !Eqn. 86
    emit(i) = real(Theta(2*i-1,2*i-1))/2.0d0/real(Lambda(2*i-1,2*i-1))  !Eqn. 91.  Theta should be real ... 
                                                                        !its imaginary part is due to rounding errors.
  enddo
end subroutine

!+
! subroutine make_PBRH(M, P, Bp, R, H)
!
! Decomposes the 1-turn transfer matrix into normal mode twiss-like parameters,
! according to Sec. IIIB of Ohmi, Hirata, and Oide paper.
!
! Note:  The Twiss parameters generated by this function are identical to those delivered
!        by mode3_mod.
!
! Input:
!   M(6,6)     -- real(rp): 1-turn transfer matrix
! Output:
!   P(6,6)     -- complex(rp):  Eqn. 97.  Phase advances.
!   Bp(6,6)    -- complex(rp):  Eqns. 89 & 101.  Beta functions.
!   R(6,6)     -- complex(rp):  Eqn. 99.  Transverse coupling.
!   H(6,6)     -- complex(rp):  Eqn. 100.  Longitudinal coupling.
!-
subroutine make_PBRH(M, P, Bp, R, H)

  use mode3_mod, only: normalize_evecs, order_evecs_by_plane_dominance
  use f95_lapack

  real(rp) M(6,6)
  complex(rp) P(6,6)
  complex(rp) Bp(6,6)
  complex(rp) R(6,6)
  complex(rp) H(6,6)
  
  complex(rp) V(6,6)
  complex(rp) Vinv(6,6)
  complex(rp) U(6,6), Up(6,6)
  complex(rp) Hx(2,2), Hy(2,2)
  complex(rp) R2(2,2)

  real(rp) a, b

  integer i

  character(*), parameter :: r_name = 'make_PBRH'

  call make_V(M,V)
  Vinv = mat_symp_conj_i(V)   
  P = matmul(V,matmul(M,Vinv))  ! Ohmi Eqn. 110
  U = matmul(mat_symp_conj_i(P),V)  ! Ohmi Eqn. 104
  a = sqrt(U(5,5)*U(6,6)-U(5,6)*U(6,5))  ! Ohmi Eqn. 105
  Hx = matmul(mat_symp_conj_i(U(5:6,1:2)),U(5:6,5:6))/a
  Hy = matmul(mat_symp_conj_i(U(5:6,3:4)),U(5:6,5:6))/a
  H(1:2,1:2) = (1.0d0-Hx(1,1)*Hx(2,2)+Hx(1,2)*Hx(2,1)/(1.0d0+a))*I2
  H(3:4,3:4) = (1.0d0-Hy(1,1)*Hy(2,2)+Hy(1,2)*Hy(2,1)/(1.0d0+a))*I2
  H(5:6,5:6) = a*I2
  H(1:2,3:4) = matmul(Hx,-mat_symp_conj_i(Hy))/(1.0d0+a)
  H(3:4,1:2) = matmul(Hy,-mat_symp_conj_i(Hx))/(1.0d0+a)
  H(1:2,5:6) = -Hx
  H(3:4,5:6) = -Hy
  H(5:6,1:2) = mat_symp_conj_i(Hx)
  H(5:6,3:4) = mat_symp_conj_i(Hy)
  Up = matmul(V,mat_symp_conj_i(H))
  b = sqrt(Up(3,3)*Up(4,4)-Up(3,4)*Up(4,3))
  R2 = matmul(mat_symp_conj_i(Up(3:4,3:4)),Up(3:4,1:2))/b
  R = 0.0d0
  R(1:2,1:2) = b*I2
  R(3:4,3:4) = b*I2
  R(5:6,5:6) = I2
  R(1:2,3:4) = -mat_symp_conj_i(R2)
  R(3:4,1:2) = R2
  Bp = matmul(Up,mat_symp_conj_i(R))
  ! write(*,*) "Bp = "
  ! do i=1,6
  !   write(*,'(6es14.5)') real(Bp(i,:))
  !   write(*,'(6es14.5)') aimag(Bp(i,:))
  !   write(*,*)
  ! enddo
  ! write(*,*) "Beta x  = ", Bp(2,2)**2*2.0d0, -2.0*(Bp(1,2)**2)
  ! write(*,*) "Alpha x = ", Bp(1,1)*Bp(2,2)*2.0d0-(1.0d0,0.0d0)
  ! write(*,*) "Alpha x = ", Bp(2,1)*Bp(2,2)*2.0d0+(0.0d0,1.0d0)
end subroutine

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! This is a sigma matrix based IBS calculation.
! It takes the beam sigma matrix and returns a matrix with changes to the 2nd order
! moments due to IBS.
!
! Use ibs_mat to change the sigma matrix like this:
! sigma_matrix_updated = sigma_matrix + ibs_mat*element_length
! See subroutine transport_with_sr_and_ibs in this module.
!
! Input:
!   sigma_mat(6,6)         - real(rp): beam sigma_matrix at element entrance
!   tail_cut               - logical:  If true, then apply tail cut to coulomb logarithm.
!   tau_a                  - real(rp): horizontal betatron damping rate.  Needed if tail_cut is true.
!   energy                 - real(rp): beam energy in eV
!   n_part                 - real(rp): number of particles in the bunch
!
! Output:
!   ibs_mat(6,6)           - real(rp): changes in 2nd order moments due to IBS are ibs_mat*element_length
!-

subroutine beam_envelope_ibs(sigma_mat, ibs_mat, tail_cut, tau_a, energy, n_part)
  use eigen_mod
  use LA_PRECISION, only: WP => DP
  use f95_lapack

  implicit none

  ! Ps.sigma_mat rearranges the sigma matrix so that all xx terms are in top left block,
  ! all px terms are in top right block, and pp terms are in bottom left block.  Lower right
  ! block is px'.
  integer, parameter :: Ps(6,6) = reshape( [1,0,0,0,0,0, 0,0,1,0,0,0, 0,0,0,0,1,0, &
                                            0,1,0,0,0,0, 0,0,0,1,0,0, 0,0,0,0,0,1],[6,6] )

  real(rp) sigma_mat(6,6)
  real(rp) ibs_mat(6,6)
  logical tail_cut
  real(rp) tau_a
  real(rp) n_part
  real(rp) energy

  real(rp) Tboost(6,6)
  real(rp) Tboost_inv(6,6)
  real(rp) Tspat(6,6)
  real(rp) Tspat_inv(6,6)
  real(rp) spatial_sigma_update(6,6)
  real(rp) sigma_update(6,6)
  real(rp) err_mat(6,6)
  real(rp) spatial_sigma_mat(6,6)
  real(rp) boosted_sigma_mat(6,6)
  real(rp) ar_sigma_mat(6,6)
  real(rp) sig_xx(3,3)
  real(rp) sig_xx_inv(3,3)
  real(rp) sig_xp(3,3)
  real(rp) sig_pp(3,3)
  real(rp) sig_pp_update(3,3)
  real(rp) sig_pl(3,3)
  real(rp) u(3)
  real(rp) R(3,3)
  integer etypes(3)
  real(rp) g1, g2, g3
  real(rp) clog
  real(rp) cI
  real(rp) Dw(3,3)
  real(rp) vol, vol1, pvol
  real(rp) ptrans, bn, bm, bmin, bmax
  real(rp) bmin1, bmin2

  integer error
  integer i, j, k

  real(rp), parameter :: pi_2 = pi/2.0d0

  !GSL for integrator
  type(fgsl_integration_workspace) :: integ_wk
  type(c_ptr) :: ptr
  type(fgsl_function) :: integrand_ready
  integer(fgsl_int) :: fgsl_status
  real(fgsl_double), target :: args(3)
  real(fgsl_double) :: abserr
  real(fgsl_double) :: integration_result

  real(rp) gamma

  logical ok

  call convert_total_energy_to(energy, -1, gamma)

  ! ! make transfer matrix from canonical to spatial coordinates
  ! Tspat = 0.0d0
  ! do i=1,6
  !   Tspat(i,i) = 1.0d0
  ! enddo
  ! !  Tspat(1,2) = element_length/4
  ! !  Tspat(3,4) = element_length/4
  ! !  Tspat(5,6) = element_length/4/gamma/gamma
  ! spatial_sigma_mat = matmul(Tspat,matmul(sigma_mat,transpose(Tspat)))
  spatial_sigma_mat = sigma_mat

  ! boost sigma matrix to COM frame of bunch
  Tboost = I6
  Tboost(5,5) = 1.0d0 + gamma*gamma/(1.0d0+gamma)
  Tboost(6,6) = 1.0d0 / gamma
  boosted_sigma_mat = matmul(Tboost,matmul(spatial_sigma_mat,transpose(Tboost)))

  ! permute sigma matrix to x,y,z,px,py,pz form
  ar_sigma_mat = matmul(matmul(transpose(Ps),boosted_sigma_mat),Ps)
  sig_xx = ar_sigma_mat(1:3,1:3)
  sig_xp = ar_sigma_mat(1:3,4:6)
  sig_pp = ar_sigma_mat(4:6,4:6)
  !call eigensys(sig_xx, u, R, etypes, 3, error)
  R=sig_xx  ! LA_SYEV destroys the contents of R

  call LA_SYEV(R,u,JOBZ='N')  !evals and evecs of symmetric real matrix
  vol1 = sqrt(u(1)*u(2)*u(3))
  vol = sqrt(4.0d0*pi)**3 * vol1
  bm = sqrt(min( u(1), u(2), u(3) ))  !minimum beam dimension
  call mat_inverse(sig_xx,sig_xx_inv,ok)
  if( .not. ok ) then
    write(*,*) "BAD: Could not invert sig_xx"
    ibs_mat = 0.0d0
    return
  endif

  !Get local momentum matrix
  sig_pl = sig_pp - matmul(transpose(sig_xp),matmul(sig_xx_inv,sig_xp))

  !Get eigen vectors of local momentum matrix
  call eigensys(sig_pl, u, R, etypes, 3, error)

  !R=sig_pl
  !call LA_SYEV(R,u,JOBZ='V',INFO=error)  !evals and evecs of symmetric real matrix
  !LA_SYEV seems to be less robust than eigensys.
  u(1) = max(u(1),1.0d-20)
  u(2) = max(u(2),1.0d-20)
  u(3) = max(u(3),1.0d-20)

  if( error .ne. 0 ) then
    write(*,'(A,I6," ",A)') "BAD: Eigenvectors of local momentum matrix not found."
    ibs_mat = 0.0d0
    return
  endif
  R=transpose(R) !R is defined as inverse of eigenvector matrix, and tr(R) = inv(r)
  ptrans = sqrt(u(1)+u(2)+u(3))
  pvol = sqrt(u(1)*u(2)*u(3))

  !- Integration using fgsl
  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(kubo_integrand, ptr)

  args = (/u(1),u(2),u(3)/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, pi_2, eps7, eps7, & 
                                     limit, 3, integ_wk, integration_result, abserr) 
  g1 = integration_result

  args = (/u(2),u(1),u(3)/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, pi_2, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  g2 = integration_result

  args = (/u(3),u(1),u(2)/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, pi_2, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  g3 = integration_result

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)
  !- End integration with fgsl

  bn = (vol/n_part)**(1.0d0/3.0d0)
  bmax = min( bm, bn )  !minimum dimension or debye radius
  if( tail_cut ) then
    !kubo's tail cut formula
    bmin1 = r_e/(ptrans*gamma)**2
    bmin2 = sqrt(abs(vol/n_part/pi/(ptrans*c_light)/tau_a))
    bmin = max( bmin1, bmin2 )
  else
    !no tail cut
    bmin = r_e/(ptrans*gamma)**2
  endif
  clog = log(bmax/bmin)

  cI = (r_e**2)*n_part*clog/4.0d0/pi/(gamma**4)/vol1/pvol   !extra factor of gamma because vol1 and pvol taken in rest frame

  Dw = 0.0d0
  Dw(1,1) = cI*(g2-g1+g3-g1)
  Dw(2,2) = cI*(g1-g2+g3-g2)
  Dw(3,3) = cI*(g1-g3+g2-g3)

  !- Build update matrix
  sig_pp_update = matmul(matmul(R,Dw),transpose(R))
  sigma_update = 0.0d0
  do i=1,3
    do j=1,3
      sigma_update(i*2,j*2) = sig_pp_update(i,j)
    enddo
  enddo

  ! boost updates to lab frame
  call mat_inverse(Tboost,Tboost_inv,ok)
  if( .not. ok ) then
    write(*,*) "BAD: Could not invert Tboost"
    ibs_mat = 0.0d0
    return
  endif
  spatial_sigma_update = matmul(Tboost_inv,matmul(sigma_update,transpose(Tboost_inv)))

  ! ! back to canonical coordinates
  ! call mat_inverse(Tspat, Tspat_inv, ok)
  ! sigma_update = matmul(Tspat_inv,matmul(spatial_sigma_update,transpose(Tspat_inv)))
  ibs_mat = spatial_sigma_update
end subroutine 

function kubo_integrand(s, params) bind(c)
    real(c_double), value :: s
    type(c_ptr), value :: params
    real(c_double) :: kubo_integrand

    real(c_double), pointer :: args(:)
    real(c_double) u1,u2,u3
    real(c_double) sins, sins2, coss, coss2

    call c_f_pointer(params,args,[3])
    u1 = args(1)
    u2 = args(2)
    u3 = args(3)

    sins = sin(s)
    coss = cos(s)
    sins2 = sins**2
    coss2 = coss**2

    kubo_integrand = (2.0d0*u1*sins2*coss) / &
                     sqrt((sins2 + u1/u2*coss2)*(sins2+u1/u3*coss2))
end function kubo_integrand

end module
























