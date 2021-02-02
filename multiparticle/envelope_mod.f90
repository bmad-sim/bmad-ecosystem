module envelope_mod

use twiss_and_track_mod
use fgsl
use, intrinsic :: iso_c_binding

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
! This subroutine is used to generate apply artificial vertical
! excitation to the beam envelope.
!
! Input:
!  ele             -- ele_struct: element to track through.
!     %mat6(6,6)   -- real(rp): element transter matrix
!     %value(l$)   -- real(rp): element length
! Output:
!  Yone(6,6)       -- real(rp): Sensitivity to envelope matrix to vertical kicks.
!-
subroutine make_Ykick_mat(ele,Yone)

type(ele_struct) ele
real(rp) Yone(:,:)

real(rp) M(6,6), Y(6,6), l

M = ele%mat6
l = ele%value(l$)

Y = 0.0d0
Y(4,4) = 1.0d0
Yone = l*matmul(M,matmul(Y,transpose(M)))

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine make_SR_mats(ele,co,M,Bone,Done)
!
! Equation references are to "From the beam-envelope matrix to synchrotron-radiation
! integrals" by K. Ohmi, K. Hirata, and K. Oide.  Phys.Rev.E v49n1, January 1994.
!
! This subroutine makes the synchrotron radiation damping and diffusion matrices.
! The damping is contained in M, and the diffusion is contained in Bone.
!
! See transport_with_sr_ele in this module for how to track through an element 
! using these matrices.
!
!
! Input:
!   ele                    -- ele_struct: 
!     %value(l$)           -- real(rp): element or slice length
!     %mat6(6,6)           -- real(rp): 6x6 transfer matrix
!     %value(b_field$)     -- real(rp): bend field
!     %value(b1_gradient$) -- real(rp): quadrupole field
!   co                     -- coord_struct:
!     %vec(6)              -- real(rp): p_z
!
! Output:
!   M(6,6)                 -- real(rp): transfer matrix with damping.
!   Bone(6,6)              -- real(rp): diffusion matrix B for one element.
!   Done(6,6)              -- real(rp): damping matrix D for one element.
!-
subroutine make_SR_mats(ele,co,M,Bone,Done)

type(ele_struct) ele
type(coord_struct) co
real(rp) M(:,:), Bone(:,:)
real(rp), optional :: Done(:,:)

real(rp) M0(6,6)
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

  M = M0 - l*matmul(M0,damping_matrix_D(gamma,rho,B0,B1,delta,co%species))  !Eqn. 80 & 81

  Bone = l*matmul(M,matmul(diffusion_matrix_B(gamma,rho,co%species),transpose(M)))  !Eqn. 88
  if(present(Done)) Done = l*matmul(M0,damping_matrix_D(gamma,rho,B0,B1,delta,co%species))
else
  M = M0
  Bone = 0.0d0
  if(present(Done)) Done = 0.0d0
endif
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Calculate diffusion matrix B from Ohmi Eqn 114 & 115.
!-
function diffusion_matrix_B(gamma,rho,species) result(mat)
real(rp) mat(6,6)
real(rp) gamma, rho
real(rp), parameter :: D_const = 55.0d0/24.0d0/sqrt(3.0d0)*(1.055e-34)/(9.11e-31)/c_light
integer species

!

mat = 0.0d0
mat(6,6) = D_const * classical_radius(species) * gamma**5 / abs(rho)**3

end function diffusion_matrix_B

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Calculate damping matrix D from Ohmi Eqn 121.
!
! This routine assumes zero scalar potential.
!-
function damping_matrix_D(gamma, rho, B0, B1, delta, species) result(mat)

real(rp) mat(6,6)
real(rp) gamma, rho, B0, B1, delta
real(rp) delta_f, energy_loss_rate
real(rp), parameter :: P_const = 2.0d0/3.0d0
integer species

!

mat = 0.0d0

delta_f = 1.0d0/(1.0d0+delta)

mat(2,2) = delta_f
mat(4,4) = delta_f
mat(6,6) = 2.0d0*delta_f

mat(6,1) = 1.0d0/rho + 2.0d0/B0*B1
mat(6,3) = 2.0d0/B0*B1

energy_loss_rate = P_const * classical_radius(species) * gamma**3 / rho**2 ! From Ohmi Eqn 43.
mat = energy_loss_rate * mat

end function damping_matrix_D

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine transport_with_sr(ele,M,Bone,Sigma_ent,Sigma_exit)
!
! Transport a 6x6 beam envelope matrix through an element including
! damping and diffusion.
!
! Only bends include damping & diffusion.
!
! Input:
!   ele            -- integer: element key
!   M(6,6)         -- real(rp): transfer matrix with damping, populated by make_SR_mats
!   Bone(6,6)      -- real(rp): diffusion matrix, populated by make_SR_mats
!   Sigma_ent(6,6) -- real(rp): incoming beam envelope matrix
! Output:
!   Sigma_exit(6,6) -- real(rp): outgoing beam envelope matrix
!-
subroutine transport_with_sr(ele,M,Bone,Yone,Sigma_ent,Sigma_exit)

type(ele_struct) ele
real(rp) Sigma_exit(6,6), Sigma_ent(6,6)
real(rp) M(:,:), Bone(:,:), Yone(:,:)

if(any(ele%key==(/sbend$,rbend$/))) then
  Sigma_exit = matmul(M,matmul(Sigma_ent,transpose(M))) + Bone + Yone
else
  Sigma_exit = matmul(ele%mat6,matmul(Sigma_ent,transpose(ele%mat6))) + Yone
endif
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine transport_with_sr_and_ibs(ele,M,Bone,Yone,Sigma_ent,Sigma_exit,tail_cut,tau,n_part, species)
!
! Transport a 6x6 beam envelope matrix through an element including
! SR damping, SR diffusion, and intrabeam scattering.
!
! Only bends contribute damping & diffusion.
! IBS contributes in all elements.
!
! Input:
!   ele            -- integer: element key
!   M(6,6)         -- real(rp): transfer matrix with damping, populated by make_SR_mats
!   Bone(6,6)      -- real(rp): diffusion matrix, populated by make_SR_mats
!   Yone(6,6)      -- real(rp): Vertical excitation matrix.  Used to introduce a source
!                               of vertical excitation.  Useful when simulating IBS
!                               in ideal lattices.
!   Sigma_ent(6,6) -- real(rp): incoming beam envelope matrix
!   tail_cut       -- logical: Apply tail cut.
!   tau            -- real(rp): Longest damping time of the three modes.  Used only
!                               if tail_cut is true.
!   n_part         -- real(rp): Number of particles in the bunch.
!   species        -- integer: Species of particle.
!
! Output:
!   Sigma_exit(6,6) -- real(rp): outgoing beam envelope matrix
!-
subroutine transport_with_sr_and_ibs(ele,M,Bone,Yone,Sigma_ent,Sigma_exit,tail_cut,tau,n_part, species)

type(ele_struct) ele
real(rp) Sigma_exit(6,6), Sigma_ent(6,6)
real(rp) ibs_mat(6,6)
real(rp) M(:,:), Bone(:,:), Yone(:,:)
real(rp) n_part, tau
integer species
logical tail_cut

!

call beam_envelope_ibs(Sigma_ent, ibs_mat, tail_cut, tau, ele%value(E_TOT$), n_part, species)

if(any(ele%key==(/sbend$,rbend$/))) then
  Sigma_exit = matmul(M,matmul(Sigma_ent,transpose(M))) + Bone + Yone + ibs_mat*ele%value(l$) 
else
  Sigma_exit = matmul(ele%mat6,matmul(Sigma_ent,transpose(ele%mat6))) + Yone + ibs_mat*ele%value(l$)
endif

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine make_V(M,V,abz_tunes)
!
! For a one-turn transfer matrix M, this routine find the eigen matrix V.
! V is ordered such that the per turn phase advance of its column pairs agree with abz_tunes.
! It is normalized to be symplectic.
!
! Input:
!   M(6,6)     - real(rp): One turn transfer matrix.
!   abz_tunes(3) - real(rp): per turn phase advance of each mode.
! Output:
!   V(6,6)     - complex(rp): Matrix that diagonalizes M. 
!-
subroutine make_V(M,V,abz_tunes)

use mode3_mod, only: normalize_evecs, order_evecs_by_tune, make_N
use f95_lapack
use sim_utils_interface

implicit none

real(rp) M(6,6)
complex(rp) V(6,6)
real(rp) abz_tunes(3)

complex(rp) Vinv(6,6)
complex(rp) temp_vec(6)
complex(rp) eval(6)

real(rp) temp_mat(6,6)
complex(rp) check_mat(6,6)
real(rp) VR(6,6)
real(rp) eval_r(6), eval_i(6)
real(rp) evec_r(6,6), evec_i(6,6)
real(rp) mat_tunes(3)

integer i
integer i_error
logical err_flag

character(*), parameter :: r_name = 'make_V'

temp_mat = M  !LA_GEEV destroys the contents of its first argument.

call la_geev(temp_mat, eval_r, eval_i, VR=VR, INFO=i_error)
eval = cmplx(eval_r, -eval_i, rp)
if ( i_error /= 0 ) THEN
  call out_io (s_fatal$, r_name, "la_geev returned error: \i0\ ", i_error)
  if (global_com%exit_on_error) call err_exit
  V = 0.0d0
  return
endif

Vinv(:,1) = cmplx(VR(:, 1),VR(:, 2), rp)
Vinv(:,2) = (0.0d0,1.0d0)*conjg(Vinv(:,1))
Vinv(:,3) = cmplx(VR(:, 3),VR(:, 4), rp)
Vinv(:,4) = (0.0d0,1.0d0)*conjg(Vinv(:,3))
Vinv(:,5) = cmplx(VR(:, 5),VR(:, 6), rp)
Vinv(:,6) = (0.0d0,1.0d0)*conjg(Vinv(:,5))
check_mat = matmul(transpose(Vinv),matmul(S6,conjg(Vinv))) !eqn. 79
if ( aimag(check_mat(1,1)) > 0.0 ) then
  temp_vec = Vinv(:,1)
  Vinv(:,1) = Vinv(:,2)
  Vinv(:,2) = temp_vec
  eval(1) = conjg(eval(1))
  eval(2) = conjg(eval(2))
endif
if ( aimag(check_mat(3,3)) > 0.0 ) then
  temp_vec = Vinv(:,3)
  Vinv(:,3) = Vinv(:,4)
  Vinv(:,4) = temp_vec
  eval(3) = conjg(eval(3))
  eval(4) = conjg(eval(4))
endif
if ( aimag(check_mat(5,5)) > 0.0 ) then
  temp_vec = Vinv(:,5)
  Vinv(:,5) = Vinv(:,6)
  Vinv(:,6) = temp_vec
  eval(5) = conjg(eval(5))
  eval(6) = conjg(eval(6))
endif

mat_tunes(1) = MyTan(aimag(eval(1)), real(eval(1)))
mat_tunes(2) = MyTan(aimag(eval(3)), real(eval(3)))
mat_tunes(3) = MyTan(aimag(eval(5)), real(eval(5)))

call order_evecs_by_tune(Vinv, eval, mat_tunes, abz_tunes, err_flag)
if (err_flag) then
  call out_io (s_fatal$, r_name, "order_evecs_by_tune failed to identify eigen modes.")
  write(*,'(a,3f14.5)') "Tunes supplied to subroutine (abz_tunes):        ", abz_tunes
  write(*,'(a,3f14.5)') "Tunes obtained from one-turn matrix (mat_tunes): ", mat_tunes
  if (global_com%exit_on_error) call err_exit
  V = 0.0d0
  return
endif

call normalize_evecs(Vinv, err_flag)
if (err_flag) then
  call out_io (s_fatal$, r_name, "Zero amplitude eigenvectors.")
  write(*,'(a,3f14.5)') "Tunes supplied to subroutine (abz_tunes):        ", abz_tunes
  write(*,'(a,3f14.5)') "Tunes obtained from one-turn matrix (mat_tunes): ", mat_tunes
  if (global_com%exit_on_error) call err_exit
  V = 0.0d0
  return
endif


V = mat_symp_conj_i(Vinv)

!-----------------------------------------
contains

function MyTan(y, x) result(arg)
  !For a complex number x+iy graphed on an xhat, yhat plane, this routine returns the angle
  !between (x, y) and the +x axis, measured counter-clockwise.  There is a branch cut along +x.
  !This routine returns a number between 0 and 2pi.
  real(rp) x, y, arg

  arg = atan2(y,x)
  if (arg .lt. 0) then
    arg = arg + twopi
  endif
  arg = twopi - arg

end function MyTan

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine integrated_mats(eles,coos,Lambda,Theta,Iota,mode)
!
! Input:
!   eles(:)                      - ele_struct: array of element structures representing ring.
!          %mat6(6,6)            - real(rp): element transfer matrix.
!          %value(l$)            - real(rp): element (slice) length.
!          %value(rho$)          - real(rp): bending radius.
!          %value(b_field$)      - real(rp): bend field.
!          %value(b1_gradient$)  - real(rp): quadrupole field.
!   coos(:)                      - coord_struct:  closed orbit at each element.
!        %vec(6)                 - real(rp): energy offset.
!   mode                         - normal_modes_struct
!       %a%tune                  - real(rp): a-mode tune
!       %b%tune                  - real(rp): b-mode tune
!       %z%tune                  - real(rp): z-mode tune (aka c-mode tune)
! Output:
!   Lambda(6,6)                  - real(rp): Integrated SR damping matrix.
!   Theta(6,6)                   - real(rp): Integrated SR diffusion matrix.
!   Iota(6,6)                    - real(rp): Integrated artificial vertical excitation matrix.  Used
!                                            to generate vertical emittance in ideal lattices.
!-
subroutine integrated_mats(eles,coos,Lambda,Theta,Iota,mode)

type(ele_struct) eles(:)
type(coord_struct) coos(:)
complex(rp) Lambda(6,6)  !integrated damping matrix
complex(rp) Theta(6,6)  !integrated diffusion matrix
complex(rp) Iota(6,6)  !integrated vertical excitation matrix
type(normal_modes_struct) mode
 
real(rp) one_turn_mat(6,6)
real(rp) delta
real(rp) gamma, l, rho, B0, B1
real(rp) Y(6,6)
real(rp) abz_tunes(3)

complex(rp) V(6,6), Vinv(6,6)

logical err_flag

integer i, j

one_turn_mat = I6
do i=1,size(eles)
  one_turn_mat = matmul(eles(i)%mat6,one_turn_mat)
enddo

Y = 0.0d0
Y(4,4) = 1.0d0

abz_tunes(1) = mode%a%tune
abz_tunes(2) = mode%b%tune
abz_tunes(3) = mode%z%tune

Lambda = (0.0d0,0.0d0)
Theta = (0.0d0,0.0d0)
Iota = (0.0d0,0.0d0)
do i=1,size(eles)
  l = eles(i)%value(l$)
  one_turn_mat = matmul(eles(i)%mat6,matmul(one_turn_mat,mat_symp_conj(eles(i)%mat6)))
  call make_V(one_turn_mat,V,abz_tunes)
  Vinv = mat_symp_conj_i(V) 
  if(any(eles(i)%key==(/sbend$,rbend$/))) then
    call convert_total_energy_to(eles(i)%value(E_TOT$), -1, gamma)
    delta = coos(i)%vec(6)
    rho = eles(i)%value(rho$)
    B0 = eles(i)%value(b_field$)
    B1 = eles(i)%value(b1_gradient$)
    Lambda = Lambda + l*matmul(V,matmul(damping_matrix_D(gamma,rho,B0,B1,delta, coos(i)%species),Vinv))
    Theta = Theta + l*matmul(V,matmul(diffusion_matrix_B(gamma,rho,coos(i)%species),conjg(transpose((V)))))
  endif
  Iota = Iota + l*matmul(V,matmul(Y,conjg(transpose((V)))))
enddo
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine envelope_radints_ibs(Lambda,Theta,Iota,eles,alpha,emit,mode,tail_cut,npart, species)
!
! Calculates damping decrement and emittance of the three
! normal modes by integrating the IBS, SR diffusion, and SR damping matrices.
!
! The IBS depends on the envelope, and so this routine iterates to
! locate the equilibrium beam envelope. This iterative process can fail to converge.
!
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
!   Lambda(6,6)                  -- real(rp): Integrated damping matrix.
!   Theta(6,6)                   -- real(rp): Integrated diffusion matrix.
!   Iota(6,6)                    -- real(rp): Integrated vertical excitation matrix.
!   eles(:)                      -- ele_struct: array of element structures representing ring.
!          %mat6(6,6)            -- real(rp): element transfer matrix.
!          %value(l$)            -- real(rp): element (slice) length.
!          %value(E_TOT$)        -- real(rp): Beam energy in element.
!   mode                         -- normal_modes_struct
!       %a%tune                  -- real(rp): tune of a-mode.
!       %b%tune                  -- real(rp): tune of b-mode.
!       %z%tune                  -- real(rp): tune of z-mode.
!   tail_cut                     -- logical: apply tail cut.
!   npart                        -- real(rp): number of particles in bunch.
!   species                      -- integer: Particle species.
!
! Output:
!   alpha(3)                     -- real(rp): Normal mode damping decrements.
!   emit(3)                      -- real(rp): Normal mode emittances.
!-
subroutine envelope_radints_ibs(Lambda,Theta,Iota,eles,alpha,emit,mode,tail_cut,npart, species)

use mode3_mod

implicit none

type(ele_struct) eles(:)
type(normal_modes_struct) mode

complex(rp) Lambda(6,6), Theta(6,6), Iota(6,6), Omega(6,6)
complex(rp) V(6,6)

real(rp) alpha(3), emit(3)
real(rp) npart
real(rp) abz_tunes(3)
real(rp) tau(3), tau_max
real(rp) emit_new(3)
real(rp) one_turn_mat(6,6), one_turn_mat0(6,6), sigma_mat(6,6), temp_mat(6,6)
real(rp) ibs_C(6,6)
real(rp) energy, res
real(rp) l

integer species
integer i, j

logical tail_cut
logical err_flag

abz_tunes(1) = mode%a%tune
abz_tunes(2) = mode%b%tune
abz_tunes(3) = mode%z%tune

one_turn_mat0 = I6
do i=1,size(eles)
  one_turn_mat0 = matmul(eles(i)%mat6,one_turn_mat0)
enddo

alpha(1) = real(Lambda(1,1))  !Eqn. 86
alpha(2) = real(Lambda(3,3))  
alpha(3) = real(Lambda(5,5))  
emit(1) = (real(Iota(1,1))+real(Theta(1,1)))/2.0d0/real(Lambda(1,1)) 
emit(2) = (real(Iota(3,3))+real(Theta(3,3)))/2.0d0/real(Lambda(3,3))  
emit(3) = (real(Iota(5,5))+real(Theta(5,5)))/2.0d0/real(Lambda(5,5))  

!alpha = abs(alpha)
!emit = abs(emit)

tau_max = maxval(eles(size(eles))%s / c_light / alpha)
mode%a%emittance = emit(1)
mode%b%emittance = emit(2)
mode%z%emittance = emit(3)
do while(.true.)
  !build integrated IBS mat using new emittances
  one_turn_mat = one_turn_mat0
  Omega = (0.0d0,0.0d0)  !Integrated IBS matrix
  do i=1,size(eles)
    call make_smat_from_abc(one_turn_mat, mode, sigma_mat, err_flag)
    l = eles(i)%value(l$)
    energy = eles(i)%value(E_TOT$)
    call make_V(one_turn_mat,V,abz_tunes)
    ibs_C = ibs_matrix_C(sigma_mat,tail_cut,tau_max,energy,npart,species)
    Omega = Omega + l*matmul(V,matmul(ibs_C,conjg(transpose((V)))))
    one_turn_mat = matmul(eles(i)%mat6,matmul(one_turn_mat,mat_symp_conj(eles(i)%mat6)))
  enddo

  !calculate new emittances
  emit_new(1) = (real(Theta(1,1))+real(Iota(1,1))+real(Omega(1,1)))/2.0d0/real(Lambda(1,1))
  emit_new(2) = (real(Theta(3,3))+real(Iota(3,3))+real(Omega(3,3)))/2.0d0/real(Lambda(3,3))
  emit_new(3) = (real(Theta(5,5))+real(Iota(5,5))+real(Omega(5,5)))/2.0d0/real(Lambda(5,5))
  res = sqrt( ((emit_new(1)-mode%a%emittance)/mode%a%emittance)**2 + ((emit_new(2)-mode%b%emittance)/mode%b%emittance)**2 + ((emit_new(3)-mode%z%emittance)/mode%z%emittance)**2 )
  mode%a%emittance = emit_new(1)
  mode%b%emittance = emit_new(2)
  mode%z%emittance = emit_new(3)
  if(res .lt. 0.00001) then
    exit
  endif
enddo
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine envelope_radints(eles,co,alpha,emit)
!
! Calculates damping decrement and emittance of the three
! normal modes from the integrate diffusion, damping, and vertical
! excitation matrices names Lambda, Theta, and Iota, respectively.
! These three matrices are obtained from the subroutine integrated_mats.
!
! The damping times can obtained from alpha using:
!    tau = lattice_length/c_light/alpha
!
! Input:
!   Lambda(6,6)                  - complex(rp): integrated damping matrix.
!   Theta(6,6)                   - complex(rp): integrated diffusion (sr excitation) matrix.
!   Iota(6,6)                    - complex(rp): integrated vertical excitation matrix.
!
! Output:
!   alpha(3)                     - real(rp): Normal mode damping decrements.
!   emit(3)                      - real(rp): Normal mode emittances.
!-
subroutine envelope_radints(Lambda,Theta,Iota,alpha,emit)
implicit none

complex(rp) Lambda(6,6), Theta(6,6), Iota(6,6)
real(rp) alpha(3), emit(3)

alpha(1) = real(Lambda(1,1))  !Eqn. 86
alpha(2) = real(Lambda(3,3))
alpha(3) = real(Lambda(5,5))
emit(1) = (real(Iota(1,1))+real(Theta(1,1)))/2.0d0/real(Lambda(1,1))  !Eqn.
emit(2) = (real(Iota(3,3))+real(Theta(3,3)))/2.0d0/real(Lambda(3,3))
emit(3) = (real(Iota(5,5))+real(Theta(5,5)))/2.0d0/real(Lambda(5,5))

!alpha = abs(alpha)
!emit = abs(emit)
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! subroutine make_PBRH(M, P, Bp, R, H, abz_tunes)
!
! Decomposes the 1-turn transfer matrix into normal mode twiss-like parameters,
! according to Sec. IIIB of Ohmi, Hirata, and Oide paper.
!
! Note:  The Twiss parameters generated by this function are identical to those delivered
!        by mode3_mod.
!
! Input:
!   M(6,6)     -- real(rp): 1-turn transfer matrix
!   abz_tunes(3) -- real(rp): tunes for a,b, and c modes.  Used to identify which eigenvector
!                             is associated with which mode.
! Output:
!   P(6,6)     -- complex(rp):  Eqn. 97.  Phase advances.
!   Bp(6,6)    -- complex(rp):  Eqns. 89 & 101.  Beta functions.
!   R(6,6)     -- complex(rp):  Eqn. 99.  Transverse coupling.
!   H(6,6)     -- complex(rp):  Eqn. 100.  Longitudinal coupling.
!-
subroutine make_PBRH(M, P, Bp, R, H, abz_tunes)

use mode3_mod, only: normalize_evecs, order_evecs_by_plane_dominance
use f95_lapack

real(rp) M(6,6)
complex(rp) P(6,6)
complex(rp) Bp(6,6)
complex(rp) R(6,6)
complex(rp) H(6,6)
real(rp) abz_tunes(3)

complex(rp) V(6,6)
complex(rp) Vinv(6,6)
complex(rp) U(6,6), Up(6,6)
complex(rp) Hx(2,2), Hy(2,2)
complex(rp) R2(2,2)

real(rp) a, b

integer i

character(*), parameter :: r_name = 'make_PBRH'

call make_V(M,V,abz_tunes)
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

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! This is a function wrapper for the beam_envelope_ibs subroutine.
! See that subroutine for documentation.
!
! Input:
!   sigma_mat(6,6)         -- real(rp): beam sigma_matrix at element entrance
!   tail_cut               -- logical:  If true, then apply tail cut to coulomb logarithm.
!   tau                    -- real(rp): horizontal betatron damping rate.  Needed if tail_cut is true.
!   energy                 -- real(rp): beam energy in eV
!   n_part                 -- real(rp): number of particles in the bunch
!   species                -- integer: Partical species.
!
! Output:
!   ibs_mat(6,6)           -- real(rp): changes in 2nd order moments due to IBS are ibs_mat*element_length
!-
function ibs_matrix_C(sigma_mat, tail_cut, tau, energy, n_part, species) result(ibs_mat)

logical tail_cut
real(rp) sigma_mat(6,6), ibs_mat(6,6)
real(rp) tau, energy, n_part
integer species

call beam_envelope_ibs(sigma_mat, ibs_mat, tail_cut, tau, energy, n_part, species)

end function ibs_matrix_C

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! Subroutine beam_envelope_ibs(sigma_mat, ibs_mat, tail_cut, tau, energy, n_part, species)
!
! This is a sigma matrix based IBS calculation.
! It takes the beam sigma matrix and returns a matrix with changes to the 2nd order
! moments due to IBS.
!
! Use ibs_mat to change the sigma matrix like this:
! sigma_matrix_updated = sigma_matrix + ibs_mat*element_length
! See subroutine transport_with_sr_and_ibs in this module.
!
! Input:
!   sigma_mat(6,6)         -- real(rp): beam sigma_matrix at element entrance
!   tail_cut               -- logical:  If true, then apply tail cut to coulomb logarithm.
!   tau                    -- real(rp): horizontal betatron damping rate.  Needed if tail_cut is true.
!   energy                 -- real(rp): beam energy in eV
!   n_part                 -- real(rp): number of particles in the bunch
!   species                -- integer: Partical species.
!
! Output:
!   ibs_mat(6,6)           -- real(rp): changes in 2nd order moments due to IBS are ibs_mat*element_length
!-

subroutine beam_envelope_ibs(sigma_mat, ibs_mat, tail_cut, tau, energy, n_part, species)

use eigen_mod
use LA_PRECISION, only: WP => DP
use f95_lapack

implicit none

! Ps.sigma_mat rearranges the sigma matrix so that all xx terms are in top left block,
! all px terms are in top right block, and pp terms are in bottom left block.  Lower right
! block is px'.
integer, parameter :: Ps(6,6) = reshape( [1,0,0,0,0,0, 0,0,1,0,0,0, 0,0,0,0,1,0, &
                                          0,1,0,0,0,0, 0,0,0,1,0,0, 0,0,0,0,0,1],[6,6] )

real(rp) sigma_mat(6,6), ibs_mat(6,6)
real(rp) tau, n_part, energy
real(rp) Tboost(6,6), Tboost_inv(6,6), Tspat(6,6), Tspat_inv(6,6)
real(rp) spatial_sigma_update(6,6), sigma_update(6,6), err_mat(6,6)
real(rp) spatial_sigma_mat(6,6), boosted_sigma_mat(6,6)
real(rp) ar_sigma_mat(6,6), sig_xx(3,3), sig_xx_inv(3,3), sig_xp(3,3)
real(rp) sig_pp(3,3), sig_pp_update(3,3), sig_pl(3,3)
real(rp) u(3), R(3,3), g1, g2, g3, clog, cI
real(rp) Dw(3,3), vol, vol1, pvol, ptrans, bn, bm, bmin, bmax
real(rp) bmin1, bmin2, gamma
real(rp), parameter :: pi_2 = pi/2.0d0

integer species, etypes(3)
integer i, j, k

logical tail_cut

!GSL for integrator
type(fgsl_integration_workspace) :: integ_wk
type(c_ptr) :: ptr
type(fgsl_function) :: integrand_ready
integer(fgsl_int) :: fgsl_status
real(fgsl_double), target :: args(3)
real(fgsl_double) :: abserr
real(fgsl_double) :: integration_result

logical ok, error

!
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
!call eigensys(sig_xx, u, R, etypes, 3, error)
sig_xx = ar_sigma_mat(1:3,1:3)
sig_xp = ar_sigma_mat(1:3,4:6)
sig_pp = ar_sigma_mat(4:6,4:6)
R=sig_xx  ! LA_SYEV destroys the contents of R

call LA_SYEV(R,u,JOBZ='N')  !evals and evecs of symmetric real matrix
u(1) = max(u(1),1.0d-20)
u(2) = max(u(2),1.0d-20)
u(3) = max(u(3),1.0d-20)
if( (sqrt(u(1)/u(2)) .gt. 1.0e10) .or. \
    (sqrt(u(2)/u(1)) .gt. 1.0e10) .or. \
    (sqrt(u(1)/u(3)) .gt. 1.0e10) .or. \
    (sqrt(u(3)/u(1)) .gt. 1.0e10) .or. \
    (sqrt(u(2)/u(3)) .gt. 1.0e10) .or. \
    (sqrt(u(3)/u(2)) .gt. 1.0e10) ) then
    write(*,'(a)') "Warning from beam_envelope_ibs: physical beam dimensions differ by more than 10 orders of magnitude."
    write(*,'(a)') "Something could be wrong, check emittances."
    write(*,'(a,3es13.3)') "sqrt(u(1)), sqrt(u(2)), sqrt(u(3)): ", sqrt(u(1)), sqrt(u(2)), sqrt(u(3))
    if (global_com%exit_on_error) call err_exit
endif
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

if( error) then
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
  bmin1 = classical_radius(species)/(ptrans*gamma)**2
  bmin2 = sqrt(abs(vol/n_part/pi/(ptrans*c_light)/tau))
  bmin = max( bmin1, bmin2 )
else
  !no tail cut
  bmin = classical_radius(species)/(ptrans*gamma)**2
endif
clog = log(bmax/bmin)

cI = (classical_radius(species)**2)*n_part*clog/4.0d0/pi/(gamma**4)/vol1/pvol   !extra factor of gamma because vol1 and pvol taken in rest frame

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

end subroutine beam_envelope_ibs

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------

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

kubo_integrand = (2.0d0*u1*sins2*coss) / sqrt((sins2 + u1/u2*coss2)*(sins2+u1/u3*coss2))
end function kubo_integrand

end module
























