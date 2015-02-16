!+
! Subroutine twiss_propagate1 (ele1, ele2, err_flag)
!
! Subroutine to propagate the twiss, coupling, and dispersion parameters from 
! the exit end of ele1 to the exit end of ele2.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele1        -- Ele_struct: Structure holding the starting parameters.
!   ele2        -- Ele_struct: Structure holding the transfer matrix.
!     %key                -- Needed since, for example, Match element are handled 
!                              differently from other elements.
!     %map_ref_orb_in     -- Important for the dispersion calc.
!     %map_ref_orb_out    -- Important for the dispersion calc.
!   bmad_status -- Common block status structure:
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   ele2     -- Ele_struct: Structure for the ending parameters.
!   err_flag -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine twiss_propagate1 (ele1, ele2, err_flag)

use bmad_interface, except_dummy => twiss_propagate1

implicit none
type (ele_struct), target :: ele1, ele2
type (twiss_struct) twiss_a
type (lat_param_struct) param
type (coord_struct), pointer :: orb, orb_out
integer key2

real(rp), pointer :: mat6(:,:)
real(rp) v_mat(4,4), v_inv_mat(4,4), y_inv(2,2), det, mat2_a(2,2), mat2_b(2,2)
real(rp) big_M(2,2), small_m(2,2), big_N(2,2), small_n(2,2)
real(rp) c_conj_mat(2,2), E_inv_mat(2,2), F_inv_mat(2,2)
real(rp) mat2(2,2), eta1_vec(6), eta_vec(6), dpz2_dpz1, rel_p1, rel_p2
real(rp) det_factor, deriv_rel, gamma2_c

logical error
logical, optional :: err_flag

character(20), parameter :: r_name = 'twiss_propagate1'

!---------------------------------------------------------------------
! init

if (present(err_flag)) err_flag = .true.

if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then

  ! For x-ray lines assume that since beta = 0 there is no interest in calculating the Twiss parameters.
  ! So this is not treated as an error.
  if (associated(ele1%branch)) then
    if (ele1%branch%param%particle == photon$) then
      if (present(err_flag)) err_flag = .false.
      return
    endif
  endif

  call out_io (s_fatal$, r_name, 'ZERO BETA DETECTED AT: ' // trim(ele1%name), &
                                 'ELEMENT # \i0\ ',  i_array = [ele1%ix_ele])
  if (global_com%exit_on_error) call err_exit
  return
endif

ele2%mode_flip = ele1%mode_flip          ! assume no flip
key2 = ele2%key

!---------------------------------------------------------------------
! Special match case

if (ele2%key == match$ .and. ele2%value(match_end$) /= 0) then
  ele2%value(beta_a0$)  = ele1%a%beta
  ele2%value(beta_b0$)  = ele1%b%beta
  ele2%value(alpha_a0$) = ele1%a%alpha
  ele2%value(alpha_b0$) = ele1%b%alpha
  ele2%value(eta_x0$)   = ele1%x%eta
  ele2%value(eta_y0$)   = ele1%y%eta
  ele2%value(etap_x0$)  = ele1%x%etap
  ele2%value(etap_y0$)  = ele1%y%etap
  ele2%value(c_11$:c_22$) = [ele1%c_mat(1,1), ele1%c_mat(1,2), ele1%c_mat(2,1), ele1%c_mat(2,2)]
  ele2%value(gamma_c$) = ele1%gamma_c
  call make_mat6 (ele2, param)
endif

!---------------------------------------------------------------------
! markers are easy

if (key2 == marker$ .or. key2 == photon_fork$ .or. key2 == fork$) then
  call transfer_twiss (ele1, ele2)
  if (present(err_flag)) err_flag = .false.
  return
endif

!---------------------------------------------------------------------
! if transfer matrix is not coupled...
! propagate c_mat coupling matrix and setup temporary element for propagation

mat6 => ele2%mat6

if (all(ele2%mat6(1:2,3:4) == 0)) then

  mat2_a = ele2%mat6(1:2,1:2)
  mat2_b = ele2%mat6(3:4,3:4) 

  call mat_symp_conj (mat2_b, y_inv) ! conj == inverse
  ele2%c_mat = matmul(matmul(mat2_a, ele1%c_mat), y_inv)
  ele2%gamma_c = ele1%gamma_c

!---------------------------------------------------------------------
! here if we are dealing with a coupled transfer matrix

else

  ! det_factor is a renormalization factor since det_original != 1

  det_factor = 1
  if (ele2%value(p0c$) /= ele2%value(p0c_start$)) then
    det_factor = sqrt(determinant (mat6(1:4,1:4)))
  endif

  big_M = mat6(1:2,1:2)
  small_m = mat6(1:2,3:4)
  big_N = mat6(3:4,3:4)
  small_n = mat6(3:4,1:2)

  call mat_symp_conj (ele1%c_mat, c_conj_mat)
  mat2 = ele1%gamma_c * big_M - matmul(small_m, c_conj_mat)
  det = determinant(mat2) / det_factor 

  ! we demand that gamma_c > 0.3 (ie det > 0.1)
  ! try to make it so that there is no net mode flip here

  if (det > 0.9 .or. (det > 0.1 .and. .not. ele1%mode_flip)) then

    gamma2_c = sqrt(det)
    mat2_a = mat2 / gamma2_c
    mat2_b = (ele1%gamma_c * big_N + matmul(small_n, ele1%c_mat)) / gamma2_c
    call mat_symp_conj (mat2_b, F_inv_mat)
    ele2%c_mat = matmul(matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m, F_inv_mat)
    ele2%gamma_c = sqrt(det)

  ! else we flip the modes

  else

    mat2 = matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m
    det = determinant(mat2) / det_factor
    if (det < 0) then
      call out_io (s_error$, r_name,  '||mat2|| < 0! (Due to roundoff?) \f10.2\ ', &
                                      'When propagating through: [\i0\]  ' // trim(ele2%name), &
                                      r_array = [det], i_array = [ele2%ix_ele])
    endif

    gamma2_c = sqrt(abs(det))
    mat2_a = (ele1%gamma_c * small_n - matmul(big_N, c_conj_mat)) / gamma2_c
    mat2_b = mat2 / gamma2_c

    call mat_symp_conj (mat2_a, E_inv_mat)
    ele2%c_mat = matmul(ele1%gamma_c * big_M - matmul(small_m, c_conj_mat), E_inv_mat)
    ele2%gamma_c = gamma2_c
    ele2%mode_flip = .not. ele1%mode_flip

  endif

endif

!---------------------------------------------------------------------
! Propagate twiss.
! If there is a mode flip, ele%a is the "b"-mode. That is, ele%a is associated with
! the lower right block of the U matrix (See Sagan & Rubin: Linear Analysis of Coupled Lattices).
! Another way of saying this: ele%a always represents the same physical mode.

if (ele1%mode_flip) then
  call twiss1_propagate (ele1%a, mat2_b, ele2%key, ele2%value(l$), ele2%a, error); if (error) return
  call twiss1_propagate (ele1%b, mat2_a, ele2%key, ele2%value(l$), ele2%b, error); if (error) return
else
  call twiss1_propagate (ele1%a, mat2_a, ele2%key, ele2%value(l$), ele2%a, error); if (error) return
  call twiss1_propagate (ele1%b, mat2_b, ele2%key, ele2%value(l$), ele2%b, error); if (error) return
endif

! Comming out of a flipped state, the calculation is often off by a factor of twopi.
! The code corrects this. However, since there is no proof that this always happens, 
! this is a bit of a kludge. Factors of twopi are not physically meaningful so this
! does not affect any calculations.

if (ele1%mode_flip .and. .not. ele2%mode_flip) then
  ele2%a%phi = ele2%a%phi - twopi
  ele2%b%phi = ele2%b%phi - twopi
endif

!----------------------------------------------------
! Dispersion calc.
! p_z2 is p_z at end of ele2 assuming p_z = 1 at end of ele1.
! This is just 1.0 (except for RF cavities).

orb  => ele2%map_ref_orb_in
orb_out => ele2%map_ref_orb_out
rel_p1 = 1 + orb%vec(6)               ! reference energy 
rel_p2 = 1 + orb_out%vec(6)

eta1_vec = [ele1%x%eta, ele1%x%etap * rel_p1, ele1%y%eta, ele1%y%etap * rel_p1, ele1%z%eta, 1.0_rp]

! For a circular ring, defining the dependence of z on pz is problematical.
! With the RF off, z is not periodic so dz/dpz is dependent upon turn number.
! With the RF on, dz/dpz, along with the other dispersion components, is not well defined.
! This being the case, eta1_vec(5) is just treated as zero for an rfcavity.

if (key2 == rfcavity$) eta1_vec(5) = 0

! Must avoid 0/0 divide at zero reference momentum. 
! If rel_p1 = 0 then total momentum is zero and orb%vec(2) and orb%vec(4) must be zero.

dpz2_dpz1 = dot_product(mat6(6,:), eta1_vec) 

if (rel_p1 == 0) then
  eta_vec(1:5) = matmul (ele2%mat6(1:5,:), eta1_vec) / dpz2_dpz1
else
  dpz2_dpz1 = dpz2_dpz1 + (mat6(6,2) * orb%vec(2) + mat6(6,4) * orb%vec(4)) / rel_p1
  deriv_rel = dpz2_dpz1 * rel_p1
  eta_vec(1) = (mat6(1,2) * orb%vec(2) + mat6(1,4) * orb%vec(4)) / deriv_rel
  eta_vec(2) = (mat6(2,2) * orb%vec(2) + mat6(2,4) * orb%vec(4)) / deriv_rel - orb_out%vec(2) / rel_p2
  eta_vec(3) = (mat6(3,2) * orb%vec(2) + mat6(3,4) * orb%vec(4)) / deriv_rel
  eta_vec(4) = (mat6(4,2) * orb%vec(2) + mat6(4,4) * orb%vec(4)) / deriv_rel - orb_out%vec(4) / rel_p2
  eta_vec(5) = (mat6(5,2) * orb%vec(2) + mat6(5,4) * orb%vec(4)) / deriv_rel
  eta_vec(1:5) = eta_vec(1:5) + matmul (ele2%mat6(1:5,:), eta1_vec) / dpz2_dpz1
endif

eta_vec(2) = eta_vec(2) / rel_p2
eta_vec(4) = eta_vec(4) / rel_p2

ele2%x%eta  = eta_vec(1)
ele2%x%etap = eta_vec(2)
ele2%y%eta  = eta_vec(3)
ele2%y%etap = eta_vec(4)
ele2%z%eta  = eta_vec(5)
ele2%z%etap = ele1%z%etap * dpz2_dpz1

call make_v_mats (ele2, v_mat, v_inv_mat)
eta_vec(1:4) = matmul (v_inv_mat, eta_vec(1:4))

ele2%a%eta  = eta_vec(1)
ele2%a%etap = eta_vec(2)
ele2%b%eta  = eta_vec(3)
ele2%b%etap = eta_vec(4)

if (present(err_flag)) err_flag = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine twiss1_propagate (twiss1, mat2, ele_key, length, twiss2, err)
!
! Subroutine to propagate the twiss parameters of a single mode.
!
! The betatron phase phi is only determined up to a factor of 2pi. 
! The length and ele_key argument is used to determine what phi should be:
!   length > 0:                           0 <= phi <  twopi
!   length = 0 or ele_key = patch$:     -pi <  phi <= pi
!   length < 0:                      -twopi <  phi <= 0
! The patch element is exceptional in that its length is defined in a somewhat 
! arbitrary manner and thus is not a good reference as to what the phase advance should be.
!
! Modules needed:
!   use bmad
!
! Input:
!   twiss1    -- Twiss_struct: Input Twiss parameters.
!   mat2(2,2) -- Real(rp): The transfer matrix.
!   ele_key   -- Integer: quadrupole$, etc.
!   length    -- Real(rp): Determines whether the phase is 
!                            increasing or decreasing.
!
! Output:
!   twiss2    -- Twiss_struct: Output Twiss parameters.
!   err       -- Logical: Set True if there is an error, false otherwise.
!-

subroutine twiss1_propagate (twiss1, mat2, ele_key, length, twiss2, err)

use bmad_interface, except_dummy => twiss1_propagate

implicit none

type (twiss_struct)  twiss1, twiss2, temp

real(rp) m11, m12, m21, m22, del_phi, length
real(rp) a1, b1, g1, a2, b2, g2, mat2(2,2), det
integer ele_key
logical err

!----------------------------------------------------
! Basic equation is given by Bovet 2.5.b page 16
! Linac rf matrices need to be renormalized.

err = .true.

det = determinant (mat2)

if (det == 0 .or. twiss1%beta == 0) return

m11 = mat2(1,1)
m12 = mat2(1,2)
m21 = mat2(2,1)
m22 = mat2(2,2)

a1 = twiss1%alpha
b1 = twiss1%beta
g1 = (1 + a1**2) / b1

b2 =       (m11**2  * b1 - 2*m11*m12 * a1 + m12**2  * g1) / det
a2 = a1 + (-m21*m11 * b1 + 2*m12*m21 * a1 - m12*m22 * g1) / det
g2 =  (1 + a2**2) /b2

del_phi = atan2(m12, m11*b1 - m12*a1)

if (ele_key /= patch$) then
  if (del_phi < 0 .and. length > 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. length < 0) del_phi = del_phi - twopi
endif

twiss2%beta = b2
twiss2%alpha = a2
twiss2%gamma = g2
twiss2%phi = twiss1%phi + del_phi

err = .false.

end subroutine
