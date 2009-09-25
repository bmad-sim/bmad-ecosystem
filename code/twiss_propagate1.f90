!+
! Subroutine twiss_propagate1 (ele1, ele2, err)
!
! Subroutine to propagate the twiss parameters from the end of ELE1 to
! the end of ELE2.
!
! Note: ele%a Twiss parameters are associated with the "A" mode and
! the ele%b Twiss parameters are associated with the "B" mode.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele1        -- Ele_struct: Structure holding the starting Twiss parameters.
!   ele2        -- Ele_struct: Structure holding the transfer matrix.
!   bmad_status -- Common block status structure:
!       %type_out -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   ele2   -- Ele_struct: Structure for the ending Twiss parameters.
!   err    -- Logical, optional: Set True if there is an error. False otherwise.
!   bmad_status -- Common block status structure:
!       %ok         -- Set False if an input beta is zero (deprecated).
!-

#include "CESR_platform.inc"

subroutine twiss_propagate1 (ele1, ele2, err)

  use bmad_struct
  use bmad_interface, except_dummy => twiss_propagate1

  implicit none
  type (ele_struct), target :: ele1, ele2
  type (twiss_struct) twiss_a
  type (lat_param_struct) param
  integer key2

  real(rp), pointer :: mat6(:,:), orb(:), orb_out(:)
  real(rp) v_mat(4,4), v_inv_mat(4,4), y_inv(2,2), det, mat2_a(2,2), mat2_b(2,2)
  real(rp) big_M(2,2), small_m(2,2), big_N(2,2), small_n(2,2)
  real(rp) c_conj_mat(2,2), E_inv_mat(2,2), F_inv_mat(2,2)
  real(rp) mat2(2,2), eta1_vec(6), eta_vec(6), dpz2_dpz1, rel_p1, rel_p2
  real(rp) mat4(4,4), det_factor

  logical error
  logical, optional :: err

  !---------------------------------------------------------------------
  ! init

  if (present(err)) err = .true.

  if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then
    print *, 'ERROR IN TWISS_PROPAGATE1: ZERO BETA DETECTED AT: ', trim(ele1%name)
    print *, '      ELEMENT #', ele1%ix_ele
    if (bmad_status%exit_on_error) call err_exit
    bmad_status%ok = .false.
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
    ele2%value(c_11$:c_22$) = &
                (/ ele1%c_mat(1,1), ele1%c_mat(1,2), ele1%c_mat(2,1), ele1%c_mat(2,2) /)
    ele2%value(gamma_c$) = ele1%gamma_c
    call make_mat6 (ele2, param)
  endif

  !---------------------------------------------------------------------
  ! markers are easy

  if (key2 == marker$ .or. key2 == photon_branch$ .or. key2 == branch$) then
    ele2%x = ele1%x
    ele2%y = ele1%y
    ele2%a = ele1%a
    ele2%b = ele1%b
    ele2%z = ele1%z
    ele2%gamma_c = ele1%gamma_c
    ele2%c_mat = ele1%c_mat
    if (present(err)) err = .false.
    return
  endif

  !---------------------------------------------------------------------
  ! if transfer matrix is not coupled...
  ! propagate c_mat coupling matrix and setup temporary element for propagation


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
    if (key2 == lcavity$) then
      det_factor = sqrt(determinant (mat4))
    endif

    mat4 = ele2%mat6(1:4,1:4)

    big_M = mat4(1:2,1:2)
    small_m = mat4(1:2,3:4)
    big_N = mat4(3:4,3:4)
    small_n = mat4(3:4,1:2)

    call mat_symp_conj (ele1%c_mat, c_conj_mat)
    mat2 = ele1%gamma_c * big_M - matmul(small_m, c_conj_mat)
    det = determinant(mat2) / det_factor 

    ! we demand that gamma_c > 0.3 (ie det > 0.1)
    ! try to make it so that there is no net mode flip here

    if (det > 0.9 .or. (det > 0.1 .and. .not. ele1%mode_flip)) then

      ele2%gamma_c = sqrt(det)
      mat2_a = mat2 / ele2%gamma_c
      mat2_b = (ele1%gamma_c * big_N + matmul(small_n, ele1%c_mat)) / ele2%gamma_c
      call mat_symp_conj (mat2_b, F_inv_mat)
      ele2%c_mat = matmul(matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m, F_inv_mat)

    ! else we flip the modes

    else

      mat2 = matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m
      det = determinant(mat2) / det_factor
      if (det < 0) then
        print *, 'TWISS_PROPAGATE1: INTERNAL ERROR! (DUE TO ROUNDOFF?)'
      endif

      ele2%gamma_c = sqrt(abs(det))
      mat2_a = (ele1%gamma_c * small_n - matmul(big_N, c_conj_mat)) / ele2%gamma_c
      mat2_b = mat2 / ele2%gamma_c

      call mat_symp_conj (mat2_a, E_inv_mat)
      ele2%c_mat = matmul(ele1%gamma_c * big_M - matmul(small_m, c_conj_mat), E_inv_mat)

      ele2%mode_flip = .not. ele1%mode_flip

    endif

  endif

  !---------------------------------------------------------------------
  ! Propagate twiss.

  call twiss1_propagate (ele1%a, mat2_a, ele2%value(l$), ele2%a, error)
  if (error) return
  call twiss1_propagate (ele1%b, mat2_b, ele2%value(l$), ele2%b, error)
  if (error) return

  if (ele2%mode_flip .neqv. ele1%mode_flip) then
    twiss_a = ele2%a
    ele2%a = ele2%b
    ele2%b = twiss_a
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

  mat6 => ele2%mat6
  orb  => ele2%map_ref_orb_in%vec
  orb_out => ele2%map_ref_orb_out%vec
  rel_p1 = 1 + orb(6)               ! reference energy 
  rel_p2 = 1 + orb_out(6)

  eta1_vec = (/ ele1%x%eta, ele1%x%etap * rel_p1, &
                ele1%y%eta, ele1%y%etap * rel_p1, ele1%z%eta, 1.0_rp /)

  if (key2 == rfcavity$) then
    eta1_vec(5) = 0
    dpz2_dpz1 = dot_product (mat6(6,1:4), eta1_vec(1:4)) + 1
  else
    dpz2_dpz1 = dot_product (mat6(6,:), eta1_vec) + &
                  (mat6(6,2) * orb(2) + mat6(6,4) * orb(4)) / rel_p1
  endif

  eta_vec(1:5) = matmul (ele2%mat6(1:5,:), eta1_vec) / dpz2_dpz1

  eta_vec(1) = eta_vec(1) + &
                (mat6(1,2) * orb(2) + mat6(1,4) * orb(4)) / (dpz2_dpz1 * rel_p1)
  eta_vec(2) = eta_vec(2) - orb_out(2) / rel_p2 + &
                (mat6(2,2) * orb(2) + mat6(2,4) * orb(4)) / (dpz2_dpz1 * rel_p1)
  eta_vec(3) = eta_vec(3) + &
                (mat6(3,2) * orb(2) + mat6(3,4) * orb(4)) / (dpz2_dpz1 * rel_p1)
  eta_vec(4) = eta_vec(4) - orb_out(4) / rel_p2 + &
                (mat6(4,2) * orb(2) + mat6(4,4) * orb(4)) / (dpz2_dpz1 * rel_p1)
  eta_vec(5) = eta_vec(5) + &
                (mat6(5,2) * orb(2) + mat6(5,4) * orb(4)) / (dpz2_dpz1 * rel_p1)

  eta_vec(2) = eta_vec(2) / rel_p2
  eta_vec(4) = eta_vec(4) / rel_p2

  ele2%x%eta  = eta_vec(1)
  ele2%x%etap = eta_vec(2)
  ele2%y%eta  = eta_vec(3)
  ele2%y%etap = eta_vec(4)
  ele2%z%eta  = eta_vec(5)

  call make_v_mats (ele2, v_mat, v_inv_mat)
  eta_vec(1:4) = matmul (v_inv_mat, eta_vec(1:4))

  ele2%a%eta  = eta_vec(1)
  ele2%a%etap = eta_vec(2)
  ele2%b%eta  = eta_vec(3)
  ele2%b%etap = eta_vec(4)

  if (present(err)) err = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine twiss1_propagate (twiss1, mat2, length, twiss2, err)
!
! Subroutine to propagate the twiss parameters of a single mode.
!
! The betatron phase phi is only determined up to a factor of 2pi. 
! The length argument is used to determine what phi should be:
!   length > 0:       0 <= phi <  twopi
!   length = 0:     -pi <  phi <= pi
!   length < 0:  -twopi <  phi <= 0
!
! Modules needed:
!   use bmad
!
! Input:
!   twiss1    -- Twiss_struct: Input Twiss parameters.
!   mat2(2,2) -- Real(rp): The transfer matrix.
!   length    -- Real(rp): Determines whether the phase is 
!                            increasing or decreasing.
!
! Output:
!   twiss2    -- Twiss_struct: Output Twiss parameters.
!   err       -- Logical: Set True if there is an error, false otherwise.
!-

subroutine twiss1_propagate (twiss1, mat2, length, twiss2, err)

  use bmad_struct
  use bmad_interface, except_dummy => twiss1_propagate

  implicit none

  type (twiss_struct)  twiss1, twiss2, temp

  real(rp) m11, m12, m21, m22, del_phi, length
  real(rp) a1, b1, g1, a2, b2, g2, mat2(2,2), det

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
  if (del_phi < 0 .and. length > 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. length < 0) del_phi = del_phi - twopi

  twiss2%beta = b2
  twiss2%alpha = a2
  twiss2%gamma = g2
  twiss2%phi = twiss1%phi + del_phi

  err = .false.

end subroutine
