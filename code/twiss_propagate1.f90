!+
! Subroutine twiss_propagate1 (ele1, ele2)
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
!   bmad_status -- Common block status structure:
!       %ok         -- Set False if an input beta is zero
!-

#include "CESR_platform.inc"

subroutine twiss_propagate1 (ele1, ele2)

  use bmad_struct
  use bmad_interface, except => twiss_propagate1

  implicit none
  type (ele_struct), target :: ele1, ele2
  type (ele_struct), save :: ele_temp

  real(rp), pointer :: mat6(:,:), orb(:), orb_out(:)
  real(rp) v_mat(4,4), v_inv_mat(4,4), y_inv(2,2), det
  real(rp) big_M(2,2), small_m(2,2), big_N(2,2), small_n(2,2)
  real(rp) c_conj_mat(2,2), E_inv_mat(2,2), F_inv_mat(2,2)
  real(rp) mat2(2,2), eta1_vec(6), eta_vec(6), dpz2_dpz1, rel_p1, rel_p2
  real(rp) mat4(4,4), det_factor, det_original

!---------------------------------------------------------------------
! init
! ELE_TEMP needs the element length for the betatron phase calculation

  if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then
    print *, 'ERROR IN TWISS_PROPAGATE1: ZERO BETA DETECTED.'
    if (bmad_status%exit_on_error) call err_exit
    bmad_status%ok = .false.
    return
  endif

  ele2%mode_flip = ele1%mode_flip          ! assume no flip
  ele_temp%value(l$) = ele2%value(l$)

!---------------------------------------------------------------------
! markers are easy

  if (ele2%key == marker$) then
    ele2%x = ele1%x
    ele2%y = ele1%y
    ele2%a = ele1%a
    ele2%b = ele1%b
    ele2%z = ele1%z
    ele2%gamma_c = ele1%gamma_c
    ele2%c_mat = ele1%c_mat
    return
  endif

!---------------------------------------------------------------------
! if transfer matrix is not coupled...
! propagate c_mat coupling matrix and setup temporary element for propagation


  if (all(ele2%mat6(1:2,3:4) == 0)) then

    ele_temp%mat6 = ele2%mat6

    call mat_symp_conj (ele2%mat6(3:4,3:4), y_inv) ! conj == inverse
    mat2 = matmul (ele2%mat6(1:2,1:2), ele1%c_mat)
    ele2%c_mat = matmul (mat2, y_inv)
    ele2%gamma_c = ele1%gamma_c

!---------------------------------------------------------------------
! here if we are dealing with a coupled transfer matrix

  else

    mat4 = ele2%mat6(1:4,1:4)

    ! det_factor is a renormalization factor since det_original != 1

    if (ele2%key == lcavity$) then  
      call mat_det (mat4, det_original)
      det_factor = sqrt(det_original)
    else
      det_factor = 1
    endif

    call mat_make_unit (ele_temp%mat6)   ! makle unit matrix

    big_M = mat4(1:2,1:2)
    small_m = mat4(1:2,3:4)
    big_N = mat4(3:4,3:4)
    small_n = mat4(3:4,1:2)

    call mat_symp_conj (ele1%c_mat, c_conj_mat)
    mat2 = ele1%gamma_c * big_M - matmul(small_m, c_conj_mat)
    call mat_det (mat2, det) 
    det = det / det_factor

! we demand that gamma_c > 0.3 (ie det > 0.1)
! try to make it so that there is no net mode flip here

    if (det > 0.9 .or. (det > 0.1 .and. .not. ele1%mode_flip)) then

      ele2%gamma_c = sqrt(det)
      ele_temp%mat6(1:2,1:2) = mat2 / ele2%gamma_c
      ele_temp%mat6(3:4,3:4) = &
            (ele1%gamma_c * big_N + matmul(small_n, ele1%c_mat)) / ele2%gamma_c
      call mat_symp_conj (ele_temp%mat6(3:4,3:4), F_inv_mat)
      ele2%c_mat = matmul(matmul(big_M, ele1%c_mat) + &
                                  ele1%gamma_c * small_m, F_inv_mat)

! else we flip the modes

    else

      mat2 = matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m
      call mat_det (mat2, det)
      det = det / det_factor
      if (det < 0) then
        print *, 'TWISS_PROPAGATE1: INTERNAL ERROR! (DUE TO ROUNDOFF?)'
      endif

      ele2%gamma_c = sqrt(abs(det))
      ele_temp%mat6(1:2,1:2) = &
            (ele1%gamma_c * small_n - matmul(big_N, c_conj_mat)) / ele2%gamma_c
      ele_temp%mat6(3:4,3:4) = mat2 / ele2%gamma_c

      call mat_symp_conj (ele_temp%mat6(1:2,1:2), E_inv_mat)
      ele2%c_mat = &
          matmul(ele1%gamma_c * big_M - matmul(small_m, c_conj_mat), E_inv_mat)

      ele2%mode_flip = .not. ele1%mode_flip

    endif

  endif

!---------------------------------------------------------------------
! Propagate twiss.

  ele_temp%mode_flip = ele2%mode_flip
  ele_temp%key       = ele2%key
  call twiss_decoupled_propagate (ele1, ele_temp)  
  ele2%a = ele_temp%a                     ! transfer twiss to ele2
  ele2%b = ele_temp%b
  ele2%x = ele_temp%x                     ! transfer twiss to ele2
  ele2%y = ele_temp%y

! Comming out of a flipped state correct phase by twopi.
! This is a kludge but factors of twopi are not physically meaningful.

  if (ele1%mode_flip .and. .not. ele2%mode_flip) then
    ele2%a%phi = ele2%a%phi - twopi
    ele2%b%phi = ele2%b%phi - twopi
  endif

!----------------------------------------------------
! Dispersion calc.
! p_z2 is p_z at end of ele2 assuming p_z = 1 at end of ele1.
! This is just 1.0 (except for RF cavities).

  mat6 => ele2%mat6
  orb => ele2%ref_orb_in
  orb_out => ele2%ref_orb_out
  rel_p1 = 1 + orb(6)               ! reference energy 
  rel_p2 = 1 + orb_out(6)

  eta1_vec = (/ ele1%x%eta, ele1%x%etap * rel_p1, &
                ele1%y%eta, ele1%y%etap * rel_p1, ele1%z%eta, 1.0_rp /)

  if (ele2%key == rfcavity$) then
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

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine twiss_decoupled_propagate (ele1, ele2)
!
! Subroutine to propagate the twiss parameters from end of ELE1 to end of ELE2.
! This routine assumes that the matrix ele2%mat6 is decoupled.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1  -- Ele_struct: Structure holding the starting Twiss parameters.
!   ele2  -- Ele_struct: Structure holding the transfer matrix.
!
! Output:
!   ele2  -- Ele_struct: Structure for the ending Twiss parameters.
!-

subroutine twiss_decoupled_propagate (ele1, ele2)

  use bmad_struct
  use bmad_interface, except => twiss_decoupled_propagate

  implicit none

  type (ele_struct)  ele1, ele2

  real(rp) m11, m12, m21, m22, a1, b1, g1, del_phi
  real(rp) a2, b2, g2, mat4(4,4), det, det_factor

!----------------------------------------------------
! Linac rf matrices need to be renormalized.

  mat4 = ele2%mat6(1:4, 1:4)

  if (ele2%key == lcavity$) then
    call mat_det (mat4, det)
    det_factor = det**0.25
  else
    det_factor = 1.0
  endif

!----------------------------------------------------
! Basic equation is given by Bovet 2.5.b page 16
! Propagate A mode ("X") of ele1

  m11 = mat4(1,1) / det_factor
  m12 = mat4(1,2) / det_factor
  m21 = mat4(2,1) / det_factor
  m22 = mat4(2,2) / det_factor

  a1 = ele1%a%alpha
  b1 = ele1%a%beta
  g1 = (1 + a1**2) / b1

  b2 =  m11**2  * b1 - 2*m11*m12     * a1 + m12**2  * g1
  a2 = -m21*m11 * b1 + (1+2*m12*m21) * a1 - m12*m22 * g1
  g2 =  (1 + a2**2) /b2

  del_phi = atan2(m12, m11*b1 - m12*a1)
  if (del_phi < 0 .and. ele2%value(l$) > 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. ele2%value(l$) < 0) del_phi = del_phi - twopi

  if (ele2%mode_flip .eqv. ele1%mode_flip) then
    ele2%a%beta = b2
    ele2%a%alpha = a2
    ele2%a%gamma = g2
    ele2%a%phi = ele1%a%phi + del_phi
  else
    ele2%b%beta = b2
    ele2%b%alpha = a2
    ele2%b%gamma = g2
    ele2%b%phi = ele1%a%phi + del_phi
  endif

!-----------------------------------------------------
! Propagate B mode ("Y") of ele1

  m11 = mat4(3,3) / det_factor
  m12 = mat4(3,4) / det_factor
  m21 = mat4(4,3) / det_factor
  m22 = mat4(4,4) / det_factor

  a1 = ele1%b%alpha
  b1 = ele1%b%beta
  g1 = (1 + a1**2) / b1

  b2 =  m11**2  * b1 - 2*m11*m12     * a1 + m12**2  * g1
  a2 = -m21*m11 * b1 + (1+2*m12*m21) * a1 - m12*m22 * g1
  g2 =  (1 + a2**2) /b2

  del_phi = atan2(m12, m11*b1 - m12*a1)
  if (del_phi < 0 .and. ele2%value(l$) > 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. ele2%value(l$) < 0) del_phi = del_phi - twopi

  if (ele2%mode_flip .eqv. ele1%mode_flip) then
    ele2%b%beta = b2
    ele2%b%alpha = a2
    ele2%b%gamma = g2
    ele2%b%phi = ele1%b%phi + del_phi
  else
    ele2%a%beta = b2
    ele2%a%alpha = a2
    ele2%a%gamma = g2
    ele2%a%phi = ele1%b%phi + del_phi
  endif

end subroutine
