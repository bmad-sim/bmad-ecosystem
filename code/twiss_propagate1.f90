!+
! Subroutine twiss_propagate1 (ele1, ele2)
!
! Subroutine to propagate the twiss parameters from the end of ELE1 to
! the end of ELE2.
!
! Note: ele%x Twiss parameters are associated with the "A" mode and
! the ele%y Twiss parameters are associated with the "B" mode.
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
  use bmad_interface

  implicit none
  type (ele_struct)  ele1, ele2, ele_temp

  integer i, j

  real(rp) v_mat(4,4), v_inv_mat(4,4), amat2(2,2), y_inv(2,2), det
  real(rp) big_M(2,2), small_m(2,2), big_N(2,2), small_n(2,2)
  real(rp) c_conj_mat(2,2), E_inv_mat(2,2), F_inv_mat(2,2)
  real(rp) mat2(2,2), vec(4)
  real(rp) eta_x(4), eta_a(4)

!---------------------------------------------------------------------
! init
! ELE_TEMP needs the element length for the betatron phase calculation

  if (ele1%x%beta == 0 .or. ele1%y%beta == 0) then
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
    ele2%gamma_c = ele1%gamma_c
    ele2%c_mat = ele1%c_mat
    return
  endif

!---------------------------------------------------------------------
! if transfer matrix is not coupled...
! propagate c_mat coupling matrix and setup temporary element for propagation


  if (all(ele2%mat6(1:2,3:4) == 0)) then

    call mat_symp_conj (ele2%mat6(3:4,3:4), y_inv) ! conj == inverse
    mat2 = matmul (ele2%mat6(1:2,1:2), ele1%c_mat)
    ele2%c_mat = matmul (mat2, y_inv)
    ele2%gamma_c = ele1%gamma_c
    ele_temp%mat6(1:4,1:4) = ele2%mat6(1:4,1:4)

!---------------------------------------------------------------------
! here if we are dealing with a coupled transfer matrix

  else

    call mat_make_unit (ele_temp%mat6)   ! makle unit matrix

    big_M = ele2%mat6(1:2,1:2)
    small_m = ele2%mat6(1:2,3:4)
    big_N = ele2%mat6(3:4,3:4)
    small_n = ele2%mat6(3:4,1:2)

    call mat_symp_conj (ele1%c_mat, c_conj_mat)
    mat2 = ele1%gamma_c * big_M - matmul(small_m, c_conj_mat)
    call mat_det (mat2, det)

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

!!      if (bmad_status%type_out) print *, 'TWISS_PROPAGATE1: MODE_FLIPPED'

      mat2 = matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m
      call mat_det (mat2, det)
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

!----------------------------------------------------
! linac rf matrices need to be renormalized.

  if (ele2%key == lcavity$) then
    call mat_det (ele_temp%mat6(1:4,1:4), det)
    ele_temp%mat6(1:4,1:4) = ele_temp%mat6(1:4,1:4) / (det**0.25)
  endif

!---------------------------------------------------------------------
! calculate components in transfer matrix for dispersion calc
! effective_mat6(i,6) = V_2^-1(i,j) * ele2.mat6(j,6)

  call make_v_mats (ele2, v_mat, v_inv_mat)
  ele_temp%mat6(1:4,6) = matmul(v_inv_mat(1:4,1:4), ele2%mat6(1:4,6))

! propagate twiss

  ele_temp%mode_flip = ele2%mode_flip
  call twiss_decoupled_propagate (ele1, ele_temp)   ! now calc new twiss
  ele2%x = ele_temp%x                     ! transfer twiss to ele2
  ele2%y = ele_temp%y

! comming out of a flipped state correct phase by twopi
! this is a kludge but factors of twopi are not physically meaningful

  if (ele1%mode_flip .and. .not. ele2%mode_flip) then
    ele2%x%phi = ele2%x%phi - twopi
    ele2%y%phi = ele2%y%phi - twopi
  endif

! calc laboratory dispersions

  call make_v_mats (ele2, v_mat, v_inv_mat)
  eta_a = (/ ele2%x%eta, ele2%x%etap, &
                    ele2%y%eta, ele2%y%etap /)
  eta_x = matmul (v_mat, eta_a)
  ele2%x%eta_lab  = eta_x(1)
  ele2%x%etap_lab = eta_x(2)
  ele2%y%eta_lab  = eta_x(3)
  ele2%y%etap_lab = eta_x(4)


end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine TWISS_DECOUPLED_PROPAGATE (ele1, ele2)
!
! Subroutine to propagate the twiss parameters from end of ELE1 to end of ELE2.
!-


subroutine twiss_decoupled_propagate (ele1, ele2)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct)  ele1, ele2

  real(rp) m11, m12, m21, m22, a1, b1, g1, del_phi
  real(rp) a2, b2, g2

! Basic equation is given by Bovet 2.5.b page 16
! Propagate A mode ("X") of ele1

  m11 = ele2%mat6(1,1)
  m12 = ele2%mat6(1,2)
  m21 = ele2%mat6(2,1)
  m22 = ele2%mat6(2,2)

  a1 = ele1%x%alpha
  b1 = ele1%x%beta
  g1 = (1 + a1**2) / b1

  b2 =  m11**2  * b1 - 2*m11*m12     * a1 + m12**2  * g1
  a2 = -m21*m11 * b1 + (1+2*m12*m21) * a1 - m12*m22 * g1
  g2 =  (1 + a2**2) /b2

  del_phi = atan2(m12, m11*b1 - m12*a1)
  if (del_phi < 0 .and. ele2%value(l$) >= 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. ele2%value(l$) < 0) del_phi = del_phi - twopi

  if (ele2%mode_flip .eqv. ele1%mode_flip) then
    ele2%x%beta = b2
    ele2%x%alpha = a2
    ele2%x%gamma = g2
    ele2%x%eta  = m11*ele1%x%eta + m12*ele1%x%etap + ele2%mat6(1,6)
    ele2%x%etap = m21*ele1%x%eta + m22*ele1%x%etap + ele2%mat6(2,6)
    ele2%x%phi = ele1%x%phi + del_phi
  else
    ele2%y%beta = b2
    ele2%y%alpha = a2
    ele2%y%gamma = g2
    ele2%y%eta  = m11*ele1%x%eta + m12*ele1%x%etap + ele2%mat6(3,6)
    ele2%y%etap = m21*ele1%x%eta + m22*ele1%x%etap + ele2%mat6(4,6)
    ele2%y%phi = ele1%x%phi + del_phi
  endif

!-----------------------------------------------------
! Propagate B mode ("Y") of ele1

  m11 = ele2%mat6(3,3)
  m12 = ele2%mat6(3,4)
  m21 = ele2%mat6(4,3)
  m22 = ele2%mat6(4,4)

  a1 = ele1%y%alpha
  b1 = ele1%y%beta
  g1 = (1 + a1**2) / b1

  b2 =  m11**2  * b1 - 2*m11*m12     * a1 + m12**2  * g1
  a2 = -m21*m11 * b1 + (1+2*m12*m21) * a1 - m12*m22 * g1
  g2 =  (1 + a2**2) /b2

  del_phi = atan2(m12, m11*b1 - m12*a1)
  if (del_phi < 0 .and. ele2%value(l$) >= 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. ele2%value(l$) < 0) del_phi = del_phi - twopi

  if (ele2%mode_flip .eqv. ele1%mode_flip) then
    ele2%y%beta = b2
    ele2%y%alpha = a2
    ele2%y%gamma = g2
    ele2%y%eta  = m11*ele1%y%eta + m12*ele1%y%etap + ele2%mat6(3,6)
    ele2%y%etap = m21*ele1%y%eta + m22*ele1%y%etap + ele2%mat6(4,6)
    ele2%y%phi = ele1%y%phi + del_phi
  else
    ele2%x%beta = b2
    ele2%x%alpha = a2
    ele2%x%gamma = g2
    ele2%x%eta  = m11*ele1%y%eta + m12*ele1%y%etap + ele2%mat6(1,6)
    ele2%x%etap = m21*ele1%y%eta + m22*ele1%y%etap + ele2%mat6(2,6)
    ele2%x%phi = ele1%y%phi + del_phi
  endif

end subroutine
