!+
! Subroutine TWISS_PROPAGATE1 (ELE1, ELE2)
!
! Subroutine to propagate the twiss parameters from the end of ELE1 to
! the end of ELE2.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ELE1        -- Ele_struct: Structure holding the starting Twiss parameters.
!   ELE2        -- Ele_struct: Structure holding the transfer matrix.
!   BMAD_STATUS -- Common block status structure
!       %TYPE_OUT -- Logical: If .true. then will type a message if the
!                    modes are flipped.
!
! Output:
!   ELE2   -- Ele_struct: Structure for the ending Twiss parameters.
!
! Note:
!
! Note: ELE%X Twiss parameters are associated with the "A" mode and
! the ELE%Y Twiss parameters are associated with the "B" mode.
!-

!$Id$
!$Log$
!Revision 1.4  2002/06/13 14:54:31  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:29  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine twiss_propagate1 (ele1, ele2)

  use bmad
  implicit none
  type (ele_struct)  ele1, ele2, ele_temp

  integer i, j

  real(rdef) v_mat(4,4), v_inv_mat(4,4), amat2(2,2), y_inv(2,2), det
  real(rdef) big_M(2,2), small_m(2,2), big_N(2,2), small_n(2,2)
  real(rdef) c_conj_mat(2,2), E_inv_mat(2,2), F_inv_mat(2,2)
  real(rdef) mat2(2,2), vec(4)

! init
! ELE_TEMP needs the element length for the betatron phase calculation

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

!! custom elements are even easier
!!  if (ele2%key == custom$) then
!!    call custom_twiss_propagate1 (ele1, ele2)
!!    return
!!  endif

!---------------------------------------------------------------------
! if transfer matrix is not coupled...
! propagate c_mat coupling matrix and setup temporary element for propagation


  if (all(ele2%mat6(1:2,3:4) == 0)) then

    call mat_symp_conj (ele2%mat6(3:4,3:4), y_inv, 2, 2) ! conj == inverse
    mat2 = matmul (ele2%mat6(1:2,1:2), ele1%c_mat)
    ele2%c_mat = matmul (mat2, y_inv)
    ele2%gamma_c = ele1%gamma_c
    ele_temp%mat6(1:4,1:4) = ele2%mat6(1:4,1:4)

!---------------------------------------------------------------------
! here if we are dealing with a coupled transfer matrix

  else

    call mat_unit (ele_temp%mat6, 6, 6)   ! makle unit matrix

    big_M = ele2%mat6(1:2,1:2)
    small_m = ele2%mat6(1:2,3:4)
    big_N = ele2%mat6(3:4,3:4)
    small_n = ele2%mat6(3:4,1:2)

    call mat_symp_conj (ele1%c_mat, c_conj_mat, 2, 2)
    mat2 = ele1%gamma_c * big_M - matmul(small_m, c_conj_mat)
    call mat_det (mat2, det, 2, 2)

! we demand that gamma_c > 0.3 (ie det > 0.1)
! try to make it so that there is no net mode flip here

    if (det > 0.9 .or. (det > 0.1 .and. .not. ele1%mode_flip)) then

      ele2%gamma_c = sqrt(det)
      ele_temp%mat6(1:2,1:2) = mat2 / ele2%gamma_c
      ele_temp%mat6(3:4,3:4) = &
            (ele1%gamma_c * big_N + matmul(small_n, ele1%c_mat)) / ele2%gamma_c
      call mat_symp_conj (ele_temp%mat6(3:4,3:4), F_inv_mat, 2, 2)
      ele2%c_mat = matmul(matmul(big_M, ele1%c_mat) + &
                                  ele1%gamma_c * small_m, F_inv_mat)

! else we flip the modes

    else

      if (bmad_status%type_out) type *, 'TWISS_PROPAGATE1: MODE_FLIPPED'

      mat2 = matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m
      call mat_det (mat2, det, 2, 2)
      if (det < 0) then
        type *, 'TWISS_PROPAGATE1: INTERNAL ERROR! (DUE TO ROUNDOFF?)'
      endif
      ele2%gamma_c = sqrt(abs(det))
      ele_temp%mat6(1:2,1:2) = &
            (ele1%gamma_c * small_n - matmul(big_N, c_conj_mat)) / ele2%gamma_c
      ele_temp%mat6(3:4,3:4) = mat2 / ele2%gamma_c

      call mat_symp_conj (ele_temp%mat6(1:2,1:2), E_inv_mat, 2, 2)
      ele2%c_mat = &
          matmul(ele1%gamma_c * big_M - matmul(small_m, c_conj_mat), E_inv_mat)

      ele2%mode_flip = .not. ele1%mode_flip

    endif

  endif

!---------------------------------------------------------------------
! calculate components in transfer matrix for dispersion calc
! effective_mat6(i,6) = V_2^-1(i,j) * ele2.mat6(j,6)

  call make_v_mats (ele2, v_mat, v_inv_mat)
  ele_temp%mat6(1:4,6) = matmul(v_inv_mat(1:4,1:4), ele2%mat6(1:4,6))

! propagate twiss

  call twiss_decoupled_propagate (ele1, ele_temp)   ! now calc new twiss
  ele2%x = ele_temp%x                     ! transfer twiss to ele2
  ele2%y = ele_temp%y

! comming out of a flipped state correct phase by twopi
! this is a kludge but factors of twopi are not physically meaningful

  if (ele1%mode_flip .and. .not. ele2%mode_flip) then
    ele2%x%phi = ele2%x%phi - twopi
    ele2%y%phi = ele2%y%phi - twopi
  endif

! calc MOBIUS_BETA and MOBIUS_ETA

  call mobius_twiss_calc (ele2, v_mat)

  return
  end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine TWISS_DECOUPLED_PROPAGATE (ele1, ele2)
!
! Subroutine to propagate the twiss parameters from end of ELE1 to end of ELE2.
!-


subroutine twiss_decoupled_propagate (ele1, ele2)

  use bmad
  implicit none

  type (ele_struct)  ele1, ele2

  real(rdef) m11, m12, m21, m22, a1, b1, g1, del_phi
  real(rdef) a2, b2, g2

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

  if (ele2%mode_flip == ele1%mode_flip) then
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

  if (ele2%mode_flip == ele1%mode_flip) then
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
