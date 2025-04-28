!+
! Subroutine transfer_mat_from_twiss (ele1, ele2, orb1, orb2, m)
!
! Subroutine to make a 6 x 6 transfer matrix from the twiss parameters
! at two points.
!
! Note: This routine assumes no RF acceleration and the value of m(5,6) will be zero. 
!
! Input:
!   ele1     -- Ele_struct: Element with twiss parameters for the starting point.
!     %a, %b       -- a-mode and b-mode Twiss paramters
!       %beta         -- Beta parameter.
!       %alpha        -- Alpha parameter.
!       %phi          -- Phase at initial point.
!     %x, %y       -- dispersion values
!       %eta          -- Dispersion at initial point.
!       %etap         -- Dispersion derivative at initial point.
!     %c_mat(2,2)  -- Coupling matrix
!   ele2     -- Ele_struct: Element with twiss parameters for the ending point.
!   orb1(6)  -- real(rp): Reference orbit at ele1 (affects m(i,6) dispersion terms).
!   orb2(6)  -- real(rp): Reference orbit at ele2 (affects m(i,6) dispersion terms).
!
! Output:
!   m(6,6) -- Real(rp): Transfer matrix between the two points.
!-

subroutine transfer_mat_from_twiss (ele1, ele2, orb1, orb2, m)

use bmad_routine_interface, dummy => transfer_mat_from_twiss

implicit none

type (ele_struct) ele1, ele2

real(rp) orb1(6), orb2(6)
real(rp) m(6,6), v_mat(4,4), v_inv_mat(4,4), det
real(rp) rel_p1, rel_p2, eta1(4), eta2(4)
character(*), parameter :: r_name = 'transfer_mat_from_twiss'

! Error check

if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then
  call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT1: ' // ele1%name)
  if (global_com%exit_on_error) call err_exit
endif

if (ele2%a%beta == 0 .or. ele2%b%beta == 0) then
  call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT2: ' // ele2%name)
  if (global_com%exit_on_error) call err_exit
endif

! Transfer matrices without coupling or dispersion

m = 0
m(5,5) = 1
m(6,6) = 1

if (ele1%mode_flip .and. ele2%mode_flip) then
  call transfer_mat2_from_twiss (ele1%b, ele2%b, m(1:2,1:2))
  call transfer_mat2_from_twiss (ele1%a, ele2%a, m(3:4,3:4))
elseif (ele1%mode_flip .and. .not. ele2%mode_flip) then
  call transfer_mat2_from_twiss (ele1%b, ele2%b, m(3:4,1:2))  ! E sub-matrix
  call transfer_mat2_from_twiss (ele1%a, ele2%a, m(1:2,3:4))  ! F sub-matrix
elseif (.not. ele1%mode_flip .and. ele2%mode_flip) then
  call transfer_mat2_from_twiss (ele1%a, ele2%a, m(3:4,1:2))  ! E sub-matrix
  call transfer_mat2_from_twiss (ele1%b, ele2%b, m(1:2,3:4))  ! F sub-matrix
elseif (.not. ele1%mode_flip .and. .not. ele2%mode_flip) then
  call transfer_mat2_from_twiss (ele1%a, ele2%a, m(1:2,1:2))
  call transfer_mat2_from_twiss (ele1%b, ele2%b, m(3:4,3:4))
endif

! Add in coupling

if (any(ele1%c_mat /= 0)) then
  det = determinant (ele1%c_mat)
  ele1%gamma_c = sqrt(1-det)
  call make_v_mats (ele1, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (m(1:4,1:4), v_inv_mat)
endif

if (any(ele2%c_mat /= 0)) then
  det = determinant (ele2%c_mat)
  ele2%gamma_c = sqrt(1-det)
  call make_v_mats (ele2, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (v_mat, m(1:4,1:4))
endif

! Add in dispersion.
! See the Bmad manual for a derivation of the equations here.

rel_p1 = 1 + orb1(6)
rel_p2 = 1 + orb2(6)

eta1 = [ele1%x%eta, rel_p1 * ele1%x%etap, ele1%y%eta, rel_p1 * ele1%y%etap]
eta2 = [ele2%x%eta, rel_p2 * ele2%x%etap, ele2%y%eta, rel_p2 * ele2%y%etap]
eta2 = eta2 - (m(1:4,2) * orb1(2) + m(1:4,4) * orb1(4)) / rel_p1 + [0.0_rp, orb2(2), 0.0_rp, orb2(4)] / rel_p2

m(1:4,6) = eta2 - matmul (m(1:4,1:4), eta1) 

! The m(5,x) terms follow from the symplectic condition.

m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)

end subroutine transfer_mat_from_twiss

