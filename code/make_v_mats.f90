!+
! Subroutine make_v_mats (ele, v_mat, v_inv_mat)
!
! Subroutine make the matrices needed to go from normal mode coords
! to X-Y coords and vice versa.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- Ele_struct: Element
!
! Output:
!   v_mat(4,4)     -- Real(rp): Normal mode to X-Y coords transformation
!   v_inv_mat(4,4) -- Real(rp): X-Y coords to Normal mode transformation
!-

#include "CESR_platform.inc"

subroutine make_v_mats (ele, v_mat, v_inv_mat)

  use bmad_struct
  use bmad_interface, except => make_v_mats
  
  implicit none

  type (ele_struct)  ele

  real(rp) v_mat(4,4), v_inv_mat(4,4), c_conj(2,2)
  integer i

!

  v_mat = 0
  v_inv_mat = 0

  do i = 1, 4
    v_mat(i,i)     = ele%gamma_c
    v_inv_mat(i,i) = ele%gamma_c
  enddo

  call mat_symp_conj (ele%c_mat, c_conj)

  v_mat(1:2,3:4) = ele%c_mat
  v_mat(3:4,1:2) = -c_conj  

  v_inv_mat(1:2,3:4) = -ele%c_mat
  v_inv_mat(3:4,1:2) = c_conj

end subroutine
