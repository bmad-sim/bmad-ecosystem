!+
! Subroutine MAKE_V_MATS (ELE, V_MAT, V_INV_MAT)
!
! Subroutine make the matrices needed to go from normal mode coords
! to X-Y coords and vice versa.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ELE        -- Ele_struct: Element
!
! Output:
!   V_MAT(4,4)     -- Real: Normal mode to X-Y coords transformation
!   V_INV_MAT(4,4) -- Real: X-Y coords to Normal mode transformation
!-


subroutine make_v_mats (ele, v_mat, v_inv_mat)

  use bmad_struct
  implicit none

  type (ele_struct)  ele

  real v_mat(4,4), v_inv_mat(4,4), c_conj(2,2)
  integer i

!

  v_mat = 0
  v_inv_mat = 0

  do i = 1, 4
    v_mat(i,i)     = ele%gamma_c
    v_inv_mat(i,i) = ele%gamma_c
  enddo

  call mat_symp_conj (ele%c_mat, c_conj, 2, 2)

  v_mat(1:2,3:4) = ele%c_mat
  v_mat(3:4,1:2) = -c_conj  

  v_inv_mat(1:2,3:4) = -ele%c_mat
  v_inv_mat(3:4,1:2) = c_conj

end subroutine
