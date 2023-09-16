!+
! Subroutine make_v_mats (ele, v_mat, v_inv_mat)
!
! Subroutine make the matrices needed to go from normal mode coords
! to X-Y coords and vice versa.
!
! Input:
!   ele        -- Ele_struct: Element
!
! Output:
!   v_mat(4,4)     -- Real(rp), optional: Normal mode to X-Y coords transformation
!   v_inv_mat(4,4) -- Real(rp), optional: X-Y coords to Normal mode transformation
!-

subroutine make_v_mats (ele, v_mat, v_inv_mat)

use bmad_interface, except_dummy => make_v_mats

implicit none

type (ele_struct)  ele

real(rp), optional :: v_mat(4,4), v_inv_mat(4,4)
real(rp) c_conj(2,2)
integer i

! Init

c_conj = mat_symp_conj (ele%c_mat)

! construct v_max

if (present(v_mat)) then
  v_mat = 0
  forall (i = 1:4) v_mat(i,i) = ele%gamma_c
  v_mat(1:2,3:4) = ele%c_mat
  v_mat(3:4,1:2) = -c_conj  
endif

! Construct v_inv_mat

if (present(v_inv_mat)) then
  v_inv_mat = 0
  forall (i = 1:4) v_inv_mat(i,i) = ele%gamma_c
  v_inv_mat(1:2,3:4) = -ele%c_mat
  v_inv_mat(3:4,1:2) = c_conj
endif

end subroutine
