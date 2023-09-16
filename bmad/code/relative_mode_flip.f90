!+
! Function relative_mode_flip (ele1, ele2)
!
! Function to see if the modes of ELE1 are flipped relative to ELE2.
! This is done by seeing which eigen planes are similar.
! This function is used with the same element of different lattices
! to help identify which mode is which.
!
! Input:
!     ele1, ele2 -- Ele_struct: Elements to compare.
!
! Output:
!     relative_mode_flip -- Logical: true if modes are relatively flipped.
!-

function relative_mode_flip (ele1, ele2) result (rel_mode)

use bmad_interface, except_dummy => relative_mode_flip

implicit none

type (ele_struct)  ele1, ele2

real(rp) mat4(4,4), det_aa, det_ab

logical rel_mode

! fill top half of 4x4 matrix with ele1 a-mode eigen axes

mat4 = 0
mat4(1,1) = ele1%gamma_c
mat4(2,2) = ele1%gamma_c
mat4(1:2,3:4) = ele1%c_mat

! fill bottom half of 4x4 matrix with ele2 a-mode eigen axes

mat4(3,1) = ele2%gamma_c
mat4(4,2) = ele2%gamma_c
mat4(3:4,3:4) = ele2%c_mat

! smallness of determinant is indicator whether the a-modes of ele1 and
! ele2 have nearly the same eigen planes

det_aa = determinant (mat4)

! fill bottom half of 4x4 matrix whith ele2 b-mode eigen axes

mat4(3,3) = ele2%gamma_c
mat4(4,4) = ele2%gamma_c
mat4(3,4) = 0
mat4(4,3) = 0
mat4(3:4,1:2) = -mat_symp_conj(ele2%c_mat)

! smallness of determinant is indicator of whether ele1 a-mode has nearly the
! same eigen plane with ele2 b-mode

det_ab = determinant (mat4)

! compare dets to see if relative mode flip

if (abs(det_aa) < abs(det_ab)) then
  rel_mode = .false.
else
  rel_mode = .true.
endif


end function
