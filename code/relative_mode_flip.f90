!+
! Function RELATIVE_MODE_FLIP (ELE1, ELE2)
!
! Function to see if the modes of ELE1 are flipped relative to ELE2.
! This is done by seeing which eigen planes are similar.
! This function is used with the same element of differing lattices
! to help identify which mode is which.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     ELE1, ELE2 -- Ele_struct: Elements to compare.
!
! Output:
!     RELATIVE_MODE_FLIP -- Logical: true if modes are relatively flipped.
!-


  function relative_mode_flip (ele1, ele2)

  use bmad_struct
  implicit none

  type (ele_struct)  ele1, ele2

  real mat4(4,4), conj_mat(2,2), det_aa, det_ab

  logical relative_mode_flip

! fill top half of 4x4 matrix with ele1 a-mode eigen axes

  mat4 = 0
  mat4(1,1) = ele1%gamma_c
  mat4(2,2) = ele1%gamma_c
  mat4(1:2,3:4) = ele1%c_mat

! fill bottom half of 4x4 matrix with ele2 a-mode eigen axes

  mat4(3,1) = ele2%gamma_c
  mat4(4,2) = ele2%gamma_c
  mat4(3:4,3:4) = ele2%c_mat

! smallness of determinate is indicator whether the a-modes of ele1 and
! ele2 have nearly the same eigen planes

  call mat_det (mat4, det_aa, 4, 4)

! fill bottom half of 4x4 matrix whith ele2 b-mode eigen axes

  mat4(3,3) = ele2%gamma_c
  mat4(4,4) = ele2%gamma_c
  mat4(3,4) = 0
  mat4(4,3) = 0
  call mat_symp_conj (ele2%c_mat, conj_mat, 2, 2)
  mat4(3:4,1:2) = -conj_mat

! smallness of determinate is indicator of whether ele1 a-mode has nearly the
! same eigen plane with ele2 b-mode

  call mat_det (mat4, det_ab, 4, 4)

! compare dets to see if relative mode flip

  if (abs(det_aa) < abs(det_ab)) then
    relative_mode_flip = .false.
  else
    relative_mode_flip = .true.
  endif

  return
  end
