!+
! Subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)
!
! Subroutine to form the 4x4 1-turn coupled matrix with the reference point
! at the end of an element.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ele   -- ele_struct: Reference element.
!     %x       -- "a" mode Twiss parameter structure.
!     %y       -- "b" mode Twiss parameter structure.
!     %c_mat   -- 2x2 C matrix.
!     %gamma_c -- gamma associated with C matrix.
!   phi_a -- Real: "a" mode tune in radians.
!   phi_b -- Real: "b" mode tune in radians.
!
! Output:
!   mat4(4,4) -- Real: 1-Turn coupled matrix.
!-

subroutine one_turn_mat_at_ele (ele, phi_a, phi_b, mat4)

  use bmad_struct
  use bmad_interface

  type (ele_struct) ele

  real phi_a, phi_b, mat4(4,4)
  real a(2,2), b(2,2), c_conj(2,2), c(2,2)
  real g

!

  call twiss_to_1_turn_mat (ele%x, phi_a, a)
  call twiss_to_1_turn_mat (ele%y, phi_b, b)
  call mat_symp_conj (ele%c_mat, c_conj, 2, 2)
  c = ele%c_mat
  g = ele%gamma_c

  mat4(1:2,1:2) = g**2 * a + matmul(matmul(c, b), c_conj)
  mat4(1:2,3:4) = g * (matmul(c, b) - matmul(a, c))
  mat4(3:4,1:2) = g * (matmul(b, c_conj) - matmul(c_conj, a))
  mat4(3:4,3:4) = g**2 * b + matmul(matmul(c_conj, a), c)

end subroutine
