!+
! Subroutine INIT_ELE (ELE)
!
! Subroutine to initialize a BMAD element. Element is initialized to be free
! (not a lord or slave) and all %VALUES set to zero.
! 
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Output:
!   ELE -- Ele_struct: Initialized element.
!-

subroutine init_ele (ele)

  use bmad_struct

  implicit none

  type (ele_struct)  ele

!

  ele%type = ' '
  ele%alias = ' '
  ele%name = '****'

  ele%value(:) = 0

  ele%ix_value = 0
  ele%ic1_lord = 0
  ele%ic2_lord = -1
  ele%n_lord = 0
  ele%ix1_slave = 0
  ele%ix2_slave = -1
  ele%n_slave = 0
  ele%control_type = free$
  ele%ix_pointer = 0
  ele%s = 0
  ele%is_on = .true.
  ele%multipoles_on = .true.
  ele%nonzero_multipoles = .false.

end subroutine
