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

!$Id$
!$Log$
!Revision 1.3  2002/01/08 21:44:39  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine init_ele (ele)

  use bmad_struct
  use bmad_interface

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
  ele%mat6_calc_method = bmad_standard$
  ele%tracking_method = bmad_standard$
  ele%num_steps = 0

end subroutine
