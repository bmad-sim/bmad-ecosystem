!+
! Subroutine INIT_ELE (ELE)
!
! Subroutine to initialize a BMAD element. Element is initialized to be free
! (not a lord or slave) and all %VALUES set to zero.
! 
! Modules needed:
!   use bmad
!
! Output:
!   ELE -- Ele_struct: Initialized element.
!-

!$Id$
!$Log$
!Revision 1.7  2002/12/04 07:46:23  dcs
!Init Twiss also
!
!Revision 1.6  2002/11/26 05:19:32  dcs
!Modified for BEGINNING floor position entry.
!
!Revision 1.5  2002/06/13 14:54:26  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:16  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:39  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine init_ele (ele)

  use bmad

  implicit none

  type (ele_struct)  ele

!

  ele%type = ' '
  ele%alias = ' '
  ele%name = '****'

  ele%key = 0
  ele%sub_key = 0

  ele%value(:) = 0


  ele%control_type = free$
  ele%ix_value = 0
  ele%ic1_lord = 0
  ele%ic2_lord = -1
  ele%n_lord = 0
  ele%ix1_slave = 0
  ele%ix2_slave = -1
  ele%n_slave = 0
  ele%ix_pointer = 0
  ele%s = 0

  ele%x_position = 0
  ele%y_position = 0
  ele%z_position = 0
  ele%theta_position = 0
  ele%phi_position = 0

  ele%mat6_calc_method = bmad_standard$
  ele%tracking_method = bmad_standard$
  ele%num_steps = 1
  ele%integration_order = 2
  ele%ptc_kind = 0

  ele%is_on = .true.
  ele%multipoles_on = .true.
  ele%symplectify = .false.
  ele%exact_rad_int_calc = .false.

  call deallocate_ele_pointers (ele)

! init Twiss

  ele%c_mat = 0
  ele%gamma_c = 1.0

  ele%x%beta  = 0
  ele%x%alpha = 0
  ele%x%gamma = 0
  ele%x%eta   = 0
  ele%x%etap  = 0
  ele%x%phi   = 0
  ele%x%sigma = 0

  ele%y%beta  = 0
  ele%y%alpha = 0
  ele%y%gamma = 0
  ele%y%eta   = 0
  ele%y%etap  = 0
  ele%y%phi   = 0
  ele%y%sigma = 0

end subroutine
