!+
! Subroutine MAKE_V_MATS (ELE, V_MAT, V_INV_MAT)
!
! Subroutine make the matrices needed to go from normal mode coords
! to X-Y coords and vice versa.
!
! Modules needed:
!   use bmad
!
! Input:
!   ELE        -- Ele_struct: Element
!
! Output:
!   V_MAT(4,4)     -- Real(rp): Normal mode to X-Y coords transformation
!   V_INV_MAT(4,4) -- Real(rp): X-Y coords to Normal mode transformation
!-

!$Id$
!$Log$
!Revision 1.6  2003/07/09 01:38:16  dcs
!new bmad with allocatable ring%ele_(:)
!
!Revision 1.5  2002/10/29 17:07:13  dcs
!*** empty log message ***
!
!Revision 1.4  2002/02/23 20:32:18  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:40  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:53  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine make_v_mats (ele, v_mat, v_inv_mat)

  use bmad_struct
  use bmad_interface
  
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

  call mat_symp_conj (ele%c_mat, c_conj, 2, 2)

  v_mat(1:2,3:4) = ele%c_mat
  v_mat(3:4,1:2) = -c_conj  

  v_inv_mat(1:2,3:4) = -ele%c_mat
  v_inv_mat(3:4,1:2) = c_conj

end subroutine
