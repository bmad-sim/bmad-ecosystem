!+
! Subroutine BETA_AVE(ring, ix_ele, betaxAve, betayAve)
!
! Calculate the vertical and horizontal average betas for a ring element
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     ring -- record/ring_struct/:  ring
!     ix_ele -- integer: ring location of element
!
! Output:
!     betaxAve -- real:	 average betax for the element
!     betayAve -- real:  average betay for the element
!-

!$Id$
!$Log$
!Revision 1.3  2002/01/08 21:44:36  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine Beta_Ave(ring, ix_ele, betaxAve, betayAve)

 	use bmad_struct
 	use bmad_interface

  implicit none

 	type (ring_struct) ring
 	type (ele_struct) ele
 	real betaxAve, betayAve, beave, kx, ky, length
 	integer ix_ele, con_index, c_index
  	
 	ele = ring%ele_(ix_ele)
  	
 	if (ele%key==rbend$ .or. ele%key==sbend$) then
 	  kx = 2.0/ele%value(l$)**2 * (cos(ele%value(angle$)) - 1.0)
 	  ky = 0
 	elseif(ele%key==quadrupole$ .or. ele%key==sol_quad$) then
 	  kx = ele%value(k1$)
 	  ky = -ele%value(k1$)
 	else
     kx = 0
 	  ky = 0
 	endif

  length = ele%value(l$)

! The following code is necessary due to limitations of the bmad controller
! implementation.  It works for now, but will probably break in the future.

  if (ele%n_slave > 0) then
    length = 0.0
    do con_index = ele%ix1_slave, ele%ix2_slave
      c_index = ring%control_(con_index)%ix_slave
      length = length + ring%ele_(c_index)%value(l$)
    enddo
  endif
  	
! Since the beta stored in the ELE is at the end, we need
! to invert alpha for BEAVE routine

  	betaxAve = beave(ele%x%beta, ele%x%alpha, kx, -length)
  	betayAve = beave(ele%y%beta, ele%y%alpha, ky, -length)
  	
end subroutine
